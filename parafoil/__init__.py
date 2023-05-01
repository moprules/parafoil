from . import builder
from . import matrix
from . import planet
from . import mass
from . import mesh
from . import aerodynamic
from . import sim
from . import autopilot
from . import store
import yaml
import os
import numpy as np


class PFSim:
    """Класс - симулятор парафойла"""

    def __init__(self, path_to_model: str = "", path_to_res: str = ""):
        self.params: dict = {}
        self.state: dict = {}
        self.model: dict = {}
        self.files: dict = {}
        self.path_to_model = os.path.abspath(path_to_model)
        self.path_to_res = os.path.abspath(path_to_res)

        if path_to_model:
            self.load(path_to_model)

    def load(self, path_to_model: str = ""):
        """Загрузка параметров из исходного файла"""
        if path_to_model:
            self.path_to_model = os.path.abspath(path_to_model)
        with open(self.path_to_model, "r") as stream:
            self.params = yaml.safe_load(stream)

    def build(self):
        self.model = builder.build(self.params)

    def init_files(self):
        store.init_files(self.files, self.path_to_res, self.path_to_model)
        store.init_coords(self.files, self.model)

    def init_state(self):
        # Новый словарь состояний
        self.state = {}
        self.state["time"] = self.model["time"]["init"]
        self.state["roll"] = self.model["body"]["roll"]
        self.state["pitch"] = self.model["body"]["pitch"]
        self.state["yaw"] = self.model["body"]["yaw"]
        self.state["rigging_angle"] = self.model["canopy"]["rigging_angle"]
        self.state["pos_east"] = self.model["body"]["east"]
        self.state["pos_north"] = self.model["body"]["north"]
        self.state["altitude"] = self.model["body"]["altitude"]

        # Снижение высоты корабля
        self.state["lost_altitude"] = 0
        # Пройденное расстояние
        self.state["traveled_distance"] = 0
        # Значение параметров на предыдущем шаге
        self.state["prev"] = {}
        self.state["prev"]["roll"] = self.state["roll"]
        self.state["prev"]["pitch"] = self.state["pitch"]
        self.state["prev"]["yaw"] = self.state["yaw"]

        self.state["ctrl"] = False

        # Состояние автопилота
        self.state["ap"] = self.model["autopilot"].copy()

        mesh.mesh_initialization(self.model, self.state)
        matrix.matrices_definition(self.state)
        planet.simulate_atmosphere(self.model, self.state)
        matrix.apparent_mass_matrices(self.model)
        mass.addmass(self.model, self.state)
        mass.inertia_initialization(self.model, self.state)
        aerodynamic.velocity_initialization(self.model, self.state)
        aerodynamic.HVM(self.model, self.state)
        mass.payload_loads_contribution(self.model, self.state)
        mass.simulate_aparent_mass(self.model, self.state)

        # Adding the forces and moments contributions in Body frame
        self.state["total_force"] = (self.state["aerodynamic_force"] +
                                     np.squeeze(self.state["force_payload_bodyframe"]) +
                                     np.squeeze(self.state["weight"]) +
                                     np.squeeze(self.state["apparent_mass_force"]))
        self.state["total_moment"] = (self.state["aerodynamic_moment"] +
                                      np.squeeze(self.state["moment_payload_bodyframe"]) +
                                      self.state["apparent_mass_moment"])

        w2b = self.state["w2b"]
        aerodynamic_force = np.zeros((3, 1))
        aerodynamic_force[:, 0] = self.state["aerodynamic_force"]
        Rho = self.state["Rho"]
        Sref_canopy = self.model["canopy"]["Sref"]
        Vrefn = self.state["Vrefn"]
        self.state["aeroforce"] = (-(w2b.transpose().dot(aerodynamic_force)) /
                                   (0.5*Rho*Sref_canopy*Vrefn**2))
        aerodynamic_moment = np.zeros((3, 1))
        aerodynamic_moment[:, 0] = self.state["aerodynamic_moment"]
        self.state["aeromom"] = w2b.transpose().dot(aerodynamic_moment)

        viner = self.state["i2b"].transpose().dot(self.state["Vcg_b"])
        self.state["trajv"] = np.arctan(viner[2, 0]/viner[0, 0])

        # Создаём файлы расчётов
        self.init_files()
        # Первая точка t = 0
        store.store(self.files, self.state)

    def start(self):
        # Собираем модель при первом запуске
        self.build()
        # Задаём начальные состояния
        self.init_state()
        # Цикл расчёта
        self.loop()

    def step(self):
        sim.simulate_state(self.model, self.state)
        planet.simulate_atmosphere(self.model, self.state)
        autopilot.control(self.model, self.state)
        # Simulate aerodynamics
        tol_angle = self.state["ap"]["tol_angle"]
        # incremental values
        d_alpha = self.state["alpha"] - self.state["prev"]["alpha"]
        d_pitch = self.state["pitch"] - self.state["prev"]["pitch"]
        d_roll = self.state["roll"] - self.state["prev"]["roll"]
        d_yaw = self.state["yaw"] - self.state["prev"]["yaw"]
        if abs(d_alpha) > tol_angle or abs(d_pitch) > tol_angle or abs(d_roll) > tol_angle or abs(d_yaw) > tol_angle:
            self.state["Uref_m"] = self.state["Vref_m"].transpose()
            # HVM
            aerodynamic.HVM(self.model, self.state)
            self.state["prev"]["alpha"] = self.state["alpha"]
            self.state["prev"]["pitch"] = self.state["pitch"]
            self.state["prev"]["roll"] = self.state["roll"]
            self.state["prev"]["yaw"] = self.state["yaw"]

        # Trapped air mass
        mass.addmass(self.model, self.state)
        # Simulate inertia
        mass.simulate_inertia(self.model, self.state)
        # Simulate payload forces
        mass.payload_loads_contribution(self.model, self.state)
        # Simulate apparent mass
        mass.simulate_aparent_mass(self.model, self.state)

        # Forces and moments calculation
        self.state["total_force"] = (self.state["aerodynamic_force"] +
                                     np.squeeze(self.state["force_payload_bodyframe"]) +
                                     np.squeeze(self.state["weight"]) +
                                     np.squeeze(self.state["apparent_mass_force"]))
        self.state["total_moment"] = (self.state["aerodynamic_moment"] +
                                      np.squeeze(self.state["moment_payload_bodyframe"]) +
                                      self.state["apparent_mass_moment"])

        # VECTORS
        # Forces and moments in the Wing Axes history
        w2b = np.zeros((3, 3))
        alpha = self.state["alpha"]
        sideslip_angle = self.state["sideslip_angle"]
        # Wing to Body axes
        w2b[0, 0] = np.cos(alpha)*np.cos(sideslip_angle)
        w2b[0, 1] = -np.cos(alpha)*np.sin(sideslip_angle)
        w2b[0, 2] = -np.sin(alpha)
        w2b[1, 0] = np.sin(sideslip_angle)
        w2b[1, 2] = np.cos(sideslip_angle)
        w2b[1, 2] = 0
        w2b[2, 0] = np.sin(alpha)*np.cos(sideslip_angle)
        w2b[2, 1] = -np.sin(alpha)*np.sin(sideslip_angle)
        w2b[2, 2] = np.cos(alpha)
        self.state["w2b"] = w2b

        aerodynamic_force = np.zeros((3, 1))
        aerodynamic_force[:, 0] = self.state["aerodynamic_force"]
        Rho = self.state["Rho"]
        Sref_canopy = self.model["canopy"]["Sref"]
        Vrefn = self.state["Vrefn"]
        self.state["aeroforce"] = (-(w2b.transpose().dot(aerodynamic_force)) /
                                   (0.5*Rho*Sref_canopy*Vrefn**2))
        aerodynamic_moment = np.zeros((3, 1))
        aerodynamic_moment[:, 0] = self.state["aerodynamic_moment"]
        self.state["aeromom"] = w2b.transpose().dot(aerodynamic_moment)

        viner = self.state["i2b"].transpose().dot(self.state["Vcg_b"])
        self.state["trajv"] = np.arctan(viner[2, 0]/viner[0, 0])

        if self.state["altitude"] < 0:
            self.state["altitude"] = 0

        # В конце каждого шага сохраняем данные
        store.store(self.files, self.state)

    def loop(self):
        while self.state["time"] < self.model["time"]["final"] and self.state["altitude"] > 0:
            self.step()
        # Закрываем все открытые файлы
        for key in self.files:
            self.files[key].close()
    
    def closeFiles(self):
        for key in self.files:
            self.files[key].close()
