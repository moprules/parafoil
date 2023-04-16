from . import builder
from . import matrix
from . import planet
from . import mass
from . import mesh
from . import aerodynamic
import yaml
import numpy as np


class PFSim:
    """Класс - симулятор парафойла"""

    def __init__(self, path_to_data: str = ""):
        self.params = {}
        self.state = {}
        self.model = {}

        if path_to_data:
            self.load(path_to_data)

    def load(self, path_to_data: str):
        """Загрузка параметров из исходного файла"""
        with open(path_to_data, "r") as stream:
            self.params = yaml.safe_load(stream)

    def build(self):
        self.model = builder.build(self.params)

    def init_state(self):
        # Новый словарь состояний
        self.state = {}
        self.state["time"] = self.model["time"]["init"]
        self.state["roll"] = self.model["body"]["roll"]
        self.state["pitch"] = self.model["body"]["pitch"]
        self.state["yaw"] = self.model["body"]["yaw"]
        self.state["rigging_angle"] = self.model["canopy"]["rigging_angle"]
        self.state["altitude"] = self.model["body"]["altitude"]
        self.state["east"] = self.model["body"]["east"]
        self.state["north"] = self.model["body"]["north"]
        # Снижение высоты корабля
        self.state["lost_altitude"] = 0
        # Пройденное расстояние
        self.state["traveled_distance"] = 0
        # Значение параметров на предыдущем шаге
        self.state["prev"] = {}
        self.state["prev"]["roll"] = self.state["roll"]
        self.state["prev"]["pitch"] = self.state["pitch"]
        self.state["prev"]["yaw"] = self.state["yaw"]

        mesh.mesh_initialization(self.model, self.state)
        matrix.matrices_definition(self.state)
        planet.simulate_atmosphere(self.model, self.state)
        matrix.apparent_mass_matrices(self.model)
        mass.addmass(self.model, self.state)
        mass.inertia_initialization(self.model, self.state)
        aerodynamic.velocity_initialization(self.model, self.state)
        aerodynamic.HVM(self.model, self.state)
        mass.payload_loads_contribution(self.model, self.state)

    def start(self):
        # Собираем модель при первом запуске
        self.build()
        # Задаём начальные состояния
        self.init_state()

    def step(self):
        pass

    def loop(self):
        while self.state["time"] < self.model["time"]["final"] and self.state["altitude"] > 0:
            self.step()

    def __str__(self):
        return str(self.params)
