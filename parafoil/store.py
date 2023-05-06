import os
import io
import datetime
import shutil
import numpy as np
# Директория пакета
PKG_DIR = os.path.dirname(os.path.abspath(__file__))


def init_files(files: dict = {}, path_to_res: str = "", path_to_model: str = ""):
    path_to_res = os.path.abspath(path_to_res)
    model_prefix = os.path.basename(path_to_model)
    model_prefix = os.path.splitext(model_prefix)[0]
    now = datetime.datetime.now()
    res_folder = model_prefix + now.strftime("%Y-%m-%d_%H.%M")
    res_folder = os.path.join(path_to_res, res_folder)
    if os.path.exists(res_folder):
        shutil.rmtree(res_folder)
    shutil.copytree(os.path.join(PKG_DIR, "templates"), res_folder)

    # Заполняем словарь файлов обхектами файлов
    # Чтобы потом через словарь можно было дозаписывать данные
    for file in os.listdir(res_folder):
        key = os.path.splitext(file)[0]
        files[key] = open(os.path.join(res_folder, file), "a")


def init_coords(files: dict[str, io.TextIOWrapper], model: dict):
    """Инициализируем узел coords и доп узлы если нужно"""
    # Нужно в файл трёх мерной траектории дописать узлы
    tr_3d = files["tr_3d"]

    # Строка разделитель
    print(file=tr_3d)
    # Узел cube
    print("cube:", file=tr_3d)
    # Точка с координатама самого тела
    x = model["body"]["north"]
    y = model["body"]["east"]
    z = model["body"]["altitude"]
    print(x, y, z, file=tr_3d)
    # Точка с координатами цели
    x = model["target"]["north"]
    y = model["target"]["east"]
    z = model["target"]["altitude"]
    print(x, y, z, file=tr_3d)
    # Пустая строка разделитель
    print(file=tr_3d)

    # Узел areas
    print("areas:", file=tr_3d)
    radius = model["target"]["radius"]
    print(radius, x, y, z, file=tr_3d)

    # Траектория в плоскости xy
    tr_xy = files["tr_xy"]
    # Строка разделитель
    print(file=tr_xy)
    # Узел areas
    print("areas:", file=tr_xy)
    radius = model["target"]["radius"]
    print(radius, x, y, z, file=tr_xy)


    # Для всех файлов записываем узел coords
    for key in files:
        # Строка разделитель
        print(file=files[key])
        # Узел coords
        print("coords:", file=files[key])


def store(files: dict[str, io.TextIOWrapper], state: dict):
    t = state["time"]
    t = f"{t:<7.3f}"
    x = state["pos_north"]
    y = state["pos_east"]
    z = state["altitude"]
    alpha = state["alpha"]
    beta = state["sideslip_angle"]
    gamma = state["trajv"]
    roll = state["roll"]
    pitch = state["pitch"]
    yaw = state["yaw"]
    velocity = state["Vrefn"]
    aeroforce = state["aeroforce"]
    C_D = aeroforce[0,0]
    C_Y = aeroforce[1,0]
    C_L = aeroforce[2,0]

    print(f"{t} -> {x} {y} {z}", file=files["tr_3d"])
    print(f"{x} -> {y}", file=files["tr_xy"])
    print(f"{x} -> {z}", file=files["tr_xz"])
    print(f"{y} -> {z}", file=files["tr_yz"])
    print(f"{t} -> {np.rad2deg(alpha)}", file=files["alpha"])
    print(f"{t} -> {np.rad2deg(beta)}", file=files["beta"])
    print(f"{t} -> {np.rad2deg(gamma)}", file=files["gamma"])
    print(f"{t} -> {np.rad2deg(roll)}", file=files["roll"])
    print(f"{t} -> {np.rad2deg(pitch)}", file=files["pitch"])
    print(f"{t} -> {np.rad2deg(yaw)}", file=files["yaw"])
    print(f"{t} -> {velocity}", file=files["velocity"])
    print(f"{t} -> {C_D}", file=files["C_D"])
    print(f"{t} -> {C_L}", file=files["C_L"])
    print(f"{t} -> {C_Y}", file=files["C_Y"])
