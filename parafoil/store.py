import os
import io
import datetime
import shutil
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
    # Пустая строка разделитель
    print(file=tr_3d)

    # Для всех файлов записываем узел coords
    for key in files:
        print("coords:", file=files[key])


def store(files: dict[str, io.TextIOWrapper], state: dict):
    tr_3d = files["tr_3d"]
    t = state["time"]
    x = state["pos_north"]
    y = state["pos_east"]
    z = state["altitude"]
    print(f"{t:<7.3f} -> {x} {y} {z}", file=tr_3d)
