import numpy as np


def mesh_initialization(model: dict, state: dict):
    # Параметры парашюта
    canopy = model["canopy"]
    mesh = canopy["mesh"].copy()

    # Кличество узлов парашюта при разделении на панели
    # Это количество панелей + 1
    mesh["nn"] = mesh["N"] + 1

    tan_sw = np.tan(canopy["AngLE"])
    mesh["coord"] = np.zeros((3, mesh["nn"]))
    c = np.zeros(mesh["nn"])
    for i in range(mesh["nn"]):
        theta = np.pi * i / mesh["N"]
        mesh["coord"][2][i] = -canopy["height"] * abs(1-np.sin(theta))
        mesh["coord"][1][i] = -0.5 * canopy["span"] * np.cos(theta)
        c[i] = canopy["root_chord"] - 2 * abs(mesh["coord"][1][i]) * (
            canopy["root_chord"] - canopy["tip_chord"])/(0.5*canopy["span"])
        mesh["coord"][0][i] = abs(mesh["coord"][1][i]) * tan_sw

    mesh["xctrl"] = np.zeros((3, mesh["N"]))
    mesh["xbound"] = np.zeros((3, mesh["N"]))
    mesh["len"] = np.zeros(mesh["N"])
    mesh["s"] = np.zeros(mesh["N"])
    c_ave = np.zeros(mesh["N"])
    for i in range(mesh["N"]):
        mesh["xbound"][:, i] = (mesh["coord"][:, i] +
                                mesh["coord"][:, i+1])*0.5
        c_ave[i] = (c[i] + c[i+1])*0.5
        mesh["xctrl"][:, i] = (mesh["xbound"][:, i] +
                               np.asarray([c_ave[i]*0.5, 0, 0]))
        delta = mesh["coord"][:, i+1] - mesh["coord"][:, i]
        mesh["len"][i] = np.sqrt(delta.dot(delta))
        mesh["s"][i] = mesh["len"][i]*c_ave[i]

    mesh["alphalo"] = np.zeros(mesh["N"])
    mesh["ypos"] = np.zeros(mesh["N"])
    a0 = canopy["a0"]
    for i in range(mesh["N"]):
        mesh["ypos"][i] = 2.0 * abs(mesh["xbound"][1,i])/canopy["span"]
        mesh["alphalo"][i] = (
            a0["root"] + (a0["tip"]-a0["root"]) * mesh["ypos"][i])

    # for param in ("coord", "xctrl", "xbound", "len", "s", "alphalo", "ypos"):
    #     mesh[param] = np.ascontiguousarray(mesh[param])

    state["mesh"] = mesh
