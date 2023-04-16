import numpy as np
from . import planet


def addmass(model, state):
    Rho = state["Rho"]
    t = model["canopy"]["thicknes"]
    c = model["canopy"]["root_chord"]
    b = model["canopy"]["span"]
    rigging_angle = state["rigging_angle"]
    R = 0.6*b
    E = np.arcsin(b/(2*R))
    zeta = 0.25
    xi = c/np.sqrt(1+rigging_angle ** 2)
    state["trapped"] = Rho*t*c*E*(2*R+(1-2*zeta)*xi)


def simulate_inertia(model, state):
    g = planet.gravity(model, state)
    moments_of_inertia: np.ndarray = state["moments_of_inertia"]
    i2b: np.ndarray = state["i2b"]
    mass = state["mass"]
    trapped = state["trapped"]
    # Вес в инерциальной системе координат
    weight = [[0.0], [0.0], [(mass+trapped)*g]]
    # Вес в системе координат тела
    weight = i2b.dot(weight)
    state["weight"] = weight
    # Матрица инерций
    inertia_matrix = np.zeros((6, 6))

    inertia_matrix[0, 0] = mass + trapped
    inertia_matrix[1, 1] = mass + trapped
    inertia_matrix[2, 2] = mass + trapped
    inertia_matrix[3, 3] = moments_of_inertia[0, 0]
    inertia_matrix[3, 4] = -moments_of_inertia[0, 1]
    inertia_matrix[3, 5] = -moments_of_inertia[0, 2]
    inertia_matrix[4, 3] = -moments_of_inertia[0, 1]
    inertia_matrix[4, 4] = moments_of_inertia[1, 1]
    inertia_matrix[4, 5] = -moments_of_inertia[1, 2]
    inertia_matrix[5, 3] = -moments_of_inertia[0, 2]
    inertia_matrix[5, 4] = -moments_of_inertia[1, 2]
    inertia_matrix[5, 5] = moments_of_inertia[2, 2]
    state["inertia_matrix"] = inertia_matrix


def inertia_initialization(model, state):
    mass = model["body"]["mass"]
    state["mass"] = mass
    I = model["body"]["iner"]
    moments_of_inertia = [[I["xx"], -I["xy"], -I["xz"]],
                          [-I["xy"], I["yy"], -I["yz"]],
                          [-I["xz"], -I["yz"], I["zz"]]]
    state["moments_of_inertia"] = np.asarray(moments_of_inertia)
    simulate_inertia(model, state)

