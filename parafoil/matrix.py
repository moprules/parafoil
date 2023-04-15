import numpy as np


def apparent_mass_matrices(model: dict):
    a = model["canopy"]["height"]
    b = model["canopy"]["span"]
    t = model["canopy"]["thicknes"]
    cref = model["canopy"]["ref_chord"]
    Sref_canopy = model["canopy"]["Sref"]
    AR = model["canopy"]["AR"]
    A = 0.66*(1+(8/3)*(a/b)**2)*t**2*b
    B = 0.267*(t**2 + 2*(a**2)*(1 - (t/cref)**2))*cref

    inv = 1/(1+AR)
    sroot = np.sqrt(1+2*((a/b)**2)*(1-(t/cref)**2))
    C = 0.785*sroot*AR*inv*b*cref**2

    P = 0.055*AR*inv*b*Sref_canopy**2
    Q = 0.0308*AR*inv*cref**3*Sref_canopy * \
        (1+(np.pi/6)*(1+AR)*AR*((a/b)**2)*((t/cref)**2))
    R = 0.0555*b**3*t**2*(1+8*(a/b)**2)

    # Write here the value of H coefficient if spanwise effect wants to be simulated
    H = 0.2

    # LAMB MATRIX I_AM
    I_AM = [[A, 0.0, 0.0],
            [0.0, B, 0.0],
            [0.0, 0.0, C]]
    model["I_AM"] = np.asarray(I_AM)

    # LAMB MATRIX I_AI
    I_AI = [[P, 0.0, 0.0],
            [0.0, Q, 0.0],
            [0.0, 0.0, R]]
    model["I_AI"] = np.asarray(I_AI)

    # LAMB MATRIX I_H
    I_H = [[0.0, H, 0.0],
           [H, 0.0, 0.0],
           [0.0, 0.0, 0.0]]
    model["I_H"] = np.asarray(I_H)

    return I_AM, I_AI, I_H


def cross_product_matrix(vector: np.ndarray) -> np.ndarray:
    vector = np.squeeze(vector)
    M = [[0.0, -vector[2], vector[1]],
         [vector[2], 0.0, -vector[0]],
         [-vector[1], vector[0], 0.0]]
    return np.asarray(M)


def matrices_definition(state: dict):
    F = state["roll"]
    T = state["pitch"]
    P = state["yaw"]
    rigging_angle = state["rigging_angle"]

    i2b = np.zeros((3, 3))
    i2b[0, 0] = np.cos(T)*np.cos(P)
    i2b[0, 1] = np.cos(T)*np.sin(P)
    i2b[0, 2] = -np.sin(T)
    i2b[1, 0] = np.sin(F)*np.sin(T)*np.cos(P) - np.cos(F)*np.sin(P)
    i2b[1, 1] = np.sin(F)*np.sin(T)*np.sin(P) + np.cos(F)*np.cos(P)
    i2b[1, 2] = np.sin(F)*np.cos(T)
    i2b[2, 0] = np.cos(F)*np.sin(T)*np.cos(P) + np.sin(F)*np.sin(P)
    i2b[2, 1] = np.cos(F)*np.sin(T)*np.sin(P) - np.sin(F)*np.cos(P)
    i2b[2, 2] = np.cos(F)*np.cos(T)
    state["i2b"] = i2b

    m2c = [[-1.0, 0.0, 0.0],
           [0.0, 1.0, 0.0],
           [0.0, 0.0, -1.0]]
    state["m2c"] = np.asarray(m2c)

    b2c = [[np.cos(rigging_angle), 0.0, -np.sin(rigging_angle)],
           [0.0,                   1.0, 0.0],
           [np.sin(rigging_angle), 0.0, np.cos(rigging_angle)]]
    state["b2c"] = np.asarray(b2c)
