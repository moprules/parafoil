import numpy as np
from . import matrix
from pprint import pprint


def velocity_initialization(model: dict, state: dict):
    state["Vcg_b"] = model["body"]["linear_velocity"].copy()
    state["Wb"] = model["body"]["angular_velocity"].copy()
    state["Vwind"] = model["wind"].copy()

    state["Vwind_b"] = state["i2b"].dot(state["Vwind"])
    state["Vref_b"] = state["Vwind_b"]-state["Vcg_b"]
    state["Vrefn"] = np.sqrt(state["Vref_b"][0]**2 +
                             state["Vref_b"][1]**2+state["Vref_b"][2]**2)
    state["Vref_m"] = state["m2c"].transpose().dot(
        state["b2c"].dot(state["Vref_b"]))
    state["alpha"] = np.arctan(state["Vref_b"][2]/state["Vref_b"][0])
    state["sideslip_angle"] = np.arcsin(state["Vref_b"][1]/state["Vrefn"])

    cr = model["canopy"]["root_chord"]
    state["Xcg_aero"] = model["canopy"]["Xref_b"] + \
        state["b2c"].transpose().dot(state["m2c"].dot(
            np.asarray([[cr/4], [0.0], [0.0]])))

    w2b = np.zeros((3, 3))

    # Wing to Body axes
    w2b[0, 0] = np.cos(state["alpha"])*np.cos(state["sideslip_angle"])
    w2b[0, 1] = -np.cos(state["alpha"])*np.sin(state["sideslip_angle"])
    w2b[0, 2] = -np.sin(state["alpha"])
    w2b[1, 0] = np.sin(state["sideslip_angle"])
    w2b[1, 2] = np.cos(state["sideslip_angle"])
    w2b[1, 2] = 0
    w2b[2, 0] = np.sin(state["alpha"])*np.cos(state["sideslip_angle"])
    w2b[2, 1] = -np.sin(state["alpha"])*np.sin(state["sideslip_angle"])
    w2b[2, 2] = np.cos(state["alpha"])
    state["w2b"] = w2b

    mesh = state["mesh"]
    # Aerodynamic speed calculation at each of the control points
    # Control point velocity
    mesh["Vinf"] = np.zeros((mesh["N"], 3))
    # Bound vortex velocity
    mesh["Vinf2"] = np.zeros((mesh["N"], 3))
    Swb = matrix.cross_product_matrix(state["Wb"])
    Xref_b = model["canopy"]["Xref_b"]
    Xref_m = model["canopy"]["Xref_m"]
    Vref_m = state["Vref_m"]
    xctrl = mesh["xctrl"]
    xbound = mesh["xbound"]
    b2c = state["b2c"]
    m2c = state["m2c"]
    for k in range(mesh["N"]):
        # r = Xref_b + b2c'*(m2c*Xref_m) + b2c'*(m2c*xctrl(:,k))
        r = (Xref_b + b2c.transpose().dot(m2c.dot(Xref_m)) +
             b2c.transpose().dot(m2c.dot(xctrl[:, [k]])))
        r2 = Xref_b + b2c.transpose().dot(m2c.dot(Xref_m)) + \
            b2c.transpose().dot(m2c.dot(xbound[:, [k]]))
        # Linear velocity due to rotation expressed in body frame
        Vrot_b = Swb.dot(r)
        # Linear velocity due to rotation expressed in body frame
        Vrot_b2 = Swb.dot(r2)
        # Transformation of Vrot to matlab axes
        Vrot_m = m2c.dot(b2c.dot(Vrot_b))
        # Transformation of Vrot to matlab axes
        Vrot_m2 = m2c.dot(b2c.dot(Vrot_b2))
        # Aerodynamic velocity seen by the current control point (matlab axes)
        mesh["Vinf"][k, :] = np.squeeze(Vref_m - Vrot_m)
        mesh["Vinf2"][k, :] = np.squeeze(Vref_m - Vrot_m2)
