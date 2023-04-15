import numpy as np


def simulate_atmosphere(model, state):
    Rho0 = model["planet"]["Rho0"]
    h = state["altitude"]
    if h <= 0:
        Rho = Rho0
    elif h >= 0 and h < 11000:
        theta = 1-22.557*(1e-6)*h
        Rho = Rho0*(theta**4.256)
    else:
        Rho = (0.297*np.exp(-157.69*(h-11000)*(1e-6)))*Rho0

    state["Rho"] = Rho


def gravity(model, state):
    R = model["planet"]["R"]
    g0 = model["planet"]["g0"]
    h = state["altitude"]
    dummy1 = 2*h/R
    dummy2 = (3*(h)**2)/R**2
    state["g"] = g0*(1 - dummy1 + dummy2)
    return state["g"]
