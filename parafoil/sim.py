from . import matrix
import numpy as np


def get_precision(f):
    str_f = str(f)
    if '.' not in str_f:
        return 0

    # Получение строки после точки и возвращение ее длины
    return len(str_f[str_f.index('.') + 1:])


def state_derivatives(Xe, state):
    U, V, W = Xe[:3]
    PHI, THETA, PSI = Xe[3:6]
    P, Q, R = Xe[6:9]
    FX, FY, FZ = state["total_force"]
    MX, MY, MZ = state["total_moment"]
    mass = state["mass"]

    # Auxiliar variables
    STH = np.sin(THETA)
    CTH = np.cos(THETA)
    SPH = np.sin(PHI)
    CPH = np.cos(PHI)
    SPSI = np.sin(PSI)
    CPSI = np.cos(PSI)

    # Force equations
    force_vector = np.zeros((6, 1))
    force_vector[0, 0] = (R*V - Q*W)*mass + FX
    force_vector[1, 0] = (P*W - R*U)*mass + FY
    force_vector[2, 0] = (Q*U - P*V)*mass + FZ

    # Angular velocity vector
    VAUX1 = np.asarray([[P], [Q], [R]])
    SWB = matrix.cross_product_matrix(VAUX1)

    VAUX1 = np.dot(SWB, state["moments_of_inertia"]).dot(VAUX1)

    # Moment equations
    force_vector[3, 0] = MX - VAUX1[0, 0]
    force_vector[4, 0] = MY - VAUX1[1, 0]
    force_vector[5, 0] = MZ - VAUX1[2, 0]

    # SOLVING THE SYSTEM
    solution_vector = np.linalg.solve(state["inertia_matrix"], force_vector)
    XDe = np.zeros(12)
    XDe[0] = solution_vector[0, 0]
    XDe[1] = solution_vector[1, 0]
    XDe[2] = solution_vector[2, 0]
    XDe[6] = solution_vector[3, 0]
    XDe[7] = solution_vector[4, 0]
    XDe[8] = solution_vector[5, 0]

    # KINEMATIC EQUATIONS
    XDe[3] = P + (STH/CTH)*(Q*SPH + R*CPH)
    XDe[4] = Q*CPH - R*SPH
    XDe[5] = (Q*SPH + R*CPH)/CTH

    # NAVIGATION EQUATIONS
    XDe[9] = (U*CTH*CPSI + V*(SPH*STH*CPSI - CPH*SPSI) +
              W*(CPH*STH*CPSI + SPH*SPSI))
    XDe[10] = (U*CTH*SPSI + V*(SPH*STH*SPSI + CPH*CPSI) +
               W*(CPH*STH*SPSI - SPH*CPSI))
    XDe[11] = U*STH - V*SPH*CTH - W*CPH*CTH

    return XDe


def EULER(T, DT, Xe: list, state: dict):
    XDe = state_derivatives(Xe, state)
    for i in range(len(Xe)):
        Xe[i] = Xe[i]+XDe[i]*DT

    T = T + DT

    return T, Xe


def simulate_state(model: dict, state: dict):
    Xe = [state["Vcg_b"][0, 0],
          state["Vcg_b"][1, 0],
          state["Vcg_b"][2, 0],
          state["roll"],
          state["pitch"],
          state["yaw"],
          state["Wb"][0, 0],
          state["Wb"][1, 0],
          state["Wb"][2, 0],
          state["pos_north"],
          state["pos_east"],
          state["altitude"]]

    T = state["time"]
    DT = model["time"]["step"]

    T, Xe = EULER(T, DT, Xe, state)

    # Trajectory integration
    state["dtrav_dist"] = ((Xe[9]-state["pos_north"])**2 +
                           (Xe[10]-state["pos_east"])**2 +
                           (Xe[11]-state["altitude"])**2)**0.5
    state["traveled_distance"] = state["traveled_distance"] + state["dtrav_dist"]

    state["dlost_altitude"] = (Xe[11] - state["altitude"])
    state["lost_altitude"] = state["lost_altitude"] + state["dlost_altitude"]

    [state["Vcg_b"][0, 0],
     state["Vcg_b"][1, 0],
     state["Vcg_b"][2, 0],
     state["roll"],
     state["pitch"],
     state["yaw"],
     state["Wb"][0, 0],
     state["Wb"][1, 0],
     state["Wb"][2, 0],
     state["pos_north"],
     state["pos_east"],
     state["altitude"]] = Xe

    state["time"] = round(T, get_precision(DT))

    matrix.matrices_definition(state)

    # Reference velocity definition
    # Vwind_b = i2b*Vwind;
    state["Vwind_b"] = state["i2b"].dot(state["Vwind"])
    state["Vref_b"] = state["Vwind_b"]-state["Vcg_b"]
    state["Vrefn"] = np.sqrt(state["Vref_b"][0]**2 +
                             state["Vref_b"][1]**2+state["Vref_b"][2]**2)
    state["Vref_m"] = state["m2c"].transpose().dot(
        state["b2c"].dot(state["Vref_b"]))
    state["alpha"] = np.arctan(state["Vref_b"][2]/state["Vref_b"][0])
    state["sideslip_angle"] = np.arcsin(state["Vref_b"][1]/state["Vrefn"])

    MAUX1 = matrix.cross_product_matrix(state["Wb"])
    VAUX = MAUX1.dot(state["Xcg_aero"]) + state["Vcg_b"]
    # velocity_canopy = state["b2c"].dot(VAUX-state["Vwind_b"])

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
        r = (Xref_b + b2c.transpose().dot(m2c.dot(Xref_m)) +
             b2c.transpose().dot(m2c.dot(xctrl[:, [k]])))
        r2 = (Xref_b + b2c.transpose().dot(m2c.dot(Xref_m)) +
              b2c.transpose().dot(m2c.dot(xbound[:, [k]])))
        # Linear velocity due to rotation expressed in body frame
        Vrot_b = Swb.dot(r)
        # Linear velocity due to rotation expressed in body frame
        Vrot_b2 = Swb.dot(r2)
        # Transformation of Vrot to matlab axes
        Vrot_m = m2c.transpose().dot(b2c.dot(Vrot_b))
        # Transformation of Vrot to matlab axes
        Vrot_m2 = m2c.transpose().dot(b2c.dot(Vrot_b2))
        # Aerodynamic velocity seen by the current control point (matlab axes)
        mesh["Vinf"][k, :] = np.squeeze(Vref_m - Vrot_m)
        mesh["Vinf2"][k, :] = np.squeeze(Vref_m - Vrot_m2)

    cr = model["canopy"]["root_chord"]
    # Aerodynamic center coordinates
    state["Xcg_aero"] = Xref_b + state["b2c"].transpose().dot(
        state["m2c"].dot(np.asarray([[cr/4], [0.0], [0.0]])))
    state["Xam"] = np.asarray([[0.0], [0.0], [0.0]])
