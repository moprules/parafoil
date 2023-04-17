import numpy as np
from . import planet
from . import matrix


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
    state["Xam"] = model["canopy"]["Xam"]
    state["mass"] = model["body"]["mass"]
    I = model["body"]["iner"]
    moments_of_inertia = [[I["xx"], -I["xy"], -I["xz"]],
                          [-I["xy"], I["yy"], -I["yz"]],
                          [-I["xz"], -I["yz"], I["zz"]]]
    state["moments_of_inertia"] = np.asarray(moments_of_inertia)
    simulate_inertia(model, state)


def payload_loads_contribution(model: dict, state: dict):
    canopy = model["canopy"]
    # Payload velicity calculation
    SWB = matrix.cross_product_matrix(state["Wb"])
    SGPL = matrix.cross_product_matrix(canopy["Xcg_payload"])
    velocity_payload = (SWB.dot(canopy["Xcg_payload"])
                        + (state["Vcg_b"] - state["i2b"].dot(state["Vwind"])))

    # Payload dynamic pressure
    [q] = 0.5*state["Rho"]*np.sqrt(velocity_payload[0] **
                                   2 + velocity_payload[1]**2 + velocity_payload[2]**2)
    pass

    # Payload force calculation
    force_payload_bodyframe = np.zeros((3, 1))
    # Payload force component X
    force_payload_bodyframe[0] = (-q*model["body"]["CD"] *
                                  model["body"]["Sref"]*velocity_payload[0])
    # Payload force component Y
    force_payload_bodyframe[1] = (-q*model["body"]["CD"] *
                                  model["body"]["Sref"]*velocity_payload[1])
    # Payload force component Z
    force_payload_bodyframe[2] = (-q*model["body"]["CD"] *
                                  model["body"]["Sref"]*velocity_payload[2])

    state["force_payload_bodyframe"] = force_payload_bodyframe

    # Payload moment calculation
    # Matrix used for the cross product
    SCGP = matrix.cross_product_matrix(canopy["Xcg_payload"])
    # r X F product
    moment_payload_bodyframe = (SCGP.dot(force_payload_bodyframe)
                                + SGPL.dot(state["weight"]))
    state["moment_payload_bodyframe"] = moment_payload_bodyframe


def simulate_aparent_mass(model: dict, state: dict):
    # Canopy velocities calculation
    angular_velocity_canopy = state["b2c"].dot(state["Wb"])

    SWB = matrix.cross_product_matrix(state["Wb"])
    VAUX = SWB.dot(state["Xcg_aero"]) + state["Vcg_b"]
    velocity_canopy = state["b2c"].dot(VAUX-(state["i2b"].dot(state["Vwind"])))

    # Calculation of the matrices used for the cross products
    SWC = matrix.cross_product_matrix(angular_velocity_canopy)
    SCGM = matrix.cross_product_matrix(state["Xam"])
    SVAC = matrix.cross_product_matrix(velocity_canopy)

    # Second rank tensors -> I_X_bis = TBC'*I_X*TBC
    IAM_bis = np.dot(state["b2c"].transpose(), model["I_AM"]).dot(state["b2c"])
    IAI_bis = np.dot(state["b2c"].transpose(), model["I_AI"]).dot(state["b2c"])
    IH_bis = np.dot(state["b2c"].transpose(), model["I_H"]).dot(state["b2c"])

    # Apparent mass force contributions
    VAUX = (model["I_AM"].dot(velocity_canopy) +
            model["I_H"].dot(angular_velocity_canopy))
    apparent_mass_force = (-state["b2c"].transpose().dot(SWC.dot(VAUX)) -
                           (np.dot(IAM_bis, SWB).dot(state["i2b"])).dot(state["Vwind"]))
    state["apparent_mass_force"] = apparent_mass_force

    # Aparent mass moment contributions
    # -(-SCGM*(b2c'*(SWC*VAUX)) - b2c'*(SWC*(I_H*velocity_canopy)) + b2c'*(((SWC*I_AI) + (SVAC*I_H))*angular_velocity_canopy)- ((((SCGM*IAM_bis)+IH_bis)*SWC)*i2b)*Vwind).*[0;1;0]
    state["apparent_mass_moment"] = 0

    # ASSEMBLING APPARENT MASS TERMS INTO GLOBAL INERTIA MATRIX
    MAUX1 = IH_bis - IAM_bis.dot(SCGM)
    MAUX2 = SCGM.dot(IAM_bis) + IH_bis
    MAUX3 = IAI_bis - IH_bis.dot(SCGM) + (SCGM.dot(IH_bis)-IAM_bis).dot(SCGM)

    iner_matrix = state["inertia_matrix"]
    for i in range(3):
        for j in range(3):
            iner_matrix[i, j] = iner_matrix[i, j] + IAM_bis[i, j]
            iner_matrix[i+3, j] = iner_matrix[i+3, j] + MAUX1[i, j]
            iner_matrix[i, j+3] = iner_matrix[i, j+3] + MAUX2[i, j]
            iner_matrix[i+3, j+3] = iner_matrix[i+3, j+3] + MAUX3[i, j]
