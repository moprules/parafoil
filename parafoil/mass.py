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
    mass = model["body"]["mass"]
    state["mass"] = mass
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
    velocity_payload = SWB.dot(
        canopy["Xcg_payload"]) + (state["Vcg_b"] - state["i2b"].dot(state["Vwind"]))

    # Payload dynamic pressure
    [q] = 0.5*state["Rho"]*np.sqrt(velocity_payload[0] **
                                   2 + velocity_payload[1]**2 + velocity_payload[2]**2)
    pass

    # Payload force calculation
    force_payload_bodyframe = np.zeros((3, 1))
    # Payload force component X
    force_payload_bodyframe[0] = -q*model["body"]["CD"] * \
        model["body"]["Sref"]*velocity_payload[0]
    # Payload force component Y
    force_payload_bodyframe[1] = -q*model["body"]["CD"] * \
        model["body"]["Sref"]*velocity_payload[1]
    # Payload force component Z
    force_payload_bodyframe[2] = -q*model["body"]["CD"] * \
        model["body"]["Sref"]*velocity_payload[2]

    state["force_payload_bodyframe"] = force_payload_bodyframe

    # Payload moment calculation
    # Matrix used for the cross product
    SCGP = matrix.cross_product_matrix(canopy["Xcg_payload"])
    # r X F product
    moment_payload_bodyframe = SCGP.dot(
        force_payload_bodyframe) + SGPL.dot(state["weight"])
    state["moment_payload_bodyframe"] = moment_payload_bodyframe


# function [apparent_mass_force,apparent_mass_moment,inertia_matrix] = simulate_aparent_mass(Vcg_b,Wb,Xam,Xcg_aero,b2c,i2b,I_AM,I_AI,I_H,Vwind,inertia_matrix)
    
#     %Canopy velocities calculation
#     angular_velocity_canopy = b2c*Wb;
    
#     SWB = cross_product_matrix(Wb);
#     VAUX = SWB*Xcg_aero + Vcg_b;   
#     velocity_canopy = b2c*(VAUX-(i2b*Vwind));
#     ANG = angular_velocity_canopy;
        
#     %Calculation of the matrices used for the cross products
#     SWC = cross_product_matrix(angular_velocity_canopy);
#     SCGM = cross_product_matrix(Xam);
#     SVAC = cross_product_matrix(velocity_canopy);
    
#     %Second rank tensors -> I_X_bis = TBC'*I_X*TBC
#     IAM_bis = b2c'*I_AM*b2c;
#     IAI_bis = b2c'*I_AI*b2c;
#     IH_bis  = b2c'*I_H*b2c;
    
#     % Apparent mass force contributions
#     VAUX = I_AM*velocity_canopy + I_H*angular_velocity_canopy;    
# 	apparent_mass_force = -b2c'*(SWC*VAUX)-((IAM_bis*SWB)*i2b)*Vwind;
	                      
#     %Aparent mass moment contributions   
#     apparent_mass_moment = 0;% -(-SCGM*(b2c'*(SWC*VAUX)) - b2c'*(SWC*(I_H*velocity_canopy)) + b2c'*(((SWC*I_AI) + (SVAC*I_H))*angular_velocity_canopy)- ((((SCGM*IAM_bis)+IH_bis)*SWC)*i2b)*Vwind).*[0;1;0];
    
#     %ASSEMBLING APPARENT MASS TERMS INTO GLOBAL INERTIA MATRIX
#     MAUX1 = IH_bis - IAM_bis*SCGM;
#     MAUX2 = SCGM*IAM_bis + IH_bis;
#     MAUX3 = IAI_bis - IH_bis*SCGM + (SCGM*IH_bis-IAM_bis)*SCGM;
    
#     for i=1:3
#         for j=1:3
#             inertia_matrix(i,j) = inertia_matrix(i,j) + IAM_bis(i,j);
#             inertia_matrix(i+3,j) = inertia_matrix(i+3,j) + MAUX1(i,j);
#             inertia_matrix(i,j+3) = inertia_matrix(i,j+3) + MAUX2(i,j);
#             inertia_matrix(i+3,j+3) = inertia_matrix(i+3,j+3) + MAUX3(i,j);
#         end
#     end
# end
