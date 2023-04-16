import numpy as np
from . import matrix


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
        r = (Xref_b + b2c.transpose().dot(m2c.dot(Xref_m)) +
             b2c.transpose().dot(m2c.dot(xctrl[:, [k]])))
        r2 = (Xref_b + b2c.transpose().dot(m2c.dot(Xref_m)) +
              b2c.transpose().dot(m2c.dot(xbound[:, [k]])))
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

    state["delta0_f"] = np.zeros(2)
    state["delta0_f"][0] = model["canopy"]["cs"]["deflection"]["L"]
    state["delta0_f"][1] = model["canopy"]["cs"]["deflection"]["R"]
    state["Uref_m"] = state["Vref_m"].transpose()


def vortxl(X1, X2, XP, gamma):
    x1 = X1[0]
    y1 = X1[1]
    z1 = X1[2]
    x2 = X2[0]
    y2 = X2[1]
    z2 = X2[2]
    xp = XP[0]
    yp = XP[1]
    zp = XP[2]
    r0 = X2-X1
    r1 = XP-X1
    r2 = XP-X2
    norm_r1 = np.linalg.norm(r1)
    norm_r2 = np.linalg.norm(r2)

    r1xr2 = np.cross(r1, r2)
    norm_r1xr2 = np.linalg.norm(r1xr2)

    inv_r1xr2 = 1.0 / norm_r1xr2
    inv_r1 = 1.0 / norm_r1
    inv_r2 = 1.0 / norm_r2

    a = r0 * inv_r1xr2
    b = r1*inv_r1 - r2*inv_r2
    c = np.dot(a, b)

    u = gamma*0.25/np.pi*c*np.dot(r1xr2, inv_r1xr2)

    return u


def HVM(model: dict, state: dict):
    angh = model["canopy"]["cs"]["angh"]
    delta0_f = state["delta0_f"]
    incalphaloflap_L = (delta0_f[0]/np.pi)*(np.pi-angh+np.sin(angh))
    incalphaloflap_R = (delta0_f[1]/np.pi)*(np.pi-angh+np.sin(angh))
    yflap = model["canopy"]["cs"]["yflap"]
    mesh = state["mesh"]
    normals = np.zeros((3, mesh["N"]))
    for i in range(mesh["N"]):
        curr_a0 = mesh["alphalo"][i]
        # Left wing
        if mesh["xbound"][1, i] <= 0.0:
            if (-mesh["ypos"][i] <= yflap[1]) and (-mesh["ypos"][i] >= yflap[0]):
                # Flap delta alpha0
                curr_a0 = curr_a0 + incalphaloflap_L
        else:
            # Rigth wing
            if (mesh["ypos"][i] >= yflap[2]) and (mesh["ypos"][i] <= yflap[3]):
                # Flap delta alpha0
                curr_a0 = curr_a0 + incalphaloflap_R
        dx = np.sin(-curr_a0)
        dy = -(mesh["coord"][2, i+1]-mesh["coord"][2, i])/mesh["len"][i]
        dz = np.cos(-curr_a0)
        aux = 1/np.sqrt(dx*dx+dy*dy+dz*dz)
        normals[0, i] = dx*aux
        normals[1, i] = dy*aux
        normals[2, i] = dz*aux

    # Equations and Biot - Savart law
    A = np.zeros((mesh["N"], mesh["N"]))
    B = np.zeros((mesh["N"], 1))
    b = model["canopy"]["span"]
    one = 1.0
    xa = np.zeros(3)
    xb = np.zeros(3)
    xc = np.zeros(3)
    xd = np.zeros(3)
    for j in range(mesh["N"]):
        xp = mesh["xctrl"][:, j]
        nunit = normals[:, j]
        for k in range(mesh["N"]):
            xb[0] = mesh["xbound"][0, k]
            xa[0] = xb[0] + 20*b
            xc[0] = xb[0]
            xd[0] = xa[0]

            xb[1] = mesh["coord"][1, k]
            xa[1] = xb[1]
            xc[1] = mesh["coord"][1, k+1]
            xd[1] = xc[1]

            xb[2] = (mesh["xbound"][2, k] -
                     (mesh["coord"][2, k+1]-mesh["coord"][2, k])/2)
            xa[2] = xb[2]
            xc[2] = (mesh["xbound"][2, k] +
                     (mesh["coord"][2, k+1]-mesh["coord"][2, k])/2)
            xd[2] = xc[2]
            xa = xa[:]
            xb = xb[:]
            xc = xc[:]
            xd = xd[:]

            # First trailing vortex
            uind1 = vortxl(xa, xb, xp, one)
            # Bounded vortex
            uind2 = vortxl(xb, xc, xp, one)
            # Second trailing vortex
            uind3 = vortxl(xc, xd, xp, one)

            uindt = uind1+uind2+uind3
            A[j, k] = np.dot(uindt, nunit)
        uinf = -mesh["Vinf"][j, :]
        B[j, 0] = -np.dot(uinf, nunit)
    # Circulation
    Circ = np.linalg.solve(A, B)
    Circ = np.squeeze(Circ)

    # Loads computation (Kutta-Joukowsky)
    local_force = np.zeros((3, mesh["N"]))
    xa = np.zeros(3)
    xb = np.zeros(3)
    xc = np.zeros(3)
    xd = np.zeros(3)
    for j in range(mesh["N"]):
        xp = mesh["xbound"][:, j]
        Wk = np.zeros(3)
        for k in range(mesh["N"]):
            xb[0] = mesh["xbound"][0, k]
            xa[0] = xb[0] + 20*b
            xc[0] = xb[0]
            xd[0] = xa[0]

            xb[1] = mesh["coord"][1, k]
            xa[1] = xb[1]
            xc[1] = mesh["coord"][1, k+1]
            xd[1] = xc[1]

            xb[2] = (mesh["xbound"][2, k] -
                     (mesh["coord"][2, k+1]-mesh["coord"][2, k])/2)
            xa[2] = xb[2]
            xc[2] = (mesh["xbound"][2, k] +
                     (mesh["coord"][2, k+1]-mesh["coord"][2, k])/2)
            xd[2] = xc[2]
            xa = xa[:]
            xb = xb[:]
            xc = xc[:]
            xd = xd[:]

            # First trailing vortex
            wkind1 = vortxl(xa, xb, xp, Circ[k])
            # Second trailing vortex
            wkind3 = vortxl(xc, xd, xp, Circ[k])

            # Down Wash
            Wk = Wk + wkind1 + wkind3
        uinf = -mesh["Vinf2"][j, :]
        # Induced + Kinematic velocity
        v_i = Wk + uinf
        # Bound vortex
        xb = mesh["coord"][:, j]
        xc = mesh["coord"][:, j+1]
        # Gamma is per unit lenght
        d_g = (xc-xb)*Circ[j]
        # KJ
        local_force[:, j] = np.cross(v_i, d_g)*state["Rho"]

    #    Aerodynamic forces and moments according to HVM
    tot_force = np.zeros(6)
    for i in range(mesh["N"]):
        r = mesh["xbound"][:, i]
        # forces
        tot_force[:3] = tot_force[:3] + local_force[:, i]
        # moments MODIFICADO
        tot_force[3:] = tot_force[3:] + np.cross(r, local_force[:, i])
        if delta0_f[0] == delta0_f[1]:
            tot_force[1] = 0
            tot_force[3] = 0
            tot_force[5] = 0

    # Airfoil polar drag along the span
    CDp = 0.0
    Uref_m = np.squeeze(state["Uref_m"])
    qq = 0.5*state["Rho"]*np.dot(Uref_m, Uref_m)
    # Integration of profile drag (from the airfoil's polar)
    for i in range(mesh["N"]):
        # Assume Cl aprox Cz
        cz = local_force[2, i] / (qq * mesh["s"][i])
        cdp = (model["canopy"]["A_p"]*cz**2 +
               model["canopy"]["B_p"]*cz+model["canopy"]["C_p"])
        CDp = CDp + cdp*mesh["s"][i]

    CDp = CDp/model["canopy"]["Sref"]
    Vrefu = Uref_m / np.sqrt(np.dot(Uref_m, Uref_m))
    tot_force[:3] = tot_force[:3] + qq*model["canopy"]["Sref"]*CDp*Vrefu
    # tot_force = np.round(tot_force,10);

    # Total aerodynamic forces and moments of the canopy
    aero_force_canopy = tot_force[:3]
    aero_moment_canopy = tot_force[3:]
    b2c = state["b2c"]
    m2c = state["m2c"]
    aerodynamic_force = b2c.transpose().dot(m2c.dot(aero_force_canopy))
    state["aerodynamic_force"] = aerodynamic_force
    aerodynamic_moment = b2c.transpose().dot(m2c.dot(aero_moment_canopy))
    aerodynamic_moment = (aerodynamic_moment +
                          np.cross(np.squeeze(state["Xcg_aero"]), aerodynamic_force))
    state["aerodynamic_moment"] = aerodynamic_moment