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
    x1 = X1[0];     y1 = X1[1];     z1 = X1[2];
    x2 = X2[0];     y2 = X2[1];     z2 = X2[2];
    xp = XP[0];     yp = XP[1];     zp = XP[2];
    r0 = X2-X1;     r1 = XP-X1;     r2 = XP-X2;
# function [u] = vortxl (X1,X2,XP,gamma) 

#     %Velocity induced by a constant strenght linear vortex
#     x1 = X1(1);     y1 = X1(2);     z1 = X1(3);
#     x2 = X2(1);     y2 = X2(2);     z2 = X2(3);
#     xp = XP(1);     yp = XP(2);     zp = XP(3);
#     r0 = X2-X1;     r1 = XP-X1;     r2 = XP-X2;

#     norm_r1 = norm(r1);
#     norm_r2 = norm(r2);

#     r1xr2 = cross(r1,r2);
#     norm_r1xr2 = norm(r1xr2);

#     u = zeros(3,1);

#     inv_r1xr2 = 1.0 / norm_r1xr2;
#     inv_r1 = 1.0 / norm_r1;
#     inv_r2 = 1.0 / norm_r2;

#     a = r0 * inv_r1xr2;
#     b = r1*inv_r1 - r2*inv_r2;
#     c = dot(a,b);

#     u = gamma.*0.25/pi.*c.*r1xr2.*inv_r1xr2;
# end

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
            if -mesh["ypos"][i] <= yflap[1] and -mesh["ypos"][i] >= yflap[0]:
                # Flap delta alpha0
                curr_a0 = curr_a0 + incalphaloflap_L
        else:
            # Rigth wing
            if mesh["ypos"][i] >= yflap[2] and mesh["ypos"][i] <= yflap[3]:
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
    A = np.zeros(mesh["N"])
    Wk = np.zeros(mesh["N"])
    B = np.zeros((mesh["N"],1))
    one = 1.0
    
    # for j in range(mesh["N"]):
    #     xp = xctrl(:,j) ;
    #     nunit = normals(:,j);
    #     for k = 1:N
    #         xb(1) = xbound(1,k);    xa(1) = xb(1) + 20*b;
    #         xc(1) = xb(1);          xd(1) = xa(1);

    #         xb(2) = coord(2,k);     xa(2) = xb(2);
    #         xc(2) = coord(2,k+1);   xd(2) = xc(2);

    #         xb(3) = xbound(3,k)-(coord(3,k+1)-coord(3,k))/2;    xa(3) = xb(3);
    #         xc(3) = xbound(3,k)+(coord(3,k+1)-coord(3,k))/2;    xd(3) = xc(3);
    #         xa=xa(:);   xb=xb(:);   xc=xc(:);   xd=xd(:);

    #         uind1 = vortxl(xa,xb,xp,one);  % First trailing vortex
    #         uind2 = vortxl(xb,xc,xp,one);  % Bounded vortex
    #         uind3 = vortxl(xc,xd,xp,one);  % Second trailing vortex

    #         uindt = uind1+uind2+uind3;
    #         A(j,k) = dot(uindt',nunit);
    #     end
    #     uinf=-Vinf(j,:);
    #     B(j,1) = -dot(uinf',nunit);
    # end
    # % Circulation
    # Circ = A\B;


# function [aerodynamic_force,aerodynamic_moment] = HVM(N,rho,b,coord,xbound,xctrl,cg_pos,len,s,S,A_p,B_p,C_p,m2c,b2c,Uinf,ypos,yflap,angh,delta0_f,alphalo,Vinf,Vinf2)
#
#     % Normals at the control ponints considering the CS deflection
#     incalphaloflap_L = (delta0_f(1)/pi)*(pi-angh+sin(angh));
#     incalphaloflap_R = (delta0_f(2)/pi)*(pi-angh+sin(angh));
#     normals = zeros(3,N);
#     for i = 1:N
#         curr_a0 = alphalo(i);
#         if ( xbound(2,i) <= 0.0 ) % Left wing
#             if ( -ypos(i) <= yflap(2) && -ypos(i) >= yflap(1) )
#                 curr_a0 = curr_a0 + incalphaloflap_L;  % Flap delta alpha0
#             end
#         else % Rigth wing
#             if ( ypos(i) >= yflap(3) && ypos(i) <= yflap(4) )
#                 curr_a0 = curr_a0 + incalphaloflap_R;  % Flap delta alpha0
#             end
#         end
#         dx = sin(-curr_a0);
#         dy = -(coord(3,i+1)-coord(3,i))/len(i);
#         dz = cos(-curr_a0);
#         aux = 1/sqrt(dx*dx+dy*dy+dz*dz);
#         normals(1,i) = dx*aux;
#         normals(2,i) = dy*aux;
#         normals(3,i) = dz*aux;
#     end

#     % Equations and Biot - Savart law
#     A = zeros(N);   Wk = zeros(N);  B = zeros(N,1); one = 1.0;
#     for j = 1:N
#         xp = xctrl(:,j) ;
#         nunit = normals(:,j);
#         for k = 1:N
#             xb(1) = xbound(1,k);    xa(1) = xb(1) + 20*b;
#             xc(1) = xb(1);          xd(1) = xa(1);

#             xb(2) = coord(2,k);     xa(2) = xb(2);
#             xc(2) = coord(2,k+1);   xd(2) = xc(2);

#             xb(3) = xbound(3,k)-(coord(3,k+1)-coord(3,k))/2;    xa(3) = xb(3);
#             xc(3) = xbound(3,k)+(coord(3,k+1)-coord(3,k))/2;    xd(3) = xc(3);
#             xa=xa(:);   xb=xb(:);   xc=xc(:);   xd=xd(:);

#             uind1 = vortxl(xa,xb,xp,one);  % First trailing vortex
#             uind2 = vortxl(xb,xc,xp,one);  % Bounded vortex
#             uind3 = vortxl(xc,xd,xp,one);  % Second trailing vortex

#             uindt = uind1+uind2+uind3;
#             A(j,k) = dot(uindt',nunit);
#         end
#         uinf=-Vinf(j,:);
#         B(j,1) = -dot(uinf',nunit);
#     end
#     % Circulation
#     Circ = A\B;

#     % Loads computation (Kutta-Joukowsky)
#     local_force = zeros(3,N) ;
#     for j = 1:N
#         xp = xbound(:,j);
#         Wk = zeros(3,1);
#         for k = 1:N
#             xb(1) = xbound(1,k);    xa(1) = xb(1) + 20*b;
#             xc(1) = xb(1);          xd(1) = xa(1);

#             xb(2) = coord(2,k);     xa(2) = xb(2);
#             xc(2) = coord(2,k+1);   xd(2) = xc(2);

#             xb(3) = xbound(3,k)-(coord(3,k+1)-coord(3,k))/2;    xa(3) = xb(3);
#             xc(3) = xbound(3,k)+(coord(3,k+1)-coord(3,k))/2;    xd(3) = xc(3);
#             xa=xa(:);   xb=xb(:);   xc=xc(:);   xd=xd(:);

#             wkind1 = vortxl(xa,xb,xp,Circ(k));  % First trailing vortex
#             wkind3 = vortxl(xc,xd,xp,Circ(k));  % Second trailing vortex

#             Wk = Wk + wkind1 + wkind3; %Down Wash
#         end
#         uinf=-Vinf2(j,:);
#         v_i = Wk' + uinf;                      % Induced + Kinematic velocity
#         xb = coord(:,j); xc = coord(:,j+1);    % Bound vortex
#         d_g = (xc-xb)*Circ(j);                 % Gamma is per unit lenght
#         local_force(:,j) = cross(v_i,d_g).*rho;  % KJ
#     end

#     % Aerodynamic forces and moments according to HVM
#     tot_force = zeros(6,1);
#     for i = 1:N
#         r = xbound(:,i);
#         tot_force(1:3) = tot_force(1:3) + local_force(:,i) ;  % forces
#         tot_force(4:6) = tot_force(4:6) + cross(r,local_force(:,i));   % moments MODIFICADO
#         if delta0_f(1) == delta0_f(2)
#             tot_force(2)=0;
#             tot_force(4)=0;
#             tot_force(6)=0;
#         end
#     end

#     % Airfoil polar drag along the span
#     CDp = 0.0;
#     qq = 0.5*rho*dot(Uinf,Uinf);
#      for i = 1:N %Integration of profile drag (from the airfoil's polar)
#          cz = local_force(3,i) / (qq * s(i));  % Assume Cl aprox Cz
#          cdp = A_p*cz^2+B_p*cz+C_p;
#          CDp = CDp + cdp*s(i);
#      end
#     CDp = CDp/S;
#     Vrefu = Uinf / sqrt(dot(Uinf,Uinf));
#     tot_force(1:3) = tot_force(1:3) + qq*S*CDp*Vrefu';
#     tot_force(:) = round(tot_force(:),10);

#     % Total aerodynamic forces and moments of the canopy
#     aero_force_canopy = tot_force(1:3);
#     aero_moment_canopy = tot_force(4:6);
#     aerodynamic_force = b2c'*(m2c*aero_force_canopy);
#     aerodynamic_moment = b2c'*(m2c*aero_moment_canopy);
#     aerodynamic_moment = aerodynamic_moment + cross(cg_pos,aerodynamic_force);
# end
