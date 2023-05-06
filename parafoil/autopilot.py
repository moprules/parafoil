from . import matrix
import numpy as np


def control(model: dict, state: dict):

    # состояние автопилота
    ap = state["ap"]

    # position_target_north = (model["target"]["north"] - state["pos_north"])
    # position_target_east = (model["target"]["east"] - state["pos_east"])

    position_target_north = model["target"]["init_north"]
    position_target_east = model["target"]["init_east"]

    i2t = matrix.inertia_track(position_target_north, position_target_east)

    # From body to inertia axes
    VAUX2 = state["i2b"].transpose().dot(state["Vcg_b"])
    # From inertia to track axes
    VAUX1 = i2t.dot(VAUX2)
    velocity_x_track = VAUX1[0, 0]
    velocity_y_track = VAUX1[1, 0]
    VAUX1 = i2t.dot(state["Vwind"])
    wind_x_track = VAUX1[0, 0]
    wind_y_track = VAUX1[1, 0]
    x_dot_track = velocity_x_track + wind_x_track
    y_dot_track = velocity_y_track + wind_y_track

    # ALONG AND CROSS-TRACK DISTANCE CALCULATION
    VAUX2 = np.zeros((3, 1))
    VAUX2[0, 0] = position_target_north - state["pos_north"]
    VAUX2[1, 0] = position_target_east - state["pos_east"]
    VAUX2[2, 0] = 0
    VAUX1 = i2t.dot(VAUX2)
    x_track = VAUX1[0, 0]
    y_track = VAUX1[1, 0]

    # ALTITUDE CONTROL
    distance_to_target = np.sqrt(x_track**2+y_track**2)
    auxiliar = abs((state["Vcg_b"][0, 0]*x_track +
                    state["Vcg_b"][1, 0]*y_track) /
                   (distance_to_target*np.linalg.norm([state["Vcg_b"][0, 0], state["Vcg_b"][1, 0]])))
    alpha = 0.1*np.arccos(auxiliar)
    distance_to_target2 = distance_to_target/(1-alpha)
    distance_to_target_ratio = distance_to_target/model["target"]["init_dist"]

    # JUST GLIDE WITHOUT ALTITUDE CONTROL BEFORE REACHING 50% OF TARGET DISTANCE
    if state["dlost_altitude"] == 0:
        glide_ratio = -0
    else:
        glide_ratio = -state["dtrav_dist"]/state["dlost_altitude"]
    glide_ratio_ideal = distance_to_target2/(state["altitude"]+10)
    err = glide_ratio - glide_ratio_ideal
    derr = err - ap["lasterr"]
    ap["lasterr"] = err
    if distance_to_target_ratio > ap["altitude_control_activation"]:
        if state["dlost_altitude"] == 0:
            glide_ratio = -0
        else:
            glide_ratio = -state["dtrav_dist"]/state["dlost_altitude"]
        ap["loiter_decision_altitude"] = distance_to_target/glide_ratio
        auto_pilot = True
    else:
        # IF TOO HIGH -> JUST TURNS, ELSE -> FLY WITH GLIDESLOPE CONTROL
        if state["altitude"] > ap["loiter_decision_altitude"]*ap["loiter_decision_offset"]:
            # TURN AND TURN
            ap["left_brake_command"] = -ap["auto_turn_brake_input"]
            auto_pilot = False
        else:
            # GLIDESLOPE CONTROL
            ap["errsum"] = ap["errsum"] + err
            ap["glide_slope_command"] = (ap["Kp_gs"]*err +
                                         ap["Ki_gs"]*ap["errsum"] +
                                         ap["Kd_gs"]*derr)
            if ap["glide_slope_command"] < ap["incidence_lower_limit"]:
                ap["glide_slope_command"] = ap["incidence_lower_limit"]
            elif ap["glide_slope_command"] > ap["incidence_upper_limit"]:
                ap["glide_slope_command"] = ap["incidence_upper_limit"]
            auto_pilot = True
            state["ctrl"] = True

    # ap["k_parameter"]*x_track/x_dot_track - y_track/y_dot_track;
    heading_error = ap["k_parameter"]*x_track*y_dot_track - y_track*x_dot_track
    ap["heading_sumerr"] = ap["heading_sumerr"] + heading_error
    heading_derr = heading_error - ap["heading_lasterr"]
    if auto_pilot:
        desired_yaw_rate = (ap["Kp"]*heading_error +
                            ap["Ki"]*ap["heading_sumerr"] +
                            ap["Kd"]*heading_derr)
        yaw_rate_command = desired_yaw_rate - state["Wb"][2, 0]
        yaw_rate_sign = np.sign(yaw_rate_command)
        if abs(yaw_rate_command) > ap["max_yaw_rate"]:
            if state["ctrl"]:
                yaw_rate_command = (yaw_rate_sign *
                                    ap["auto_turn_brake_input"]*1.5)
            else:
                yaw_rate_command = yaw_rate_sign*abs(ap["max_yaw_rate"])
        if yaw_rate_command <= 0:
            ap["left_brake_command"] = ap["roll_gain"]*yaw_rate_command
        elif yaw_rate_command > 0:
            ap["right_brake_command"] = -ap["roll_gain"]*yaw_rate_command
    else:
        yaw_rate_command = 0

    ap["heading_lasterr"] = heading_error

    # FLIGHT CONTROLS SIMULATE
    # Flare
    if state["altitude"] < 15 or state["altitude"] > model["body"]["altitude"]:
        state["rigging_angle"] = np.radians(-7)
        # left brake input
        state["delta0_f"][0] = np.radians(-25)
        # right brake input
        state["delta0_f"][1] = np.radians(-25)
    else:
        state["delta0_f"][0] = ap["left_brake_command"]
        state["delta0_f"][1] = ap["right_brake_command"]
        state["rigging_angle"] = (state["rigging_angle"] +
                                  ap["incidence_gain"]*ap["glide_slope_command"])
        if state["rigging_angle"] > -np.radians(0):
            state["rigging_angle"] = -np.radians(0)
        if state["rigging_angle"] < -np.radians(30):
            state["rigging_angle"] = -np.radians(30)

    # New rotation matrix from body to canopy axes
    b2c = [[np.cos(state["rigging_angle"]), 0.0, -np.sin(state["rigging_angle"])],
           [0.0,                            1.0, 0.0],
           [np.sin(state["rigging_angle"]), 0.0, np.cos(state["rigging_angle"])]]
    state["b2c"] = np.asarray(b2c)
