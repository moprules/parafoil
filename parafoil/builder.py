import numpy as np
from . import matrix


def build(params: dict) -> dict:
    """
    На основе исходных параметров строим удобную модель для расчётов
    переводим градусы в радианы и добавляем параметры по умолчанию,
    которые не нужно указать в файле параметров.
    """
    # Копируем параметры в новый словарь
    model = params.copy()
    # Подготавливаем модель для расчётов по каждой секции
    prepare_body(model)
    prepare_canopy(model)
    prepare_wind(model)
    prepare_autopilot(model)
    prepare_target(model)
    return model


def prepare_body(model: dict):
    # Параметры твёрдого тела
    body = model["body"]

    # переводим начальные углы ориентации в радианы
    # Угол крена
    body["roll"] = np.radians(body["roll"])
    # Угол тангажа
    body["pitch"] = np.radians(body["pitch"])
    # Угол рысканья
    body["yaw"] = np.radians(body["yaw"])
    # Линейные скорости центра гравитации преобразуем в вектор столбец
    lv = body["linear_velocity"]
    body["linear_velocity"] = np.asarray([[lv["u"]], [lv["v"]], [lv["w"]]])

    # Угловые скорости преобразуем в вектор столбец
    av = body["angular_velocity"]
    body["angular_velocity"] = np.asarray([[av["p"]], [av["q"]], [av["r"]]])
    # также приводим градусы в радианы
    body["angular_velocity"] = np.radians(body["angular_velocity"])


def prepare_canopy(model: dict):
    # Параметры парашюта
    canopy = model["canopy"]

    # Переводим углы в радианы
    canopy["sweep"] = np.radians(canopy["sweep"])
    canopy["dihedral"] = np.radians(canopy["dihedral"])
    canopy["rigging_angle"] = np.radians(canopy["rigging_angle"])

    canopy["Xref_m"] = np.asarray([[0.25*canopy["root_chord"]], [0], [0]])

    # преобразуем координаты передней кромки корневой хорды в вектор столбец
    le = canopy["leading_edge"]
    canopy["leading_edge"] = np.asarray([[le["x"]], [le["y"]], [le["z"]]])
    canopy["Xref_b"] = canopy["leading_edge"].copy()

    # преобразуем координаты полезной нагрузки в вектор столбец
    pc = canopy["paylod_coord"]
    canopy["paylod_coord"] = np.asarray([[pc["x"]], [pc["y"]], [pc["z"]]])
    canopy["Xcg_payload"] = canopy["paylod_coord"].copy()

    # преобразуем координаты центра кажущейся массы в вектор столбец
    amc = canopy["apparent_mass_center"]
    canopy["apparent_mass_center"] = np.asarray(
        [[amc["x"]], [amc["y"]], [amc["z"]]])
    canopy["Xam"] = canopy["apparent_mass_center"].copy()

    # Преобразуем начальные отклонения управляющей поверхности в радианы
    deflection = canopy["cs"]["deflection"]
    deflection["R"] = np.radians(deflection["R"])
    deflection["L"] = np.radians(deflection["L"])

    # угор шарнира приводм в радианы
    canopy["cs"]["angh"] = np.radians(canopy["cs"]["angh"])

    # соотношение сторон
    canopy["AR"] = (canopy["span"] ** 2) / canopy["Sref"]
    # крайняя хорда парашюта
    canopy["tip_chord"] = canopy["root_chord"] * canopy["taper"]

    # Sweep angle LE
    canopy["AngLE"] = canopy["sweep"]
    # sweep angle 1/4 chord
    canopy["AngC"] = np.arctan(np.tan(
        canopy["AngLE"])-2*(1/4)*(canopy["root_chord"]-canopy["tip_chord"])/canopy["span"])

    canopy["a0"]["root"] = np.radians(canopy["a0"]["root"])
    canopy["a0"]["tip"] = np.radians(canopy["a0"]["tip"])


def prepare_wind(model: dict):
    # Преобразуем скорост ветра в вектор столбец
    w = model["wind"]
    model["wind"] = np.asarray([[w["x"]], [w["y"]], [w["z"]]])


def prepare_autopilot(model: dict):
    ap = model["autopilot"]
    dt = model["time"]["step"]

    ap["Ki_gs"] *= dt
    ap["Kd_gs"] /= dt
    ap["Ki"] *= dt
    ap["Kd"] /= dt

    ap["incidence_lower_limit"] = np.radians(ap["incidence_lower_limit"])
    ap["incidence_upper_limit"] = np.radians(ap["incidence_upper_limit"])

    ap["lasterr"] = 0
    ap["errsum"] = 0
    ap["heading_sumerr"] = 0
    ap["heading_lasterr"] = 0
    ap["loiter_decision_altitude"] = 100
    ap["left_brake_command"] = 0
    ap["right_brake_command"] = 0
    ap["glide_slope_command"] = 0

    ap["tol_angle"] = np.radians(ap["tol_angle"])


def prepare_target(model: dict):
    target = model["target"]
    body = model["body"]

    pos_north = (target["north"] - body["north"])
    pos_east = (target["east"] - body["east"])
    target["init_dist"] = (pos_north**2+pos_east**2)**0.5
