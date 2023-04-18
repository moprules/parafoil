"""
Модуль обёртка для загрузки библиотек на си
и задания параметров функций написанных на си
"""


from ..cliber import load_lib
import ctypes
from numpy.ctypeslib import ndpointer

# Тип данных для трансляции двумерных массивов numpy в массивы на си
с_arr2d = ndpointer(dtype=ctypes.c_double, ndim=2, flags='C_CONTIGUOUS')

# Загружаем библиотеку aero - расчёт аэродинамики в СИ
aero = load_lib("aero", __file__)

# Назначаем параметры функции расчёта по закону Биота-Саварта
# Функция ничего не возвращает
aero.biot_savart.restype = None
aero.biot_savart.argtypes = [ctypes.c_int,      # N
                             ctypes.c_double,   # b
                             с_arr2d,           # A
                             с_arr2d,           # B
                             с_arr2d,           # normals
                             с_arr2d,           # xctrl
                             с_arr2d,           # xbound
                             с_arr2d,           # coord
                             с_arr2d]           # Vinf

# Назначаем параметры функции расчёта Kutta-Joukowsky
# Функция ничего не возвращает
aero.kutta_joukowsky.restype = None
aero.kutta_joukowsky.argtypes = [ctypes.c_int,      # N
                                 ctypes.c_double,   # Rho
                                 с_arr2d,           # local_force
                                 с_arr2d,           # xbound
                                 с_arr2d,           # coord
                                 с_arr2d]           # Circ

__all__ = ["aero"]
