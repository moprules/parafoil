"""
Модуль обёртка для загрузки библиотек на си
и задания параметров функций написанных на си
"""


from ..cliber import load_lib
import ctypes
from numpy.ctypeslib import ndpointer

# Тип данных для трансляции двумерных массивов numpy в массивы на си
c_arr2d = ndpointer(dtype=ctypes.c_double, ndim=2, flags='C_CONTIGUOUS')
c_arr = ndpointer(dtype=ctypes.c_double, ndim=1, flags='C_CONTIGUOUS')

# Загружаем библиотеку aero - расчёт аэродинамики в СИ
aero = load_lib("aero", __file__)


# Описание параметров для функции расчёта нормалей
aero.set_normals.restype = None
aero.set_normals.argtypes = [ctypes.c_int,      # N
                             ctypes.c_double,   # angh
                             c_arr2d,           # normals
                             c_arr,           # delta0_f
                             c_arr,           # yflap
                             c_arr,           # aphalo
                             c_arr2d,           # xbound
                             c_arr,           # ypos
                             c_arr2d,           # coord
                             c_arr]           # len

# Назначаем параметры функции расчёта по закону Биота-Саварта
# Функция ничего не возвращает
aero.biot_savart.restype = None
aero.biot_savart.argtypes = [ctypes.c_int,      # N
                             ctypes.c_double,   # b
                             c_arr2d,           # A
                             c_arr2d,           # B
                             c_arr2d,           # normals
                             c_arr2d,           # xctrl
                             c_arr2d,           # xbound
                             c_arr2d,           # coord
                             c_arr2d]           # Vinf

# Назначаем параметры функции расчёта Kutta-Joukowsky
# Функция ничего не возвращает
aero.kutta_joukowsky.restype = None
aero.kutta_joukowsky.argtypes = [ctypes.c_int,      # N
                                 ctypes.c_double,   # b
                                 ctypes.c_double,   # Rho
                                 c_arr2d,           # local_force
                                 c_arr2d,           # xbound
                                 c_arr2d,           # coord
                                 c_arr2d,           # Circ
                                 c_arr2d]           # Vinf2

__all__ = ["aero"]
