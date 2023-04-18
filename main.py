# from parafoil import PFSim


# import ctypes
# import numpy as np
# from numpy.ctypeslib import ndpointer

# # Определение типа элементов массива
# dtype = np.double

# # Загрузка библиотеки
# my_lib = ctypes.cdll('./lol.so')

# # Определение входного и выходного типов функции
# test = my_lib.test
# test.argtypes = [ndpointer(
#     dtype=dtype, ndim=2, flags='C_CONTIGUOUS'), ctypes.c_int, ctypes.c_int]
# test.restype = None

# # Создание и передача массива в функцию
# m = [[1, 2], [3, 4]]
# my_array = np.asarray(m)
# my_array = np.ascontiguousarray(my_array, dtype=dtype)
# test(my_array, 2, 2)
# print(my_array)

from parafoil import PFSim


def main():
    lander = PFSim("data/space_rider.yaml")
    lander.start()


if __name__ == "__main__":
    main()
