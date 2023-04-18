import os
import sys
import ctypes

# Логическая переменная, запускается код на Windows или нет
IS_WIN = "win" in sys.platform.lower()


def load_lib(lib_name: str, path_to_lib: str):
    ext = "so"
    if IS_WIN:
        ext = "dll"

    # Названеи библиотеки с расширением
    lib_name = f"{lib_name}.{ext}"
    # На всякий случай путь до библиотеки должен быть абсолютным
    path_to_lib = os.path.abspath(path_to_lib)
    # если вместо директории предали переменную __file__
    if os.path.isfile(path_to_lib):
        # Обрубаем путь до папки без названия файла
        path_to_lib = os.path.dirname(path_to_lib)

    # полный путь до библиотеки
    lib_path = os.path.join(path_to_lib, lib_name)
    # Загружаем библиотеку в код на питоне
    lib = ctypes.CDLL(lib_path)

    # Возвращаем загруженную бибилиотеку в основную программу
    return lib
