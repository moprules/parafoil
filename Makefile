ifeq ($(OS),Windows_NT)
    LIB_EXT = dll
else
    LIB_EXT = so
endif

aero:
	gcc -shared -Ofast -o parafoil/aerodynamic/aero.$(LIB_EXT) -fPIC parafoil/aerodynamic/aero.c
