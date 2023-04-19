# aero:
# 	gcc -shared -o simple_add.dll simple_add.c

aero:
	gcc -shared -o parafoil/aerodynamic/aero.dll -fPIC parafoil/aerodynamic/aero.c
	gcc -shared -o parafoil/aerodynamic/aero.so -fPIC parafoil/aerodynamic/aero.c
