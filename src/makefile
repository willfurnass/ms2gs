all:	kind.o ALliball.o ebv.o aux_sub11.o ms2gs.o
	gfortran -O3 kind.o ALliball.o ebv.o aux_sub11.o ms2gs.o -o ms2gs
clean:
	rm *.mod *.o
kind.o:	kind.f90
	gfortran -O3 -c kind.f90
ALliball.o:	ALliball.f90 kind.f90
	gfortran -O3 -c ALliball.f90
ebv.o:	ebv.f90
	gfortran -O3 -c ebv.f90
aux_sub11.o:	aux_sub11.f90
	gfortran -O3 -c aux_sub11.f90
ms2gs.o:	ms2gs.f90
	gfortran -O3 -c ms2gs.f90
install:	ms2gs README.ms2gs
	cp ms2gs /opt/.
	cp README.ms2gs /opt/.
# 	gfortran -fcheck=all kind.o ALliball.o ebv.o aux_sub11.o ms2gs.o -o ms2gs
