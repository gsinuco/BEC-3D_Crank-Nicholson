CC = mpic++
GF = mpif90
GFS= gfortran
GFFLAGS = -L/usr/lib64/atlas/lib -llapack -lf77blas -lcblas -latlas -lgsl -lgslcblas -O3 
CFLAGS =  -lgsl -lgslcblas -O3  


all: BEC-ground  BEC-TE BEC-gs_serie BEC-te_serie

clean:
	rm  BEC-ground  BEC-TE BEC-gs_serie BEC-te_serie *.mod

BEC-ground: CN_paralelo_ground.f90  
	$(GF)  -o $@ CN_paralelo_ground.f90 $(GFFLAGS) 

BEC-TE: CN_paralelo_time_evol_test.f90  
	$(GF)  -o $@ CN_paralelo_time_evol_test.f90 $(GFFLAGS) 

BEC-gs_serie:  BEC-Crank-Nicolson-3D-groundV2.f90
	$(GFS)  -o $@ BEC-Crank-Nicolson-3D-groundV2.f90 $(GFFLAGS) 

BEC-te_serie:  BEC-Crank-Nicolson-3D-timeEvol.f90
	$(GFS)  -o $@ BEC-Crank-Nicolson-3D-timeEvol.f90 $(GFFLAGS) 
