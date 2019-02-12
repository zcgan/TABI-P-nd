# makefile for bimpb
#F90= ifort
#flag= -fast -c  
F90 = gfortran
flag= -O -c 
bimpb.out: var_modules.o treecode3d_pb.o main.o readin.o parameter.o solvers.o quadrature.o kirkwood.o dgmres_dep.o
	$(F90)  -o bimpb.out *.o 
var_modules.o:	var_modules.f90
	$(F90) $(flag) var_modules.f90
main.o:		main.f90
	$(F90) $(flag) main.f90
readin.o:	readin.f90
	$(F90) $(flag) readin.f90
parameter.o:	parameter.f90
	$(F90) $(flag) parameter.f90
solvers.o:	solvers.f90
	$(F90) $(flag) solvers.f90
quadrature.o:	quadrature.f90
	$(F90) $(flag) quadrature.f90
kirkwood.o:     kirkwood.f90
	$(F90) $(flag) kirkwood.f90
treecode3d_pb.o:	treecode3d_pb.f
	$(F90) $(flag) treecode3d_pb.f
dgmres_dep.o:	dgmres_dep.f
	$(F90) $(flag) dgmres_dep.f
molecule.mod:   var_modules.f90
	$(F90) $(flag) var_modules.f90
comdata.mod:   var_modules.f90
	$(F90) $(flag) var_modules.f90
bicg.mod:   	var_modules.f90
	$(F90) $(flag) var_modules.f90
treecode.mod:	var_modules.f90
	$(F90) $(flag) var_modules.f90
treecode3d_procedures.mod:	treecode3d_pb.f
	$(F90) $(flag) treecode3d_pb.f
clean: 
	rm *.o *.mod bimpb.out
