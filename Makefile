CC=gcc
CFLAGS=-lm -I. -Imatrix -Imaterial-data -Isolver -Ioutput -Isolver/mesh -Isolver/integration -Isolver/ode -Iscaling -Igui/heating -Wall -g3 -O0 
VPATH=problems problems/slab-drying gui solver/mesh solver/ode solver/integration matrix material-data scaling solver output

OBJECTFILES=integrate.o basis.o mesh2d.o finite-element.o isoparam.o finite-element1d.o mesh1d.o solution.o auxsoln.o scaling_ht.o linsolve1d.o nlinsolve1d.o predict1d.o linsolve2d.o nlinsolve2d.o output.o material-data.a matrix.a

all: diffusion

doc:
	doxygen DoxyFile
	make -C doc/latex
	cp doc/latex/refman.pdf doc/Reference.pdf

2dlaplace: 2dlaplace.o $(OBJECTFILES)
	$(CC) -o $@ $^ $(CFLAGS)

ce675p1: ce675p1.o $(OBJECTFILES)
	$(CC) -o $@ $^ $(CFLAGS)

spheroid: spheroid.o $(OBJECTFILES)
	$(CC) -o $@ $^ $(CFLAGS)

ce675p2: ce675p2.o $(OBJECTFILES)
	$(CC) -o $@ $^ $(CFLAGS)

heat-explicit: heat-explicit.o $(OBJECTFILES)
	$(CC) -o $@ $^ $(CFLAGS)

heat-cyl: heat-cyl.o heat-gui.o $(OBJECTFILES)
	$(CC) -o $@ $^ $(CFLAGS)

heat-transfer: heat-transfer.o ht-main.o common.o $(OBJECTFILES)
	$(CC) -o $@ $^ $(CFLAGS)

ht-mt: diffusion.o heat-transfer.o main.o common.o $(OBJECTFILES)
	$(CC) -o $@ $^ $(CFLAGS)

diffusion: diffusion.o mt-main.o common.o $(OBJECTFILES)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
#rm -rf spheroid 2dlaplace ce675p1 ce675p2 heat-explicit heat-cyl meshtest
	rm -rf heat-transfer
	rm -rf *.o *.a
	$(MAKE) -C matrix clean
	$(MAKE) -C material-data clean
	rm -rf doc


freezing.o: material-data/freezing/freezing.c material-data/freezing/freezing.h
	$(CC) -c material-data/freezing/freezing.c $(CFLAGS)

heat-implicit.o: problems/heat-implicit.c
	$(CC) -c problems/heat-implicit.c $(CFLAGS)

heat-gui.o: heating/heat-gui.c heating/heat-gui.h
	$(CC) -c gui/heating/heat-gui.c $(CFLAGS)

2dlaplace.o: 2dlaplace.c
ce675p1.o: ce675p1.c
ce675p2.o: ce675p2.c
heat-cyl.o: problems/heat-cyl.c gui/heating/heat-gui.h
heat-explicit.o: heat-explicit.c
spheroid.o: spheroid.c
heat-transfer.o: heat-transfer.h
main.o: heat-transfer.h diffusion.h common.h
mt-main.o: diffusion.h common.h
ht-main.o: heat-transfer.h common.h
common.o: common.h
diffusion.o: diffusion.h

# Mesh-related files
mesh2d.o: mesh/mesh2d.c mesh/mesh2d.h
mesh1d.o: mesh/mesh1d.h mesh/mesh1d.c

# Finite element-specific files
finite-element.o: finite-element.c finite-element.h
finite-element1d.o: finite-element1d.c finite-element1d.h
solution.o: solution.c solution.h
auxsoln.o: auxsoln.c auxsoln.h

# Solvers
linsolve1d.o: solve.h
linsolve2d.o: solve.h
nlinsolve1d.o: solve.h
nlinsolve2d.o: solve.h
predict1d.o: solve.h

# Scaling stuff
scaling_ht.o: scaling_ht.h
scaling_pasta.o: scaling_pasta.h

# Random files that don't have a home
about.o: about.h
basis.o: basis.h
isoparam.o: isoparam.h
integrate.o: integrate.c integrate.h
output.o: output.h

matrix.a:
	$(MAKE) -C matrix
	cp matrix/matrix.a .

material-data.a:
	$(MAKE) -C material-data
	cp material-data/material-data.a .

