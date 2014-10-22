CC=gcc
CFLAGS=-lm -I. -Imatrix -Isolver -Isolver/mesh -Isolver/integration -Isolver/ode -Imaterial-data/choi-okos -Iscaling -Igui/heating -Wall -g3 -O0
VPATH=problems gui solver/mesh solver/ode solver/integration matrix material-data scaling solver

OBJECTFILES=integrate.o basis.o mesh2d.o finite-element.o isoparam.o finite-element1d.o mesh1d.o solution.o auxsoln.o scaling_ht.o solve.o material-data.a matrix.a

all: heat-cyl

2dlaplace: 2dlaplace.o $(OBJECTFILES)
	gcc -o 2dlaplace 2dlaplace.o $(OBJECTFILES) $(CFLAGS)

ce675p1: ce675p1.o $(OBJECTFILES)
	$(CC) -o ce675p1 ce675p1.o $(OBJECTFILES) $(CFLAGS)

spheroid: spheroid.o $(OBJECTFILES)
	$(CC) -o spheroid spheroid.o $(OBJECTFILES) $(CFLAGS)

ce675p2: ce675p2.o $(OBJECTFILES)
	$(CC) -o ce675p2 ce675p2.o $(OBJECTFILES) $(CFLAGS)

heat-explicit: heat-explicit.o $(OBJECTFILES)
	$(CC) -o heat-explicit heat-explicit.o $(OBJECTFILES) $(CFLAGS)

heat-cyl: heat-cyl.o heat-gui.o $(OBJECTFILES)
	$(CC) -o heat-cyl heat-cyl.o heat-gui.o $(OBJECTFILES) $(CFLAGS)

clean:
	rm -rf spheroid 2dlaplace ce675p1 ce675p2 heat-explicit heat-cyl meshtest *.o *~
	$(MAKE) -C matrix clean
	$(MAKE) -C material-data clean


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

# Mesh-related files
mesh2d.o: mesh/mesh2d.c mesh/mesh2d.h
mesh1d.o: mesh/mesh1d.h mesh/mesh1d.c

# Finite element-specific files
finite-element.o: finite-element.c finite-element.h
finite-element1d.o: finite-element1d.c finite-element1d.h
solution.o: solution.c solution.h
auxsoln.o: auxsoln.c auxsoln.h

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

