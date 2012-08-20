CC=gcc

CFLAGS=-lm -I. -Imatrix -Ifem -Imesh -Imaterial-data/freezing -Wall -ggdb -O0

OBJECTFILES=integrate.o basis.o matrix.o mesh2d.o mtxsolver.o vector.o finite-element.o isoparam.o finite-element1d.o mesh1d.o solution.o freezing.o 

all: heat-explicit heat-implicit

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

heat-implicit: heat-implicit.o $(OBJECTFILES)
	$(CC) -o heat-implicit heat-implicit.o $(OBJECTFILES) $(CFLAGS)

meshtest: problems/meshtest.c $(OBJECTFILES)
	$(CC) -o meshtest problems/meshtest.c $(OBJECTFILES) $(CFLAGS)

clean:
	rm -rf spheroid 2dlaplace ce675p1 ce675p2 heat-implicit heat-explicit meshtest *.o *~


freezing.o: material-data/freezing/freezing.c material-data/freezing/freezing.h
	$(CC) -c material-data/freezing/freezing.c $(CFLAGS)

ce675p1.o: problems/ce675p1.c
	$(CC) -c problems/ce675p1.c $(CFLAGS)

ce675p2.o: problems/ce675p2.c
	$(CC) -c problems/ce675p2.c $(CFLAGS)

2dlaplace.o: problems/2dlaplace.c
	$(CC) -c problems/2dlaplace.c $(CFLAGS)

spheroid.o: problems/spheroid.c
	$(CC) -c problems/spheroid.c $(CFLAGS)

heat-implicit.o: problems/heat-implicit.c
	$(CC) -c problems/heat-implicit.c $(CFLAGS)

heat-explicit.o: problems/heat-explicit.c
	$(CC) -c problems/heat-explicit.c $(CFLAGS)


matrix.o: matrix/matrix.c matrix/matrix.h
	$(CC) -c matrix/matrix.c $(CFLAGS)

vector.o: matrix/vector.c matrix/matrix.h
	$(CC) -c matrix/vector.c $(CFLAGS)


mesh2d.o: mesh/mesh2d.c mesh/mesh2d.h
	$(CC) -c mesh/mesh2d.c $(CFLAGS)

mesh1d.o: mesh/mesh1d.h mesh/mesh1d.c
	$(CC) -c mesh/mesh1d.c $(CFLAGS)


finite-element.o: fem/finite-element.c fem/finite-element.h
	$(CC) -c fem/finite-element.c $(CFLAGS)

finite-element1d.o: fem/finite-element1d.c fem/finite-element1d.h
	$(CC) -c fem/finite-element1d.c $(CFLAGS)


integrate.o: integrate.c integrate.h
	$(CC) -c integrate.c $(CFLAGS)

basis.o: basis.c basis.h
	$(CC) -c basis.c $(CFLAGS)

mtxsolver.o: mtxsolver.c mtxsolver.h
	$(CC) -c mtxsolver.c $(CFLAGS)

solution.o: fem/solution.c fem/solution.h
	$(CC) -c fem/solution.c $(CFLAGS)

isoparam.o: isoparam.c isoparam.h
	$(CC) -c isoparam.c $(CFLAGS)

