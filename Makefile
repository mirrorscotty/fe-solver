CC=gcc

CFLAGS=-lm -I. -Imatrix -Ifem -Imesh -Wall -ggdb -O0

OBJECTFILES=integrate.o basis.o matrix.o mesh2d.o mtxsolver.o vector.o finite-element.o isoparam.o finite-element1d.o mesh1d.o solution.o time.o

all: heat

2dlaplace: 2dlaplace.o $(OBJECTFILES)
	gcc -o 2dlaplace 2dlaplace.o $(OBJECTFILES) $(CFLAGS)

ce675p1: ce675p1.o $(OBJECTFILES)
	$(CC) -o ce675p1 ce675p1.o $(OBJECTFILES) $(CFLAGS)

spheroid: spheroid.o $(OBJECTFILES)
	$(CC) -o spheroid spheroid.o $(OBJECTFILES) $(CFLAGS)

ce675p2: ce675p2.o $(OBJECTFILES)
	$(CC) -o ce675p2 ce675p2.o $(OBJECTFILES) $(CFLAGS)

heat: heat.o $(OBJECTFILES)
	$(CC) -o heat heat.o $(OBJECTFILES) $(CFLAGS)

clean:
	rm -rf spheroid 2dlaplace ce675p1 ce675p2 heat *.o *~



ce675p1.o: problems/ce675p1.c
	$(CC) -c problems/ce675p1.c $(CFLAGS)

ce675p2.o: problems/ce675p2.c
	$(CC) -c problems/ce675p2.c $(CFLAGS)

2dlaplace.o: problems/2dlaplace.c
	$(CC) -c problems/2dlaplace.c $(CFLAGS)

spheroid.o: problems/spheroid.c
	$(CC) -c problems/spheroid.c $(CFLAGS)

heat.o: problems/heat.c
	$(CC) -c problems/heat.c $(CFLAGS)


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

time.o: time.c time.h
	$(CC) -c time.c $(CFLAGS)
