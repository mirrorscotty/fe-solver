CC=gcc

CFLAGS=-lm -I. -I..-Wall -ggdb -O2

OBJECTFILES=integrate.o basis.o matrix.o mesh.o mtxsolver.o vector.o finite-element.o isoparam.o finite-element1d.o mesh1d.o

all: ce675p1 

2dlaplace: 2dlaplace.o $(OBJECTFILES)
	gcc -o 2dlaplace 2dlaplace.o $(OBJECTFILES) $(CFLAGS)

ce675p1: ce675p1.o $(OBJECTFILES)
	$(CC) -o ce675p1 ce675p1.o $(OBJECTFILES) $(CFLAGS)

spheroid: spheroid.o $(OBJECTFILES)
	$(CC) -o spheroid spheroid.o $(OBJECTFILES) $(CFLAGS)

ce675p2: ce675p2.o $(OBJECTFILES)
	$(CC) -o ce675p2 ce675p2.o $(OBJECTFILES) $(CFLAGS)

clean:
	rm -rf spheroid 2dlaplace ce675p1 ce675p2 *.o *~



ce675p1.o: problems/ce675p1.c
	$(CC) -c problems/ce675p1.c $(CFLAGS)

ce675p2.o: problems/ce675p2.c
	$(CC) -c problems/ce675p2.c $(CFLAGS)

2dlaplace.o: problems/2dlaplace.c
	$(CC) -c problems/2dlaplace.c $(CFLAGS)

spheroid.o: problems/spheroid.c
	$(CC) -c problems/spheroid.c $(CFLAGS)

integrate.o: integrate.c integrate.h
	$(CC) -c integrate.c $(CFLAGS)

basis.o: basis.c basis.h
	$(CC) -c basis.c $(CFLAGS)

matrix.o: matrix.c matrix.h
	$(CC) -c matrix.c $(CFLAGS)

mesh.o: mesh.c mesh.h
	$(CC) -c mesh.c $(CFLAGS)

mtxsolver.o: mtxsolver.c mtxsolver.h
	$(CC) -c mtxsolver.c $(CFLAGS)

vector.o: vector.c matrix.h
	$(CC) -c vector.c $(CFLAGS)

finite-element.o: finite-element.c finite-element.h
	$(CC) -c finite-element.c $(CFLAGS)

isoparam.o: isoparam.c isoparam.h
	$(CC) -c isoparam.c $(CFLAGS)

finite-element1d.o: finite-element1d.c finite-element1d.h
	$(CC) -c finite-element1d.c $(CFLAGS)

mesh1d.o: mesh1d.h mesh1d.c
	$(CC) -c mesh1d.c $(CFLAGS)
