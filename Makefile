CC=gcc

CFLAGS=-lm -Wall -ggdb

all: spheroid 

2dlaplace: integrate.o basis.o matrix.o mesh.o mtxsolver.o vector.o finite-element.o 2dlaplace.o isoparam.o
	gcc -o 2dlaplace 2dlaplace.o integrate.o basis.o matrix.o mesh.o mtxsolver.o vector.o finite-element.o isoparam.o -lm -ggdb

ce675p1: integrate.o basis.o matrix.o mesh.o mtxsolver.o vector.o finite-element.o ce675p1.o
	$(CC) -o ce675p1 ce675p1.o integrate.o basis.o matrix.o mesh.o mtxsolver.o vector.o finite-element.o $(CFLAGS)

spheroid: integrate.o basis.o matrix.o mesh.o mtxsolver.o vector.o finite-element.o isoparam.o spheroid.o
	$(CC) -o spheroid integrate.o basis.o matrix.o mesh.o mtxsolver.o vector.o finite-element.o isoparam.o spheroid.o $(CFLAGS)

clean:
	rm -rf spheroid 2dlaplace ce675p1 *.o *~



ce675p1.o: ce675p1.c

2dlaplace.o: 2dlaplace.c
	$(CC) -c 2dlaplace.c $(CFLAGS)

spheroid.o: spheroid.c
	$(CC) -c spheroid.c $(CFLAGS)

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

