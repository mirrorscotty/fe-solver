CC=gcc

CFLAGS=-lm -Wall

all: 2dlaplace

2dlaplace: integrate.o basis.o matrix.o mesh.o mtxsolver.o vector.o finite-element.o 2dlaplace.o
	gcc -o 2dlaplace 2dlaplace.c integrate.c basis.c matrix.c mesh.c mtxsolver.c vector.c finite-element.c -lm -ggdb

ce675p1: integrate.o basis.o matrix.o mesh.o mtxsolver.o vector.o finite-element.o ce675p1.o
	$(CC) -o ce675p1 ce675p1.o integrate.o basis.o matrix.o mesh.o mtxsolver.o vector.o finite-element.o $(CFLAGS)
#	gcc -o ce675p1 ce675p1.c integrate.c basis.c matrix.c mesh.c mtxsolver.c vector.c finite-element.c -lm -ggdb

clean:
	rm -rf 2dlaplace ce675p1 *.o *~

ce675p1.o: ce675p1.c

2dlaplace.o: 2dlaplace.c
	$(CC) -c 2dlaplace.c $(CFLAGS)

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
