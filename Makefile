all: 2dlaplace ce675p1

2dlaplace:
	gcc -o 2dlaplace 2dlaplace.c integrate.c basis.c matrix.c mesh.c mtxsolver.c vector.c finite-element.c -lm -ggdb

ce675p1:
	gcc -o ce675p1 ce675p1.c integrate.c basis.c matrix.c mesh.c mtxsolver.c vector.c -lm -ggdb

clean:
	rm 2dlaplace ce675p1