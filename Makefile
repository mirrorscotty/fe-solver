all:
	gcc -o fe -ggdb integrate.c matrix.c basis.c fe-solver.c mtxsolver.c -lm

