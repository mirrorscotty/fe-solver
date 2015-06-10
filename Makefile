CC=gcc
CFLAGS=-lm -I. -Imatrix -Imaterial-data -Isolver -Ioutput -Isolver/mesh -Isolver/integration -Isolver/ode -Iscaling -Igui/heating -Wall -g3 -O0 
VPATH=solver/mesh solver/ode solver/integration matrix scaling solver output

OBJECTFILES=integrate.o basis.o mesh2d.o finite-element.o isoparam.o finite-element1d.o mesh1d.o solution.o auxsoln.o scaling_ht.o linsolve1d.o nlinsolve1d.o predict1d.o linsolve2d.o nlinsolve2d.o output.o stepsize.o

all: fe-solver.a

doc:
	doxygen DoxyFile
	make -C doc/latex
	cp doc/latex/refman.pdf doc/Reference.pdf

clean:
	rm -rf *.o *.a
	$(MAKE) -C matrix clean
	$(MAKE) -C material-data clean
	rm -rf doc

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
stepsize.o: solve.h

# Scaling stuff
scaling_ht.o: scaling_ht.h
scaling_pasta.o: scaling_pasta.h

# Random files that don't have a home
about.o: about.h
basis.o: basis.h
isoparam.o: isoparam.h
integrate.o: integrate.c integrate.h
output.o: output.h

fe-solver.a: $(OBJECTFILES)
	ar -cvr $@ $?

matrix.a:
	$(MAKE) -C matrix
	cp matrix/matrix.a .

