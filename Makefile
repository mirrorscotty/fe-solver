CC=gcc
CFLAGS=-lm -I. -Imatrix -Isolver -Ioutput -Isolver/mesh -Isolver/integration -Isolver/ode -Iscaling -Igui/heating -Wall -g3 -O0 

SRC=$(wildcard solver/*.c) \
    $(wildcard solver/mesh/*.c) \
    $(wildcard solver/ode/*.c) \
    $(wildcard solver/integration/*.c) \
    scaling/scaling_ht.c \
    $(wildcard output/*.c)

all: fe-solver.a

doc:
	doxygen DoxyFile
	make -C doc/latex
	cp doc/latex/refman.pdf doc/Reference.pdf

clean:
	rm -rf *.a
	rm -rf doc
	rm -rf $(SRC:.c=.o)
	rm -rf $(SRC:.c=.d)
	$(MAKE) -C matrix clean

fe-solver.a: $(SRC:.c=.o)
	ar -cvr $@ $?

matrix.a:
	$(MAKE) -C matrix
	cp matrix/matrix.a .

%.o: %.c
	$(CC) -c $(CFLAGS) $*.c -o $*.o
	$(CC) -MM $(CFLAGS) $*.c > $*.d
	@mv -f $*.d $*.d.tmp
	@sed -e 's|.*:|$*.o:|' < $*.d.tmp > $*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
          sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

-include $(OBJ:.o=.d)

