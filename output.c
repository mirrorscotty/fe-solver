#include <stdio.h>

#include "output.h"

#include "finite-element.h"
#include "freezing-gui.h"
#include "solution.h"
#include "auxsoln.h"
#include "scaling_ht.h"

#include "material-data/choi-okos/choi-okos.h"

void CSVOutFixedNode(struct fe1d *p, int row, char *filename)
{
    int i;
    FILE *fp;
    solution *s;
    double T;

    fp = fopen(filename, "w+");

    if(!fp) {
        fprintf(stderr, "Unable to open file for writing.");
        return;
    }

    /* Print out the column headers */
    fprintf(fp, "Time,Temperature,ProdConc,IceMassFrac,Density,HeatCapacity,ThermalConductivity\n");

    /* Print out the values */
    for(i=0; i<p->maxsteps; i++) {
        s = FetchSolution(p, i);
        T = uscaleTemp(p->charvals, val(s->val, row, 0));
        fprintf(fp, "%g,%g,", uscaleTime(p->charvals, i*p->dt), T);
        s = FetchAuxSoln(p, 0, i);
        fprintf(fp, "%g,", val(s->val, row, 0));
        fprintf(fp, "%g,",  IceMassFrac(T));
        fprintf(fp, "%g,%g,%g\n", rho(T), Cp(T), k(T));
    }
    fprintf(fp, "\n");

    fclose(fp);

    return;
}

void CSVOutFixedTime(struct fe1d *p, int tstep, char *filename)
{
    int i;
    FILE *fp;
    solution *T, *c1;

    vector *defmesh;
#ifdef CALC_ICE_FORMATION
    defmesh = deformMesh(p, tstep);
#else
    defmesh = p->mesh->nodes;
#endif

    fp = fopen(filename, "w+");

    if(!fp) {
        fprintf(stderr, "Unable to open file for writing.");
        return;
    }

    T = FetchSolution(p, tstep);
    c1 = FetchAuxSoln(p, 0, tstep);

    fprintf(fp, "Radius,Temperature,ProdConc,IceMassFrac,Density,HeatCapacity,ThermalConductivity\n");
    for(i=0; i<p->nrows; i++)
        fprintf(fp, "%g,%g,%g,%g,%g,%g,%g\n",
                uscaleLength(p->charvals, valV(defmesh, i)),
                uscaleTemp(p->charvals, val(T->val, i, 0)),
                val(c1->val, i, 0),
                IceMassFrac(uscaleTemp(p->charvals, val(T->val, i, 0))),
                rho(uscaleTemp(p->charvals, val(T->val, i, 0))),
                Cp(uscaleTemp(p->charvals, val(T->val, i, 0))),
                k(uscaleTemp(p->charvals, val(T->val, i, 0))));
    fprintf(fp, "\n");

    fclose(fp);

    return;
}

