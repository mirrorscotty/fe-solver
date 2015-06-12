/**
 * @file output.c
 * Spit out simulation results into a CSV file!
 */

#include <stdio.h>


#include "solver/finite-element.h"
#include "solver/finite-element1d.h"

#include "output/output.h"
#ifdef CALC_ICE_FORMATION
//#include "freezing-gui.h"
#else
//#include "heat-gui.h"
#endif
#include "solver/solution.h"
#include "solver/ode/auxsoln.h"
#include "scaling/scaling_ht.h"

//#include "material-data.h"

//extern choi_okos *comp_global;

/**
 * @brief Output bunches of data for a single node
 *
 * Spits out the time, temperature, etc. for one node n a CSV file. This needs
 * to be changed, depending on whether the freezing model or the sterilization
 * model is built.
 * @param p The problem with the data in it
 * @param row The number of the row (node) to output data for
 * @param filename The name of the file to spit stuff out into
 */
//void CSVOutFixedNode(struct fe1d *p, int row, char *filename)
//{
//    int i;
//    FILE *fp;
//    solution *s;
//    double T;
//
//    fp = fopen(filename, "w+");
//
//    if(!fp) {
//        fprintf(stderr, "Unable to open file for writing.");
//        return;
//    }
//
//    /* Print out the column headers */
//    fprintf(fp, "Time,Temperature,ProdConc,IceMassFrac,Density,HeatCapacity,ThermalConductivity\n");
//
//    /* Print out the values */
//    for(i=0; i<p->maxsteps; i++) {
//        s = FetchSolution(p, i);
//        T = uscaleTemp(p->charvals, val(s->val, row, 0));
//        fprintf(fp, "%g,%g,", uscaleTime(p->charvals, i*p->dt), T);
//        //s = FetchAuxSoln(p, 0, i);
//        //fprintf(fp, "%g,", val(s->val, row, 0));
//       // fprintf(fp, "%g,",  IceMassFrac(T));
//        fprintf(fp, "%g,%g,%g\n", rho(comp_global, T), Cp(comp_global, T), k(comp_global, T));
//    }
//    fprintf(fp, "\n");
//
//    fclose(fp);
//
//    return;
//}

/**
 * @brief Output bunches of data for a single node
 *
 * Spits out the time, temperature, etc. for one node n a CSV file. This needs
 * to be changed, depending on whether the freezing model or the sterilization
 * model is built.
 * @param p The problem with the data in it
 * @param row The number of the row (node) to output data for
 * @param filename The name of the file to spit stuff out into
 */
void CSVOutFixedNode2(struct fe1d *p, int row, char *filename)
{
    int i, n;
    FILE *fp;
    solution *s;
    double T, X;
    Mesh1D *mesh;

    fp = fopen(filename, "w+");

    if(!fp) {
        fprintf(stderr, "Unable to open file for writing.");
        return;
    }

    /* Print out the column headers */
    //fprintf(fp, "Time,Position,Temperature,Density,HeatCapacity,ThermalConductivity\n");
    fprintf(fp, "Time,Position,Temperature\n");

    /* Print out the values */
    for(i=0; i<p->maxsteps; i++) {
        s = FetchSolution(p, i);

        T = val(s->val, row, 0);

        mesh = p->mesh;
        for(n=0; n<p->maxsteps-i-1; n++)
            mesh = mesh->prev;
        X = valV(mesh->nodes, row);

        fprintf(fp, "%g,%g,%g\n", i*p->dt, X, T);
        //fprintf(fp, "%g,%g,%g\n", rho(comp_global, T), Cp(comp_global, T), k(comp_global, T));
    }
    fprintf(fp, "\n");

    fclose(fp);

    return;
}

/**
 * @brief Output bunches of data for a single node
 *
 * Spits out the time, temperature, etc. for one node n a CSV file. This needs
 * to be changed, depending on whether the freezing model or the sterilization
 * model is built.
 * @param p The problem with the data in it
 * @param row The number of the row (node) to output data for
 * @param filename The name of the file to spit stuff out into
 */
void CSVOutFixedNodeDx(struct fe1d *p, int row, char *filename)
{
    int i, n;
    FILE *fp;
    solution *s;
    double dT, X;
    Mesh1D *mesh;

    fp = fopen(filename, "w+");

    if(!fp) {
        fprintf(stderr, "Unable to open file for writing.");
        return;
    }

    /* Print out the column headers */
    //fprintf(fp, "Time,Position,Temperature,Density,HeatCapacity,ThermalConductivity\n");
    fprintf(fp, "Time,Position,Temperature\n");

    /* Print out the values */
    for(i=0; i<p->maxsteps; i++) {
        s = FetchSolution(p, i);


        mesh = p->mesh;
        for(n=0; n<p->maxsteps-i-1; n++)
            mesh = mesh->prev;
        X = valV(mesh->nodes, row);
        dT = EvalDSoln1DG(p, 0, s, valV(p->mesh->orig->nodes, row), 0);

        fprintf(fp, "%g,%g,%g\n", i*p->dt, X, dT);
        //fprintf(fp, "%g,%g,%g\n", rho(comp_global, T), Cp(comp_global, T), k(comp_global, T));
    }
    fprintf(fp, "\n");

    fclose(fp);

    return;
}

/**
 * @brief Spits out data for all the nodes in the simulation at a single time
 * step.
 * @param p Problem structure
 * @param tstep The time step to output data for
 * @param filename Where to save all the data.
 */
//void CSVOutFixedTime(struct fe1d *p, int tstep, char *filename)
//{
//    int i;
//    FILE *fp;
//    solution *T, *c1;
//
//    vector *defmesh;
//#ifdef CALC_ICE_FORMATION
//    defmesh = deformMesh(p, tstep);
//#else
//    defmesh = p->mesh->nodes;
//#endif
//
//    fp = fopen(filename, "w+");
//
//    if(!fp) {
//        fprintf(stderr, "Unable to open file for writing.");
//        return;
//    }
//
//    T = FetchSolution(p, tstep);
//    c1 = FetchAuxSoln(p, 0, tstep);
//
//    fprintf(fp, "Radius,Temperature,ProdConc,IceMassFrac,Density,HeatCapacity,ThermalConductivity\n");
//    for(i=0; i<p->nrows; i++)
//        fprintf(fp, "%g,%g,%g,%g,%g,%g\n",
//                uscaleLength(p->charvals, valV(defmesh, i)),
//                uscaleTemp(p->charvals, val(T->val, i, 0)),
//                val(c1->val, i, 0),
//               // IceMassFrac(uscaleTemp(p->charvals, val(T->val, i, 0))),
//                rho(comp_global, uscaleTemp(p->charvals, val(T->val, i, 0))),
//                Cp(comp_global, uscaleTemp(p->charvals, val(T->val, i, 0))),
//                k(comp_global, uscaleTemp(p->charvals, val(T->val, i, 0))));
//    fprintf(fp, "\n");
//
//    fclose(fp);
//
//    return;
//}
//
