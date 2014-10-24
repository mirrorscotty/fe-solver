/**
 * @file heat-cyl.c
 * Console-only front end to the can sterilization model. The rest of the code
 * is in gui/heating
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "basis.h"
#include "integrate.h"
#include "mesh1d.h"
#include "isoparam.h"
#include "finite-element1d.h"
#include "auxsoln.h"
#include "solve.h"

#include "output.h"
#include "material-data/choi-okos/choi-okos.h"

#include "heat-transfer.h"

#define TREF 500 // K
#define THICKNESS 2.0
#define HCONV 50 

choi_okos *comp_global;

int main(int argc, char *argv[])
{
    Mesh1D *mesh;
    basis *b;
    matrix *IC;
    struct fe1d* problem;

    /* Load a data file if one is supplied. */
    //if(argc == 2)
    //    init(argv[1]);
    comp_global = CreateChoiOkos(0, 0, 0, 1, 0, 0, 0);

    /* Make a linear 1D basis */
    b = MakeLinBasis(1);

    /* Create a uniform mesh */
    mesh = GenerateUniformMesh1D(b, 0.0, THICKNESS, 20);
    
    problem = CreateFE1D(b, mesh,
                         &CreateDTimeMatrix,
                         &CreateElementMatrix,
                         &CreateElementLoad,
                         &ApplyAllBCs,
                         5000);
    problem->nvars = 1; /* Number of simultaneous PDEs to solve */
    problem->dt = 0.001; /* Dimensionless time step size */
    problem->charvals = SetupScaling(alpha(comp_global, TREF), TREF, THICKNESS, k(comp_global, TREF), HCONV);

    /* Set the initial temperature */
    IC = GenerateInitCondConst(problem, 0, scaleTemp(problem->charvals, 273.0));

    /* Apply the initial condition to the problem and set up the transient
     * solver. */
    FE1DTransInit(problem, IC);

    while(problem->t<problem->maxsteps) {
        NLinSolve1DTransImp(problem, NULL);
    }

    PrintScalingValues(problem->charvals);

    printf("Solution at t = %g:\n", uscaleTime(problem->charvals, problem->t*problem->dt));
    PrintSolution(problem, problem->t-1);

    CSVOutFixedNode(problem, 19, "output.csv");

    /* Clean up */
    DestroyFE1D(problem);

    return 0;
}

