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

#include "material-data/choi-okos/choi-okos.h"

#include "heat-gui.h"

double EaA = 10000,
       EaB = 20000,
       AA = 0.9,
       AB = 0.8,
       Text_hot = 500,
       Text_cold = 290,
       t_heat = 500;

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
    mesh = GenerateUniformMesh1D(b, 0.0, 2.0, 6);
    
    problem = CreateFE1D(b, mesh,
                         &CreateDTimeMatrix,
                         &CreateElementMatrix,
                         &CreateElementLoad,
                         &ApplyAllBCs,
                         600);
    problem->nvars = 1;
    problem->dt = .001;

    IC = GenerateInitCondConst(problem, 0, 273.0); /* Initial temperature */
    //ApplyInitialCondition(problem, 1, &InitC); /* Initial concentration */
    //
    FE1DTransInit(problem, IC);

    while(problem->t<problem->maxsteps) {
        LinSolve1DTransImp(problem);
    }
    printf("t = %g\n", (problem->maxsteps-1) * problem->dt);

    puts("Solutions:");
    //PrintSolution(problem, 0);
    //puts("");
    //PrintSolution(problem, 1);
    //puts("");
    PrintSolution(problem, problem->t-1);

    FE1DInitAuxSolns(problem, 2);
    SolveODE(problem, 0, 0, &react1, 1);
    SolveODE(problem, 0, 1, &react2, 2);
    puts("");
    PrintAuxSoln(problem, 0, 1);
    puts("");
    PrintAuxSoln(problem, 1, 1);

    /* Clean up */
    DestroyFE1D(problem);

    return 0;
}

