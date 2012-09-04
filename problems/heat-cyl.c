#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix.h"
#include "basis.h"
#include "mtxsolver.h"
#include "integrate.h"
#include "mesh1d.h"
#include "isoparam.h"
#include "finite-element1d.h"
#include "auxsoln.h"

#include "material-data/freezing/freezing.h"

#include "heat-gui.h"

int main(int argc, char *argv[])
{
    Mesh1D *mesh;
    basis *b;
    matrix *IC;
    struct fe1d* problem;

    /* Load a data file if one is supplied. */
    //if(argc == 2)
    //    init(argv[1]);

    /* Make a linear 1D basis */
    b = MakeLinBasis(1);

    /* Create a uniform mesh */
    mesh = GenerateUniformMesh1D(b, 0.0, 2.0, 6);
    
    problem = CreateFE1D(b, mesh,
                         &CreateDTimeMatrix,
                         &CreateElementMatrix,
                         &CreateElementLoad,
                         &ApplyAllBCs,
                         5);
    problem->nvars = 1;
    problem->dt = .001;

    IC = GenerateInitCondConst(problem, 0, 273.0); /* Initial temperature */
    //ApplyInitialCondition(problem, 1, &InitC); /* Initial concentration */
    //
    FE1DTransInit(problem, IC);

    while(problem->t<problem->maxsteps) {
        LinSolve1DTransImp(problem);
    }
    //printf("t = %g\n", (problem->maxsteps-1) * problem->dt);

    puts("Solutions:");
    PrintSolution(problem, 0);
    puts("");
    PrintSolution(problem, 1);
    puts("");
    PrintSolution(problem, 2);

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

