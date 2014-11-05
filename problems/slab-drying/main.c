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
#include "material-data.h"

#include "common.h"

choi_okos *comp_global;

int main(int argc, char *argv[])
{
    Mesh1D *mesh;
    basis *b;
    matrix *IC, *IC_heat, *IC_mass;
    struct fe1d* problem;
    scaling_ht scale_heat;
    scaling_ht scale_mass;

    comp_global = CreateChoiOkos(0, 0, 0, 1, 0, 0, 0);
    //scale_heat = SetupScaling(alpha(comp_global, TINIT), TINIT, TAMB, THICKNESS, k(comp_global, TINIT), HCONV);
    scale_mass = SetupScaling(DIFF(CINIT, TINIT), CINIT, CAMB, THICKNESS, DIFF(CINIT, TINIT), KC_CONV);

    /* Make a linear 1D basis */
    b = MakeLinBasis(1);

    /* Create a uniform mesh */
    //mesh = GenerateUniformMesh1D(b, 0.0, scaleLength(scale_heat, THICKNESS), 10);
    mesh = GenerateUniformMesh1D(b, 0.0, scaleLength(scale_mass, THICKNESS), 10);
    
    problem = CreateFE1D(b, mesh,
                         &CreateDTimeMatrix,
                         &CreateElementMatrix,
                         &CreateElementLoad,
                         &ApplyAllBCs,
                         1000);
    problem->nvars = 1; /* Number of simultaneous PDEs to solve */
    problem->dt = 0.001; /* Dimensionless time step size */
    //problem->charvals = scale_heat;
    problem->chardiff = scale_mass;

    /* Set the initial temperature */
    //IC_heat = GenerateInitCondConst(problem, TVAR, scaleTemp(problem->charvals, TINIT));
    IC_mass = GenerateInitCondConst(problem, CVAR, scaleTemp(problem->chardiff, CINIT));
    //IC = mtxadd(IC_heat, IC_mass);
    //DestroyMatrix(IC_heat, IC_mass);

    /* Apply the initial condition to the problem and set up the transient
     * solver. */
    //FE1DTransInit(problem, IC_heat);
    FE1DTransInit(problem, IC_mass);

    while(problem->t<problem->maxsteps) {
        NLinSolve1DTransImp(problem, NULL);
        if(problem->t-1 > 0)
            problem->mesh = 
                MoveMeshF(problem, problem->mesh->orig,
                          problem->t-1, &DeformationGrad);
    }

    //PrintScalingValues(problem->charvals);
    PrintScalingValues(problem->chardiff);

    //printf("Solution at t = %g:\n", uscaleTime(problem->charvals, problem->t*problem->dt));
    printf("Solution at t = %g:\n", uscaleTime(problem->chardiff, problem->t*problem->dt));
    PrintSolution(problem, problem->t-1);

    CSVOutFixedNode2(problem, 0, "output00.csv");
    CSVOutFixedNode2(problem, 1, "output01.csv");
    CSVOutFixedNode2(problem, 2, "output02.csv");
    CSVOutFixedNode2(problem, 3, "output03.csv");
    CSVOutFixedNode2(problem, 4, "output04.csv");
    CSVOutFixedNode2(problem, 5, "output05.csv");
    CSVOutFixedNode2(problem, 6, "output06.csv");
    CSVOutFixedNode2(problem, 7, "output07.csv");
    CSVOutFixedNode2(problem, 8, "output08.csv");
    CSVOutFixedNode2(problem, 9, "output09.csv");
    CSVOutFixedNode2(problem, 10, "output10.csv");

    PrintVector(problem->mesh->orig->nodes);
    PrintVector(problem->mesh->nodes);

    /* Clean up */
    DestroyFE1D(problem);

    return 0;
}

