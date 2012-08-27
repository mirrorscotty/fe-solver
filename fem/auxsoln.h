#ifndef AUXSOLN_H
#define AUXSOLN_H

struct fe1d;

void FE1DInitAuxSolns(struct fe1d*, int);
void SolveODE(struct fe1d*, int, int, double (*)(double, double, double), double);

void PrintAuxSoln(struct fe1d*, int, int);

#endif
