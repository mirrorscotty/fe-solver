#ifndef STEPSIZE_H
#define STEPSIZE_H

#define ERR 1e-4 /* Desired error */

double CorrectorError(struct fe1d *, matrix *, matrix *);
double StepSizeBase(struct fe1d *, matrix *, matrix *);
double StepSize(struct fe1d *, matrix *, matrix *);
double CurrentTime(struct fe1d *, int);

#endif

