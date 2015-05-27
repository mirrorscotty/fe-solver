#ifndef STEPSIZE_H
#define STEPSIZE_H

#define ERR 1e-3 /* Desired error */

double CorrectorError(struct fe1d *, matrix *, matrix *);
double StepSizeBase(struct fe1d *, matrix *, matrix *);
double StepSize(struct fe1d *, matrix *, matrix *);

#endif

