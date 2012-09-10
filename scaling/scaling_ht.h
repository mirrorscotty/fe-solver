#ifndef SCALING_HT_H
#define SCALING_HT_H

typedef struct {
    double alpha; /* Thermal Diffusivity */
    double Tc; /* Characteristic Temperature */
    double Lc; /* Characteristic Length */
    double k; /* Thermal Conductivity */
    double h; /* Convective heat transfer coefficient */
} scaling_ht;

scaling_ht SetupScaling(double, double, double, double, double);

double scaleTemp(scaling_ht, double);
double uscaleTemp(scaling_ht, double);

double scaleLength(scaling_ht, double);
double uscaleLength(scaling_ht, double);

double scaleTime(scaling_ht, double);
double uscaleTime(scaling_ht, double);

double BiotNumber(scaling_ht);

#endif
