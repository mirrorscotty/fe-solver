/**
 * @file scaling_ht.h
 * Characteristic values for heat transfer
 */

#ifndef SCALING_HT_H
#define SCALING_HT_H

/**
 * @struct scaling_ht
 * @brief A data structure to hold all of the values needed to make a heat
 * transfer problem dimensionless
 */
typedef struct {
    double alpha; /**< Thermal Diffusivity */
    double Tc; /**< Initial Temperature */
    double Te; /**< Equilibrium Temperature */
    double Lc; /**< Characteristic Length */
    double k; /**< Thermal Conductivity */
    double h; /**< Convective heat transfer coefficient */
} scaling_ht;

scaling_ht SetupScaling(double, double, double, double, double, double);

double scaleTemp(scaling_ht, double);
double uscaleTemp(scaling_ht, double);

double scaleLength(scaling_ht, double);
double uscaleLength(scaling_ht, double);

double scaleTime(scaling_ht, double);
double uscaleTime(scaling_ht, double);

double BiotNumber(scaling_ht);

void PrintScalingValues(scaling_ht);

#endif
