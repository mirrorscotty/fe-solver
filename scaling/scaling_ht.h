#ifndef SCALING_HT_H
#define SCALING_HT_H

/**
 * @struct scaling_ht
 * @brief A data structure to hold all of the values needed to make a heat
 * transfer problem dimensionless
 * @var scaling_ht::alpha
 * Thermal Diffusivity
 * @var scaling_ht::Tc
 * Characteristic Temperature
 * @var scalin_ht::Lc
 * Characteristic length
 * @var scaling_ht::k
 * Thermal Conductivity
 * @var scaling_ht::h
 * Convective heat transfer coefficient
 */
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

void PrintScalingValues(scaling_ht);

#endif
