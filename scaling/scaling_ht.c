/**
 * @file scaling_ht.c
 * All the stuff for making the heat equation dimensionless. This file can also
 * be used to make the diffusion equation dimensionless. Just let alpha be the
 * diffusion coefficient, h be the convective mass transfer coefficient, and Tc
 * be the characteristic concentration.
 */

#include <stdio.h>

#include "scaling_ht.h"

/**
 * @brief Store all of the values used to make the problem dimensionless
 *
 * @param alpha The characteristic thermal diffusivity
 * @param Tc Characteristic temperature
 * @param Lc Characteristic length
 * @param k Characteristic thermal conductivity
 * @param h Heat transfer coefficient at the surface
 * @returns A data structure containing all of the supplied values
 */
scaling_ht SetupScaling(double alpha, double T0, double Te, double Lc, double k, double h)
{
    scaling_ht stuff;

    stuff.alpha = alpha;
    stuff.Tc = T0; 
    stuff.Te = Te;
    stuff.Lc = Lc;
    stuff.k = k;
    stuff.h = h;

    return stuff;
}

/**
 * @brief Convert temperature to dimensionless temperature
 *
 * The supplied value of T should be in the same units as Tc
 *
 * @param stuff The data structure containing the characteristic values
 * @param T Temperature
 * @returns Dimensionless temperature
 */
double scaleTemp(scaling_ht stuff, double T)
{
    return (T-stuff.Te)/(stuff.Tc-stuff.Te);
    //return T/stuff.Te;
}
/**
 * @brief Convert dimensionless temperature back to normal units
 *
 * @param stuff The data structure containing the characteristic values
 * @param theta Dimensionless temperature
 * @returns Normal temperature
 */
double uscaleTemp(scaling_ht stuff, double theta)
{
    return theta*(stuff.Tc-stuff.Te)+stuff.Te;
    //return theta*stuff.Te;
}

/**
 * @brief Convert length to dimensionless length 
 * @param stuff The data structure containing the characteristic values
 * @param x Length in the same units as Lc
 * @returns Dimensionless length
 */
double scaleLength(scaling_ht stuff, double x)
{
    return x/stuff.Lc;
}

/**
 * @brief Convert dimensionless length back to normal
 * @param stuff The data structure containing the characteristic values
 * @param l Dimensionless length
 * @returns Length in normal units
 */
double uscaleLength(scaling_ht stuff, double l)
{
    return l*stuff.Lc;
}

/**
 * @brief Convert time to dimensionless time
 * @param stuff The data structure containing the characteristic values
 * @param t Time in normal units
 * @returns Dimensionless time
 */
double scaleTime(scaling_ht stuff, double t)
{
    return stuff.alpha*t/(stuff.Lc*stuff.Lc);
}

/** 
 * @brief Convert dimensionless time back to normal
 *
 * @param stuff The data structure containing the characteristic values
 * @param XFo Fourier number (dimensionless time)
 * @returns Time in normal units
 */
double uscaleTime(scaling_ht stuff, double XFo)
{
    return XFo*stuff.Lc*stuff.Lc/stuff.alpha;
}

/**
 * @brief Calculate the Biot number
 *
 * @param stuff The data structure containing the characteristic values
 * @returns The Biot number
 */
double BiotNumber(scaling_ht stuff)
{
    return stuff.h*stuff.Lc/stuff.k;
}

/**
 * @brief Print out everything in the data structure
 * @param stuff The struct of values to print
 */
void PrintScalingValues(scaling_ht stuff)
{
    printf("alpha = %g, T0 = %g, Te = %g, Lc = %g, k = %g, h = %g, Bi = %g\n",
           stuff.alpha, stuff.Tc, stuff.Te, stuff.Lc, stuff.k, stuff.h, BiotNumber(stuff));
}

