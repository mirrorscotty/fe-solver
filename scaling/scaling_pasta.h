#ifndef SCALING_PASTA_H
#define SCALING_PASTA_H

typedef struct {
    double rhoc; /**< Density/mass concentration */
    double Cpc; /**< Heat capacity */
    double Tc; /**< Temperature */
    double Kc; /**< Thermal conductivity */
    double kwc; /**< Permeability */
    double Lc; /**< Length */
    double Pc; /**< Pressure */
    double muc; /**< Viscosity */
} scaling_pasta;

#endif

