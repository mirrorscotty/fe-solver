scaling_pasta SetupScaling()
{
    scaling_pasta s;

    choi_okos = CreateChoiOkos(WATERCOMP);
    s.Tc = 60+273.15; /* 60 C */

//    s.rhoc = rho(co, s.Tc);
//    s.Cpc = Cp(co, s.Tc);
//    s.Kc = k(co, s.Tc);

    s.muc = visc_wat(T);

    s.Lc = 0;//???

    s.kwc = 0;//???

    s.Pc = 101300; /* 101.3 kPa */

    return s;
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
double scaleTemp(scaling_pasta s, double T)
{
    return T/f.Tc;
}

/**
 * @brief Convert dimensionless temperature back to normal units
 *
 * @param stuff The data structure containing the characteristic values
 * @param theta Dimensionless temperature
 * @returns Normal temperature
 */
double uscaleTemp(scaling_pasta s, double theta)
{
    return theta*s.Tc;
}

/**
 * @brief Convert length to dimensionless length
 * @param stuff The data structure containing the characteristic values
 * @param x Length in the same units as Lc
 * @returns Dimensionless length
 */
double scaleLength(scaling_ht s, double x)
{
    return x/s.Lc;
}

/**
 * @brief Convert dimensionless length back to normal
 * @param stuff The data structure containing the characteristic values
 * @param l Dimensionless length
 * @returns Length in normal units
 */
double uscaleLength(scaling_pasta s, double l)
{
    return l*s.Lc;
}

/**
 * @brief Convert time to dimensionless time
 * @param stuff The data structure containing the characteristic values
 * @param t Time in normal units
 * @returns Dimensionless time
 */
double scaleTime(scaling_pasta s, double t)
{
    double alpha = s.Kc/(s.rhoc*s.Cpc);
    return alpha*t/(s.Lc*s.Lc);
}

/**
 * @brief Convert dimensionless time back to normal
 *
 * @param stuff The data structure containing the characteristic values
 * @param XFo Fourier number (dimensionless time)
 * @returns Time in normal units
 */
double uscaleTime(scaling_pasta s, double XFo)
{
    double alpha = s.Kc/(s.rhoc*s.Cpc);
    return XFo*s.Lc*s.Lc/alpha;
}

double scaleConc(scaling_pasta s, double c)
{
    return c/s.rhoc;
}

double uscaleConc(scaling_pasta s, double xi)
{
    return xi*s.rhoc;
}

double scaleDiff(scaling_pasta s, double D)
{
}
double uscaleDiff(scaling_pasta s, double D) {}

double scalePerm(scaling_pasta s, double k) {}
double uscalePerm(scaling_pasta s, double K) {}

/**
 * @brief Print out everything in the data structure
 * @param stuff The struct of values to print
 */
void PrintScalingValues(scaling_pasta s)
{
    printf("Rho = %g, Cp = %g, T = %g, K = %g, L = %g, P = %g, mu = %g\n",
            s.rhoc, s.Cpc, s.Tc, s.Kc, s.Lc, s.Pc, s.muc);
}

double DarcyNumber(scaling_pasta s)
{
    return s.kwc/(s.Lc*s.Lc);
}

