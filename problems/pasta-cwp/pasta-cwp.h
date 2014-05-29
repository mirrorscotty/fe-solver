#ifndef PASTA_CWP_H
#define PASTA_CWP_H

#define D_LIQ(CW, T) \
    scaleDiff((S), DiffCh10(uscaleConc((S), (CW)), uscaleTemp((S), (T))))
#define D_GAS(CW, T) \
    scaleDiff((S), VaporDiffCh10(uscaleConc((S), (CW)), uscaleTemp((S), (T))))
#define K_WAT(CW, PHI, T) \
    scalePerm((S), perm_wat(uscaleConc((S), (CW)), (PHI), uscaleTemp((S), (T))))
#define K_GAS(CW, PHI, T) \
    scalePerm((S), perm_gas(uscaleConc((S), (CW)), (PHI), uscaleTemp((S), (T))))
#define S_GAS(CW, PHI, T) \
    sat_gas(uscaleConc((S), (CW)), (PHI), uscaleTemp((S), (T)))
#define X_VAP(WV) molefrac_vap((wv))
#define RHO_GAS(WV, T, P) \
    scaleConc((S), rho_gas((WV), uscalePres((S), (P)), uscaleTemp((S), (T))))
#define CG_STAR(P) 0

#define C_VAP(CW, WV, PHI, T, P)\
    scaleConc((S) conc_vap(uscaleConc((S), (CW)),\
                           (WV),\
                           (PHI),\
                           uscaleTemp((S), (T)),\
                           uscalePress((S), (T))))
                
#define C_GAS(CW, WV, PHI, T, P)\
    scaleConc((S) conc_gas(uscaleConc((S), (CW)),\
                           (WV),\
                           (PHI),\
                           uscaleTemp((S), (T)),\
                           uscalePress((S), (T))))

#define VARCW 0
#define VARWV 1
#define VARP 2

#endif

