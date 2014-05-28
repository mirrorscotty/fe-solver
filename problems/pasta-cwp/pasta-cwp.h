#ifndef PASTA_CWP_H
#define PASTA_CWP_H

#define Dliq(CW, T) DiffCh10((CW), (T))
#define Kw(CW, PHI, T) perm_wat((CW), (PHI), (T))
#define Kg(CW, PHI, T) perm_gas((CW), (PHI), (T))

#define VARCW 0
#define VARWV 1
#define VARP 2


#endif

