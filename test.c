/*
 * Series of test functions for figuring stuff out. Pretty much nothing here is
 * used for actual analysis.
 */

#include <stdio.h>
#include "basis.h"

int main(int argc, char *argv[])
{
    basis *b;
    b = MakeLinBasis(2);

    PrintNodeLocations(b);

    DestroyBasis(b);
    return 0;
}

/* Verify that the node locations are where they should be. Just prints out the
 * values of the four basis functions for a bilinear element at each of the four
 * nodes. All the values should be 1.
 */
void PrintNodeLocations(basis *b)
{
    printf("Phi0(%d, %d) = %g\n", 0, 0, EvalLin2D(b, 0, 0, 0));
    printf("Phi1(%d, %d) = %g\n", 0, 1, EvalLin2D(b, 1, 0, 1));
    printf("Phi2(%d, %d) = %g\n", 1, 0, EvalLin2D(b, 2, 1, 0));
    printf("Phi3(%d, %d) = %g\n", 1, 1, EvalLin2D(b, 3, 1, 1));
}
