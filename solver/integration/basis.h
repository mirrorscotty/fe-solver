#ifndef BASIS_H
#define BASIS_H

/**
 * @struct _basis
 * @brief A data structure to hold interpolation functions
 * @var _basis::n
 * Number of functions in the basis
 * @var _basis::overlap
 * Amount of overlap when assembling the global matricies
 * @var _basis::phi
 * Array of pointers to the basis functions
 * @var _basis::dphi
 * Array of pointers to the first derivatives of the basis functions
 * @var _basis::phi2d
 * Only used for 2D triangular meshes
 * @var _basis::dim
 * The dimension of the problem
 * @var _basis::Eval2D
 * Pointer to the function used to evaluate the basis functions in 2D
 * @var _basis::Eval2Dx
 * Same thing as Eval2D only for the derivative with respect to x
 * @var _basis::Eval2Dy
 * Same as Eval2D, only for the y derivative
 */
struct _basis {
    int n; /* Number of functions in the basis */
    int overlap; /* Amount of overlap when assembling the global matricies */
    double (**phi)(double); /* Array of basis functions */
    double (**dphi)(double); /* First derivatives of the basis functions */
    double (**phi2d)(double, double); /* Only for 2d triangles */
    int dim; /* Dimension of the problem */

    double (*Eval2D)(struct _basis*, int, double, double);
    double (*Eval2Dx)(struct _basis*, int, double, double);
    double (*Eval2Dy)(struct _basis*, int, double, double);
};
typedef struct _basis basis;

basis* MakeLinBasis(int);
basis* MakeQuadBasis(int);
basis* MakeCubicBasis(int);

void DestroyBasis(basis*);

double EvalLin2D(basis*, int, double, double);
double EvalLin2Dx(basis*, int, double, double);
double EvalLin2Dy(basis*, int, double, double);
//double EvalBasis(basis*, ... );

double EvalQuad2D(basis*, int, double, double);
double EvalQuad2Dx(basis*, int, double, double);
double EvalQuad2Dy(basis*, int, double, double);

double EvalLinTri2D(basis*, int, double, double);
double EvalLinTri2Dx(basis*, int, double, double);
double EvalLinTri2Dy(basis*, int, double, double);

#endif
