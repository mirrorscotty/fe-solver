#include "arclength.h"
#include <math.h>

double ALength(struct fe *p, matrix *guess, matrix *prevguess, Elem2D *elem,
               double x, double y, double l, double lprev)
{
    basis *b;
    b = p->b;
    int i;
    double phi;
    double g, pg;
    double result = 0;
    
    for(i=0; i<b->n*p->nvars; i++) {
        phi = b->Eval2D(p->b, i/p->nvars, x, y);
        
        g = val(guess, i, 0);
        pg = val(prevguess, i, 0);
        result += pow((g-pg)*phi, 2);
    }
    
    return result;
}

double ALengthx(struct fe *p, matrix *guess, matrix *prevguess, Elem2D *elem,
                double x, double y, double l, double lprev)
{
    basis *b;
    b = p->b;
    int i;
    double phi;
    double g, pg;
    double result = 0;
    
    for(i=0; i<b->n*p->nvars; i++) {
        phi = IEvalLin2Dx(p->b, i/p->nvars, x, y);
        
        g = val(guess, i, 0);
        pg = val(prevguess, i, 0);
        result += pow((g-pg)*phi, 2);
    }
    
    return result;
}

double ALengthy(struct fe *p, matrix *guess, matrix *prevguess, Elem2D *elem,
                double x, double y, double l, double lprev)
{
    basis *b;
    b = p->b;
    int i;
    double phi;
    double g, pg;
    double result = 0;
    
    for(i=0; i<b->n*p->nvars; i++) {
        phi = IEvalLin2Dy(p->b, i/p->nvars, x, y);
        
        g = val(guess, i, 0);
        pg = val(prevguess, i, 0);
        result += pow((g-pg)*phi, 2);
    }
    
    return result;
}

matrix *AssembleBottomRow(struct fe *p, matrix *guess, matrix *prevguess)
{
    matrix *botrow;
    matrix *tmp1, *tmp2;
    int cols = p->nvars*p->nrows;

    botrow = CreateMatrix(1, cols);

    int i, j;
    int col;
    double l, lprev, value;
    for(i=0; i<p->mesh->nelemx*p->mesh->nelemy; i++) {
        tmp1 = GetLocalGuess(p, guess, i);
        tmp2 = GetLocalGuess(p, prevguess, i);

        for(j=0; j<p->b->n; j++) {
            /* Determine the global node number */
            col = valV(p->mesh->elem[i]->map, j);

            /* Set the value for the x derivative */
            value = quad2d3arc(p, guess, prevguess, p->mesh->elem[i], &ALengthx, l, lprev);
            addval(botrow, value, 0, 2*col);

            /* Do the same for the y derivative */
            value = quad2d3arc(p, guess, prevguess, p->mesh->elem[i], &ALengthy, l, lprev);
            addval(botrow, value, 0, 2*col+1);
        }

        DestroyMatrix(tmp1);
        DestroyMatrix(tmp2);
    }
    return botrow;
}
        
matrix *AssembleRightCol(struct fe *p, matrix *guess, matrix *prevguess)
{
    matrix *rightcol;
    matrix *tmp1, *tmp2;
    int rows = p->nvars*p->nrows;

    rightcol = CreateMatrix(rows, 1);

    int i, j;
    int row;
    double l, lprev, value;

    for(i=0; i<p->mesh->nelemx*p->mesh->nelemy; i++) {
        tmp1 = GetLocalGuess(p, guess, i);
        tmp2 = GetLocalGuess(p, prevguess, i);

        for(j=0; j<p->b->n; j++) {
            row = valV(p->mesh->elem[i]->map, j);

            value = quad2d3arc(p, guess, prevguess, p->mesh->elem[i], &ALengthx, l, lprev);
            addval(rightcol, value, 2*row, 0);
            addval(rightcol, value, 2*row+1, 0);

        }
        DestroyMatrix(tmp1);
        DestroyMatrix(tmp2);
    }
    return rightcol;
}

double AssembleBotCorner(struct fe *p, matrix *guess, matrix *prevguess)
{
    return 0;
}

/* Todo: fix this! */
double AssembleResidual(struct fe *p, matrix *guess, matrix *prevguess)
{
    return 0;
}

void AddArcLengthConstraint(struct fe *p, matrix *guess, matrix *prevguess)
{
    matrix *Jnew, Fnew;
    matrix *rightcol, *botrow;
    int i, j;
    int nrows = mtxlen2(p->J);
    double corner, res;

    p->nconstr = 1;

    Jnew = CreateMatrix(nrows+p->nconstr, nrows+p->nconstr);
    Fnew = CreateMatrix(nrows+p->nconstr, 1);

    botrow = AssembleBottomRow(p, guess, prevguess);
    rightcol = AssembleRightCol(p, guess, prevguess);
    corner = AssembleBotCorner(p, guess, prevguess);
    res = AssembleResidual(p, guess, prevguess);
    
    /* Copy the values of the original Jacobian to the new matrix. */
    for(i=0; i<nrows; i++) {
        for(j=0; j<nrows; j++) {
            setval(Jnew, val(p->J, i, j), i, j);
        }
        setval(Fnew, val(p->F, i, 0), i, 0);
    }

    /* Copy over the bottom row and right column. */
    for(i=0; i<nrows; i++) {
        setval(Jnew, val(botrow, 0, i), nrows, i);
        setval(Jnew, val(rightcol, i, 0), i, nrows);
    }

    /* Set the value of the bottom right corner. */
    setval(Jnew, corner, nrows, nrows);

    /* Set the value in the new residual matrix. */
    setval(Fnew, res, nrows, 0);

    /* Free the unneeded matricies */
    DestroyMatrix(rightcol);
    DestroyMatrix(botrow);

    p->Jconstr = Jnew;
    p->Fconstr = Fnew;
    
    return;
}

