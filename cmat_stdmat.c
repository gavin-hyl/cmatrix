#include <math.h>

#include "cmat_stdmat.h"
#include "cmat_basics.h"

Matrix M_identity(int d)
{
    Matrix iden = M_new(d, d);
    for (int i = 0; i < d; i++)
    {
        iden->elements[i][i] = 1;
    }
    return iden;
}

Matrix M_combine_row(int d, int fr, int to)
{
    Matrix combrow = M_identity(d);
    combrow->elements[to][fr] = 1;
    return combrow;
}

Matrix M_swap_row(int d, int r1, int r2)
{
    Matrix swaprow = M_identity(d);
    swaprow->elements[r1][r1] = 0;
    swaprow->elements[r2][r2] = 0;
    swaprow->elements[r1][r2] = 1;
    swaprow->elements[r2][r1] = 1;
    return swaprow;
}

Matrix M_multiply_row(int d, int r, elem_t mult)
{
    Matrix multrow = M_identity(d);
    multrow->elements[r][r] = mult;
    return multrow;
}

Matrix M_2d_rotation(elem_t angle)
{
    Matrix mat = M_new(2, 2);
    
    mat->elements[0][0] = cos(angle);
    mat->elements[1][0] = sin(angle);
    mat->elements[0][1] = -sin(angle);
    mat->elements[1][1] = cos(angle);
    return mat;
}