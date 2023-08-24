#ifndef _CMAT_STDMAT_H_
#define _CMAT_STDMAT_H_

#include "cmat_defs.h"

Matrix M_identity(int);
Matrix M_swap_row(int, int, int);
Matrix M_combine_row(int, int, int);
Matrix M_multiply_row(int, int, elem_t);
Matrix M_2d_rotation(elem_t);
Matrix M_2d_flip(Matrix);
Matrix M_3d_rotation(elem_t, elem_t);
Matrix M_3d_axis_rotation(elem_t, Matrix);
Matrix M_3d_flip(Matrix, Matrix);

#endif