#ifndef _CMATRIX_ALGEBRA_H_
#define _CMATRIX_ALGEBRA_H_

#include "cmatrix_defs.h"

void swap_row(Matrix, int, int);
void multiply_row(Matrix, int, elem_t);
void combine_row(Matrix, int, int);

elem_t row_echelon_form(Matrix);
elem_t determinant(Matrix);
Vector solve_linear_system(Matrix, Vector);
Matrix basis_transform(Matrix, Matrix);
int linear_independence(Vector *, const int);
Matrix *QR_decomposition(const Matrix);
Matrix QR_decomposition_R_only(const Matrix);


#endif