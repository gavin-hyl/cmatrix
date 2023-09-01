/**
 * @file cmatrix_geometry.h
 * @author Gavin Hua (139950129+GavinHYL@users.noreply.github.com)
 * @brief More advanced algebratic operations not included in cmatrix_basics.h
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef _CMATRIX_ALGEBRA_H_
#define _CMATRIX_ALGEBRA_H_

#include "cmatrix_defs.h"

void swap_row(Matrix, int, int);
void multiply_row(Matrix, int, flt_t);
void combine_row(Matrix, int, int);
Matrix row_swap_matrix(int d, int r1, int r2);
Matrix row_multiply_matrix(int d, int r, flt_t mult);
Matrix row_combine_matrix(int d, int fr, int to);

flt_t row_echelon_form(Matrix);
flt_t determinant(Matrix);
Vector solve_linear_system(Matrix, Vector);
Matrix basis_transform(Matrix, Matrix);
int linear_independence(Vector *, const int);
Matrix *QR_decomposition(Matrix);
Matrix QR_decomposition_R_only(Matrix);
Vector eigenvalues(Matrix);

#endif