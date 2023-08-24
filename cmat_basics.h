#ifndef _CMAT_BASICS_H_
#define _CMAT_BASICS_H_

#include "cmat_defs.h"

enum {LEFT, RIGHT, UP, DOWN, UR, UL, DL, DR};
enum {BY_ROW, BY_COLUMN};

Matrix M_new(int, int);
Matrix M_copy(Matrix);
void set_by_value(Matrix, elem_t);
void set_by_function(Matrix, elem_t (*)(int, int));
void set_by_line(Matrix, elem_t *, char);
void set_row(Matrix, int, elem_t *);
void set_column(Matrix, int, elem_t *);
void clear(Matrix);

Matrix M_get_row(Matrix, int);
Matrix M_get_column(Matrix, int);
elem_t *flatten(Matrix, char);
void print_matrix(Matrix, int);
void print_row(Matrix, int, int);

void swap_row(Matrix, int, int);
void multiply_row(Matrix, int, elem_t);
void combine_row(Matrix, int, int);

Matrix M_add(Matrix, Matrix);
Matrix M_multiply(Matrix, Matrix);
void scalar_multiply(Matrix, elem_t);
Matrix M_transpose(Matrix);
Matrix M_append(const Matrix, const Matrix, const char);

elem_t row_echelon_form(Matrix);
elem_t determinant(Matrix);
Matrix M_inverse(Matrix);
Matrix M_solve(Matrix, Matrix);
Matrix M_basis_transform(Matrix, Matrix);

Matrix M_eigenvectors(Matrix);
elem_t *eigenvalues(Matrix);

#endif