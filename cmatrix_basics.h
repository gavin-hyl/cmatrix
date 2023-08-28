#ifndef _CMATRIX_BASICS_H_
#define _CMATRIX_BASICS_H_

#include "cmatrix_defs.h"

Matrix new_matrix(int, int);
Matrix copy_matrix(Matrix);
int compare_matrix(Matrix, Matrix);
void set_matrix_by_value(Matrix, elem_t);
void set_matrix_by_function(Matrix, elem_t (*)(int, int));
void set_matrix_by_line(Matrix, elem_t *, char);
void set_matrix_row(Matrix, int, elem_t *);
void set_matrix_column(Matrix, int, elem_t *);
Matrix get_row_matrix(Matrix, int);
Matrix get_column_matrix(Matrix, int);
void clear_matrix(Matrix);
elem_t *flatten(Matrix, char);
void print_matrix(Matrix);
void print_row(Matrix, int);
Matrix add_matrix(Matrix, Matrix);
Matrix multiply_matrix(Matrix, Matrix);
void scale_matrix(Matrix, elem_t);
Matrix matrix_power(Matrix, int);
Matrix transpose(Matrix);
Matrix inverse(Matrix);
Matrix append_horizontal(const Matrix, const Matrix);

Vector new_vector(int);
Vector copy_vector(Vector);
int compare_vector(Vector, Vector);
void set_vector_by_value(Vector, elem_t);
void set_vector_by_function(Vector, elem_t (*)(int));
void set_vector(Vector, elem_t *);
void clear_vector(Vector);
void print_vector(Vector);
Vector add_vector(Vector, Vector);
Vector multiply_matrix_vector(Matrix, Vector);
void scale_vector(Vector, elem_t);
elem_t vector_length(Vector);
void normalize(Vector);
elem_t dot_product(Vector, Vector);
Vector cross_product(Vector, Vector);

Matrix vector_to_row_matrix(Vector);
Matrix vector_to_column_matrix(Vector);
Vector matrix_to_vector(Matrix);
Vector get_row_vector(Matrix, int);
Vector get_column_vector(Matrix, int);

Vector std_unit_vector(int, int);
Matrix identity(int);

#endif