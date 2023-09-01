/**
 * @file cmatrix_basics.h
 * @author Gavin Hua (139950129+GavinHYL@users.noreply.github.com)
 * @brief Basic unctions for creating and manipulating matrices and vectors.
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef _CMATRIX_BASICS_H_
#define _CMATRIX_BASICS_H_

#include "cmatrix_defs.h"

Matrix new_matrix(int, int);
void free_matrix(Matrix);
Matrix copy_matrix(Matrix);
int matrix_is_equal(Matrix, Matrix);
void set_matrix_by_value(Matrix, flt_t);
void set_matrix_by_function(Matrix, flt_t (*)(int, int));
void set_matrix_by_line(Matrix, flt_t *, char);
void set_matrix_row(Matrix, int, flt_t *);
void set_matrix_column(Matrix, int, flt_t *);
Matrix get_row_matrix(Matrix, int);
Matrix get_column_matrix(Matrix, int);
void clear_matrix(Matrix);
flt_t *flatten(Matrix, char);
void print_matrix(Matrix);
void print_row(Matrix, int);
Matrix add_matrix(Matrix, Matrix);
Matrix multiply_matrix(Matrix, Matrix);
Matrix scale_matrix(Matrix, flt_t);
Matrix matrix_power(Matrix, int);
Matrix transpose(Matrix);
Matrix inverse(Matrix);
Matrix append_horizontal(const Matrix, const Matrix);
void set_submatrix(Matrix, Matrix, int, int);
Matrix get_submatrix(Matrix, int, int, int, int);
int is_symmetric(Matrix);
int is_orthogonal(Matrix);
int is_upper_triangular(Matrix);

Vector new_vector(int);
void free_vector(Vector);
Vector copy_vector(Vector);
int vector_is_equal(Vector, Vector);
void set_vector_by_value(Vector, flt_t);
void set_vector_by_function(Vector, flt_t (*)(int));
void set_vector(Vector, flt_t *);
void clear_vector(Vector);
void print_vector(Vector);
Vector add_vector(Vector, Vector);
Vector multiply_matrix_vector(Matrix, Vector);
Vector scale_vector(Vector, flt_t);
flt_t norm(Vector);
flt_t normalize(Vector);
void normalize_columns(Matrix A);
flt_t dot_product(Vector, Vector);
Vector cross_product(Vector, Vector);

Matrix vector_to_row_matrix(Vector);
Matrix vector_to_column_matrix(Vector);
Vector matrix_to_vector(Matrix);
Vector get_row_vector(Matrix, int);
Vector get_column_vector(Matrix, int);
Matrix combine_row_vectors(Vector *, int);
Matrix combine_column_vectors(Vector *, int);
Vector *get_row_vectors(Matrix);
Vector *get_column_vectors(Matrix);

Vector std_unit_vector(int, int);
Matrix identity(int);
Matrix householder_reflection(const Vector v);
Matrix vector_projection_matrix(const Vector v);

#endif