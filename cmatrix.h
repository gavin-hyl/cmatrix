/**
 * CMatrix - a lightweight linear algebra library for embedded systems.
 * 
 * Includes functions and types for creating and manipulating matrices, as well as various
 * linear algebra operations. Also includes a set of generators for common
 * matrices.
 * 
 * Author: Gavin Hua (ghua@caltech.edu)
 * 
*/

#define EPSILON 0.01
#define ELEMENTWIDTH 7  // for printing the matrix, has no influence on calculations.
#define ELEMENTPREC 2

typedef double elem_t;

typedef struct MaTrix {
    int rows;
    int cols;
    elem_t **elements;
} *Matrix;

typedef struct VecTor {
    int dim;
    elem_t *elements;
} *Vector;

enum {LEFT, RIGHT, UP, DOWN, UR, UL, DL, DR};
enum {BY_ROW, BY_COLUMN};

Matrix new_matrix(int, int);
Matrix copy_matrix(Matrix);
Vector new_vector(int);
Vector copy_vector(Vector);
int compare_matrix(Matrix, Matrix);
int compare_vector(Vector, Vector);
void set_matrix_by_value(Matrix, elem_t);
void set_vector_by_value(Vector, elem_t);
void set_matrix_by_function(Matrix, elem_t (*)(int, int));
void set_vector_by_function(Vector, elem_t (*)(int));
void set_matrix_by_line(Matrix, elem_t *, char);
void set_matrix_row(Matrix, int, elem_t *);
void set_matrix_column(Matrix, int, elem_t *);
void clear_matrix(Matrix);
void clear_vector(Vector);

Matrix vector_to_row_matrix(Vector);
Matrix vector_to_column_matrix(Vector);
Vector matrix_to_vector(Matrix);

Matrix get_row_matrix(Matrix, int);
Matrix get_column_matrix(Matrix, int);
Vector get_row_vector(Matrix, int);
Vector get_column_vector(Matrix, int);
elem_t *flatten(Matrix, char);
void print_matrix(Matrix);
void print_row(Matrix, int);
void print_vector(Vector);

void swap_row(Matrix, int, int);
void multiply_row(Matrix, int, elem_t);
void combine_row(Matrix, int, int);

Matrix add_matrix(Matrix, Matrix);
Vector add_vector(Vector, Vector);
Matrix multiply_matrix(Matrix, Matrix);
Vector multiply_matrix_vector(Matrix, Vector);
void scalar_multiply_matrix(Matrix, elem_t);
void scalar_multiply_vector(Vector, elem_t);
Matrix transpose(Matrix);
Matrix append(const Matrix, const Matrix, const char);

elem_t row_echelon_form(Matrix);
elem_t determinant(Matrix);
Matrix inverse(Matrix);
Matrix solve_linear_system(Matrix, Vector);
Matrix basis_transform(Matrix, Matrix);
// unfinished.
Matrix *qr_decomposition(Matrix);
Matrix eigenvectors(Matrix);
elem_t *eigenvalues(Matrix);

Matrix identity(int);
Matrix row_swap_matrix(int, int, int);
Matrix row_combine_matrix(int, int, int);
Matrix row_multiply_matrix(int, int, elem_t);
Matrix vector_projection_matrix(Vector);

Matrix M2_rotation(elem_t);
Matrix M2_flip(Matrix);
Matrix M3_rotation(elem_t, elem_t);
Matrix M3_axis_rotation(elem_t, Vector);
Matrix M3_flip(Matrix, Matrix);

Vector standard_unit_vector(int, int);
void normalize(Vector);
elem_t vector_length(Vector);
elem_t dot_product(Vector, Vector);
Vector cross_product(Vector, Vector);