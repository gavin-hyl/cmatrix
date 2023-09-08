# cmatrix - C-based Linear Algebra
A lightweight C-based linear algebra library. Work in progress, please fork/raise an issue if you'd like to contribute!

# File Structure (How to Use)

This is intended as an overview, ***not an exhaustive list***. please refer to the documentation for more detailed descriptions of the functions included in this library.
## `cmatrix.h`
Includes everything that the library has to offer, essentially a one-stop approach for using this library. 
## `cmatrix_defs.h`
Includes type definitions for vectors/matrices and basic error checking macros.

### `flt_t`
The floating point type to be used in this library. Only `float` and `double` are allowed, `long double` is not supported.

### `Matrix`
Whenever we declare a matrix, we declare a pointer to a struct that represents a generic matrix. The struct itself contains the number of rows, number of columns, and the element array. Storage allocation is taken care of by the `new_matrix` function.
```
typedef struct MaTrix {
    int rows;
    int cols;
    flt_t **elements;
} *Matrix;
```
Access an element on row `i` and column `j` as such:
```
m->elements[i][j]
```

### `ROW_MAJOR` and `COLUMN_MAJOR`
These two variables are used in modes for setting and getting matrix elements. 
### `Vector`
Due to the pervasivenss of vectors in linear algebra, a seperate object is defined specially for vectors. The struct itself contains the dimension of the vector and the element array.
```
typedef struct VecTor {
    int dim;
    flt_t *elements;
} *Vector;
```
Access the `i`'th element as such:
```
v->elements[i]
```
Vectors and matrices may be converted to one another, see `cmatrix_basics.h`.

## `cmatrix_basics.h`
Includes all the basic opertions needed for `Matrix` and `Vector` creation, conversion, manipulation, and property checks.

### Creation
Creating and freeing `Matrix` objects:
```
Matrix new_matrix(int, int);
void free_matrix(Matrix);
Matrix copy_matrix(Matrix);
```
Creating and freeing `Vector` objects:
```
Vector new_vector(int);
void free_vector(Vector);
Vector copy_vector(Vector);
```
Common matrices and vectors:
```
#define ex2 std_unit_vector(2, 0)
#define ey2 std_unit_vector(2, 1)
#define ex3 std_unit_vector(3, 0)
#define ey3 std_unit_vector(3, 1)
#define ez3 std_unit_vector(3, 2)
// ...
Vector std_unit_vector(int, int);
Matrix identity(int);
Matrix householder_reflection(Vector v);
Matrix vector_projection_matrix(Vector v);
```

### Conversion
Converting from `Matrix` to `Vector`:
```
Vector matrix_to_vector(Matrix);
Vector get_row_vector(Matrix, int);
Vector get_column_vector(Matrix, int);
Vector *get_row_vectors(Matrix);
Vector *get_column_vectors(Matrix);
```
Converting from `Vector` to `Matrix`:
```
Matrix vector_to_row_matrix(Vector);
Matrix vector_to_column_matrix(Vector);
Matrix combine_row_vectors(Vector *, int);
Matrix combine_column_vectors(Vector *, int);
```
It should be noted that when passing in an argument for `Vector *`, it should be casted first to `(Vector [])`, as such:
```
Matrix M = combine_row_vectors((Vector []) {ex2, ex3}, 2)
```

### Manipulation
We can set the values of the `elements` array in `Matrix` objects in bulk by using the following functions:
```
void set_matrix_by_value(Matrix, flt_t);
void set_matrix_by_function(Matrix, flt_t (*)(int, int));
void set_matrix_by_line(Matrix, flt_t *, char);
void set_matrix_row(Matrix, int, flt_t *);
void set_matrix_column(Matrix, int, flt_t *);
```

Note that for the operations listed below, all will return a new `Matrix` or `Vector`, if applicable.

For `Matrix` objects, some common operations are
```
Matrix add_matrix(Matrix, Matrix);
Matrix multiply_matrix(Matrix, Matrix);
Matrix scale_matrix(Matrix, flt_t);
Matrix matrix_power(Matrix, int);
Matrix transpose(Matrix);
Matrix inverse(Matrix);
```

For `Vector` objects, some common operations are
```
Vector add_vector(Vector, Vector);
Vector multiply_matrix_vector(Matrix, Vector);
Vector scale_vector(Vector, flt_t);
flt_t norm(Vector);
flt_t normalize(Vector);
flt_t dot_product(Vector, Vector);
Vector cross_product(Vector, Vector);
```
### Property checks
Mainly for `Matrix` objects. Returns 1 for true.
```
int is_symmetric(Matrix);
int is_orthogonal(Matrix);
int is_upper_triangular(Matrix);
```

## `cmatrix_algebra.h`
Includes more advanced algebraic opertions, such as constructing Gram-Schmidt orthogonal bases and performing QR decompositions.
```
flt_t row_echelon_form(Matrix);
flt_t determinant(Matrix);
Vector solve_linear_system(Matrix, Vector);
Matrix basis_transform(Matrix, Matrix);
int linear_independence(Vector *, int);
Matrix GS_orthogonal_basis(Vector *, int);
Matrix *QR_decomposition(Matrix);
Matrix QR_decomposition_R_only(Matrix);
Vector eigenvalues(Matrix);
```

## `cmatrix_geometry.h`
Includes generators for various rotation, reflection, and projection matrices in 2D and 3D. Also includes functions for coordinate system transformations.
```
Matrix rotation_matrix_2d(flt_t);
Matrix rotation_matrix_3d(Vector, flt_t);
Matrix reflect_about_plane(Vector, Vector);
Matrix reflect_about_normal(Vector);
Vector unit_normal(Vector, Vector);
Vector polar_to_cartesian(flt_t, flt_t);
Vector spherical_to_cartesian(flt_t, flt_t, flt_t);
Vector cylindrical_to_cartesian(flt_t, flt_t, flt_t);
Matrix attitude_description(char *, flt_t *);
```
