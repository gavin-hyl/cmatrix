/**
 * @file cmatrix_basics.c
 * @author Gavin Hua (139950129+GavinHYL@users.noreply.github.com)
 * @brief Implementation of cmatrix_basics.h
 * 
 * @copyright Copyright (c) 2023
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cmatrix_basics.h"

/**
 * @brief Allocates memory for a new matrix.
 * 
 * @param r the rows of the matrix
 * @param c the columns of the matrix
 * @return Matrix: the new matrix
 */
Matrix new_matrix(int r, int c)
{
    Matrix A = (Matrix) malloc(sizeof(struct MaTrix));
    check_null(A);

    A->rows = r;
    A->cols = c;
    A->elements = (elem_t **) calloc(r, sizeof(elem_t *));
    for (int i = 0; i < r; i++)
    {
        A->elements[i] = (elem_t *) calloc(c, sizeof(elem_t));
    }
    return A;
}

/**
 * @brief Allocates memory for a new matrix, and copies the old matrix to it.
 * 
 * @param A the matrix to be copied
 * @return Matrix: the copy
 */
Matrix copy_matrix(Matrix A)
{
    int r=A->rows, c=A->cols;
    Matrix cpy = new_matrix(r, c);

    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            cpy->elements[i][j] = A->elements[i][j];
        }
    }
    return cpy;
}

/**
 * @brief Compares two matrices. 
 * 
 * @param A the first matrix.
 * @param B the second matrix.
 * @return int: 1 if they have the same size and same element values, 0 if not.
 */
int matrix_is_equal(Matrix A, Matrix B)
{
    check_null(A);
    check_null(B);
    if (A->rows!=B->rows || A->cols!=B->cols)
    {
        return 1;
    }
    int r=A->rows, c=A->cols;
    
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            if (!nearly_zero(A->elements[i][j] - B->elements[i][j]))
            {
                return 0;
            }
        }
    }
    return 1;
}

/**
 * @brief Set all elements of a matrix to a value.
 * 
 * @param A the matrix
 * @param val the value
 */
void set_matrix_by_value(Matrix A, elem_t val)
{
    check_null(A);
    int r=A->rows, c=A->cols;

    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            A->elements[i][j] = val;
        }
    }
}

/**
 * @brief Sets all elements of the matrix with a custom function.
 * 
 * @param A the matrix to be set
 * @param gen the function, takes two arguments (the row and column of the current element)
 */
void set_matrix_by_function(Matrix A, elem_t (*gen)(int, int))
{
    check_null(A);
    int r=A->rows, c=A->cols;
    
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            A->elements[i][j] = (*gen)(i, j);
        }
    }
}

/**
 * @brief Sets all elements of the matrix using an array of values.
 * 
 * @param A the matrix to be set
 * @param elements the elements
 * @param mode either ROW_MAJOR or COLUMN_MAJOR
 */
void set_matrix_by_line(Matrix A, elem_t *elements, char mode)
{
    check_null(A);
    int r=A->rows, c=A->cols;

    switch(mode)
    {
        case ROW_MAJOR:
            for (int i = 0; i < r; i++)
            {
                for (int j = 0; j < c; j++)
                {
                    A->elements[i][j] = *elements++;
                }
            }
            break;
        case COLUMN_MAJOR:
            for (int i = 0; i < c; i++)
            {
                for (int j = 0; j < r; j++)
                {
                    A->elements[j][i] = *elements++;
                }
            }
            break;
        default:
            perror("Unknown mode (set_by_line): must be BY_ROW or BY_COLUMN.\n");
            break;
    }
}

/**
 * @brief Sets elements of one row in a matrix using an array of values.
 * 
 * @param A the matrix to be set
 * @param r the row to be set (index from 0)
 * @param rowvals the values
 */
void set_matrix_row(Matrix A, int r, elem_t *rowvals)
{
    int c = A->cols;
    for (int i = 0; i < c; i++)
    {
        A->elements[r][i] = *rowvals++;
    }
}

/**
 * @brief Sets elements of one column in a matrix using an array of values.
 * 
 * @param A the matrix to be set
 * @param r the column to be set (index from 0)
 * @param colvals the values
 */
void set_matrix_column(Matrix A, int c, elem_t *colvals)
{
    int r = A->rows;
    for (int i = 0; i < r; i++)
    {
        A->elements[i][r] = *colvals++;
    }
}

/**
 * @brief Extracts a row matrix from a matrix.
 * 
 * @param A the matrix
 * @param c the row number (index from 0)
 * @return Matrix: the row matrix
 */
Matrix get_row_matrix(Matrix A, int r)
{
    int c = A->cols;
    Matrix row = new_matrix(1, c);

    for (int i = 0; i < c; i++)
    {
        row->elements[0][i] = A->elements[r][i];
    }
    return row;
}

/**
 * @brief Extracts a column matrix from a matrix.
 * 
 * @param A the matrix
 * @param c the column number (index from 0)
 * @return Matrix: the column matrix
 */
Matrix get_column_matrix(Matrix A, int c)
{
    int r = A->rows;
    Matrix col = new_matrix(r, 1);
    
    for (int i = 0; i < r; i++)
    {
        col->elements[i][0] = A->elements[i][c];
    }
    return col;
}

/**
 * @brief Sets all elements of the matrix to 0.
 * 
 * @param v the matrix to be cleared
 */
void clear_matrix(Matrix A)
{
    set_matrix_by_value(A, 0);
}

/**
 * @brief Extracts the elements of a matrix into a 1D element array.
 * 
 * @param A the matrix to be flattened
 * @param mode either ROW_MAJOR or COLUMN_MAJOR
 * @return elem_t*: the element array
 */
elem_t *flatten(Matrix A, char mode)
{
    int r=A->rows, c=A->cols;
    elem_t *flat = (elem_t *) calloc(r*c, sizeof(elem_t));

    switch(mode)
    {
        case ROW_MAJOR:
            for (int i = 0; i < r; i++)
            {
                for (int j = 0; j < c; j++)
                {
                    flat[i*c + j] = A->elements[i][j];
                }
            }
            break;
        case COLUMN_MAJOR:
            for (int i = 0; i < c; i++)
            {
                for (int j = 0; j < r; j++)
                {
                    flat[i*r + j] = A->elements[j][i];
                }
            }
            break;
        default:
            perror("[ERROR] Unknown mode (flatten): must be BY_ROW or BY_COLUMN.\n");
    }

    return flat;
}

/**
 * @brief Formats and prints a matrix to the standard output.
 * 
 * @param A the matrix to be printed
 */
void print_matrix(Matrix A)
{
    check_null(A);
    int r=A->rows, c=A->cols;
    if (r == 1)
    {
        fprintf(stdout, "( ");
        print_row(A, 0);
        fprintf(stdout, " )\n");
        return;
    }

    fprintf(stdout, "/ ");
    print_row(A, 0);
    fprintf(stdout, " \\\n");
    for (int i = 1; i < r-1; i++)
    {
        fprintf(stdout, "| ");
        print_row(A, i);
        fprintf(stdout, " |\n");
    }
    fprintf(stdout, "\\ ");
    print_row(A, r-1);
    fprintf(stdout, " /\n\n");
}

/**
 * @brief Formats and prints a row of a matrix to the standard output.
 * 
 * @param A the matrix from which a row is to be printed
 * @param r the row to be printed
 */
void print_row(Matrix A, int r)
{
    check_null(A);
    int c = A->cols;

    for (int i = 0; i < c; i++)
    {
        fprintf(stdout, "%*.*f", ELEMENT_PRINT_WIDTH, ELEMENT_PRINT_PRECISION, (float) A->elements[r][i]);
        fprintf(stdout, (i == c-1) ? "" : " ");
    }
}

/**
 * @brief Adds 2 matrices.
 * 
 * @param A the first matrix
 * @param B the second matrix
 * @return Matrix: the sum
 */
Matrix add_matrix(Matrix A, Matrix B)
{
    check_null(A);
    check_null(B);
    int r1=A->rows, c1=A->cols, r2=B->rows, c2=B->cols;
    
    if (r1!=r2 || c1!=c2)
    {
        fprintf(stdout, "ERROR (A_add): Arices must be same size.\n");
    }
    Matrix result = new_matrix(r1, c1);
    for (int i = 0; i < r1; i++)
    {
        for (int j = 0; j < c1; j++)
        {
            result->elements[i][j] = A->elements[i][j] + B->elements[i][j];
        }
    }
    return result;
}

/**
 * @brief Multiplies 2 matrices.
 * 
 * @param A the first matrix
 * @param B the second matrix
 * @return Matrix: the product
 */
Matrix multiply_matrix(Matrix A, Matrix B)
{
    check_null(A);
    check_null(B);
    check_multipliable(A, B);
    int r1=A->rows, c1=A->cols, r2=B->rows, c2=B->cols;
    
    Matrix result = new_matrix(r1, c2);
    elem_t temp;
    for (int i = 0; i < r1; i++)
    {   
        for (int j = 0; j < c2; j++)
        {   
            temp = 0;
            for (int k = 0; k < c1; k++)
            {
                temp += A->elements[i][k] * B->elements[k][j];
            }
            result->elements[i][j] = temp;
        }
    }
    return result;
}

/**
 * @brief Scales a matrix with a scalar
 * 
 * @param A the matrix to be scaled
 * @param n the scaling factor
 * @return Matrix: the scaled matrix
 */
Matrix scale_matrix(Matrix A, elem_t n)
{
    check_null(A);
    Matrix Acpy = copy_matrix(A);
    int r=A->rows, c=A->cols;

    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            Acpy->elements[i][j] *= n;
        }
    }
    return Acpy;
}

/**
 * @brief Raises a square matrix to an arbitrary integer power.
 * 
 * @param A the matrix base of the power
 * @param pow the power for A to be raised. If 0, then the identity is returned.
 * If <0, then we take the power of inv(A) instead.
 * @return Matrix: the resultant matrix
 */
Matrix matrix_power(Matrix A, int pow)
{
    check_square(A);
    Matrix result = identity(A->cols);

    if (pow < 0)
    {
        Matrix inv = inverse(A);
        while (pow-- > 0)
        {
            result = multiply_matrix(result, inv);
        }
        return result;
    }
    else
    {
        while (pow-- > 0)
        {
            result = multiply_matrix(result, A);
        }
        return result;
    }
}

/**
 * @brief Computes the transpose of a matrix.
 * 
 * @param A the matrix
 * @return Matrix: the transpose
 */
Matrix transpose(Matrix A)
{
    int r=A->rows, c=A->cols;
    Matrix result = new_matrix(c, r);

    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            result->elements[j][i] = A->elements[i][j];
        }
    }
    return result;
}

/**
 * @brief Computes the inverse of a square matrix.
 * 
 * @param A the matrix
 * @return Matrix: the inverse
 */
Matrix inverse(Matrix A)
{
    check_square(A);
    int d = A->rows;
    Matrix Acpy = copy_matrix(A);
    Matrix I = identity(d);

    for (int fd = 0; fd < d; fd++)
    {
        elem_t fd_scaler = 1 / Acpy->elements[fd][fd];
    
        for (int j = 0; j < d; j++)
        {
            Acpy->elements[fd][j] *= fd_scaler;
            I->elements[fd][j] *= fd_scaler;
        }

        for (int i = 0; i < fd; i++)
        {
            elem_t cr_scaler = Acpy->elements[i][fd];
            for (int j = 0; j < d; j++)
            {
                Acpy->elements[i][j] -= cr_scaler * Acpy->elements[fd][j];
                I->elements[i][j] -= cr_scaler * I->elements[fd][j];
            }
        }

        for (int i = fd+1; i < d; i++)
        {
            elem_t cr_scaler = Acpy->elements[i][fd];
            for (int j = 0; j < d; j++)
            {
                Acpy->elements[i][j] -= cr_scaler * Acpy->elements[fd][j];
                I->elements[i][j] -= cr_scaler * I->elements[fd][j];
            }
        }
    }

    return I;
}

/**
 * @brief Appends two matrices horizontally.
 * 
 * @param A the left matrix
 * @param B the right matrix
 * @return Matrix: the result
 */
Matrix append_horizontal(Matrix A, Matrix B)
{
    int Acols=A->cols, Arows=A->rows, Bcols=B->cols, Brows=B->rows;
    int Mrows=Arows, Mcols=Acols+Bcols;
    Matrix M = new_matrix(Mrows, Mcols);

    for (int i = 0; i < Mrows; i++)
    {
        for (int j = 0; j < Acols; j++)
        {
            M->elements[i][j] = A->elements[i][j];
        }
        for (int j = Acols; j < Mcols; j++)
        {
            M->elements[i][j] = B->elements[i][j-Acols];
        }
    }

    return M;
}

void set_submatrix(Matrix A, Matrix a, int r1, int c1)
{
    int r=a->rows, c=a->cols, R=A->rows, C=A->cols;
    if (r1+r>R || c1+c>C)
    {
        fprintf(stderr, "Boundaries of A exceeded.\n");
        exit(1);
    }
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            A->elements[r1+i][c1+j] = a->elements[i][j];
        }
    }
}

/**
 * @brief Get the submatrix object
 * 
 * @param A 
 * @param r1 
 * @param c1 
 * @param r2 
 * @param c2 
 * @return Matrix 
 */
Matrix get_submatrix(Matrix A, int r1, int c1, int r2, int c2)
{
    int R=A->rows, C=A->cols, r=r2-r1, c=c2-c1;
    Matrix submatrix = new_matrix(r, c);
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            submatrix->elements[i][j] = A->elements[r1+i][c1+j];
        }
    }
    return submatrix;
}

/**
 * @brief Checks if a matrix is orthogonal
 * 
 * @param A the matrix
 * @return int: 1 if A is orthogonal, 0 if not
 */
int is_orthogonal(Matrix A)
{
    if (A->rows != A->cols)
    {
        return 0;
    }
    Matrix At = transpose(A);
    return matrix_is_equal(identity(A->cols), multiply_matrix(A, At));
}

/**
 * @brief Checks if a matrix is symmetric
 * 
 * @param A the matrix
 * @return int: 1 if A is symmetric, 0 if not
 */
int is_symmetric(Matrix A)
{
    if (A->rows != A->cols)
    {
        return 0;
    }
    return matrix_is_equal(transpose(A), A);
}


/*
--------------------------- VECTOR FUNCTIONS BELOW -----------------------------
*/

/**
 * @brief Allocates memory for a new vector.
 * 
 * @param d the dimension of the new vector
 * @return Vector: the new vector
 */
Vector new_vector(int d)
{
    Vector v = (Vector) malloc(sizeof(struct VecTor));
    check_null(v);

    v->dim = d;
    v->elements = (elem_t *) calloc(d, sizeof(elem_t));
    return v;
}

/**
 * @brief Allocates memory for a new vector, and copies the old vector to it.
 * 
 * @param A the vector to be copied
 * @return Vector: the copy
 */
Vector copy_vector(Vector v)
{
    int d = v->dim;
    Vector cpy = new_vector(d);
    
    for (int i = 0; i < d; i++)
    {
        cpy->elements[i] = v->elements[i];
    }
    return cpy;
}

/**
 * @brief Compares two vectors.
 * 
 * @param v the first vector
 * @param w the second vector
 * @return int: 1 if the vectors have equal dimension and elements, 0 if not
 */
int vector_is_equal(Vector v, Vector w)
{
    if (v->dim != w->dim)
    {
        return 0;
    }
    int d = v->dim;
    
    for (int i = 0; i < d; i++)
    {
        if (!nearly_zero(v->elements[i] - w->elements[i]))
        {
            return 0;
        }
    }
    return 1;
}

/**
 * @brief Sets all elements of a vector to be a constant value.
 * 
 * @param v the vector to be set
 * @param e the value
 */
void set_vector_by_value(Vector v, elem_t e)
{
    int d = v->dim;
    for (int i = 0; i < d; i++)
    {
        v->elements[i] = e;
    }
}

/**
 * @brief Sets all elements of the vector with a custom function.
 * 
 * @param v the vector to be set
 * @param gen the function, takes one argument (the index of the current element)
 */
void set_vector_by_function(Vector v, elem_t (*gen)(int))
{
    check_null(v);
    int d = v->dim;

    for (int i = 0; i < d; i++)
    {
        v->elements[i] = gen(i);
    }
}

/**
 * @brief Sets elements of the vector with an array
 * 
 * @param v the vector to be set
 * @param elements the element array
 */
void set_vector(Vector v, elem_t *elements)
{
    check_null(v);
    int d = v->dim;

    for (int i = 0; i < d; i++)
    {
        v->elements[i] = *elements++;
    }
}

/**
 * @brief Sets all elements of the vector to 0.
 * 
 * @param v the vector to be cleared
 */
void clear_vector(Vector v)
{
    set_vector_by_value(v, 0);
}

/**
 * @brief Formats and prints a vector onto the standard output.
 * 
 * @param v the vector to be printed
 */
void print_vector(Vector v)
{
    check_null(v);
    int d = v->dim;

    if (d == 1)
    {
        fprintf(stdout, "(%*.*f)\n", ELEMENT_PRINT_WIDTH, ELEMENT_PRINT_PRECISION, v->elements[0]);
        return;
    }
    fprintf(stdout, "/%*.*f\\\n", ELEMENT_PRINT_WIDTH, ELEMENT_PRINT_PRECISION, v->elements[0]);
    for (int i = 1; i < d-1; i++)
    {
        fprintf(stdout, "|%*.*f|\n", ELEMENT_PRINT_WIDTH, ELEMENT_PRINT_PRECISION, v->elements[i]);
    }
    fprintf(stdout, "\\%*.*f/\n", ELEMENT_PRINT_WIDTH, ELEMENT_PRINT_PRECISION, v->elements[d-1]);
}

/**
 * @brief Adds 2 vectors.
 * 
 * @param v the first vector
 * @param w the second vector
 * @return the sum
 */
Vector add_vector(Vector v, Vector w)
{
    check_equal_dimension(v, w);
    int d = v->dim;
    Vector result = new_vector(d);

    while (d-- > 0)
    {
        result->elements[d] = v->elements[d] + w->elements[d];
    }
    return result;
}

/**
 * @brief Multiplies a matrix by a vector, in the form Ax.
 * 
 * @param A the matrix to be multiplied
 * @param v the vector to be multipled
 * @return the product
 */
Vector multiply_matrix_vector(Matrix A, Vector v)
{
    check_null(A);
    check_null(v);
    return matrix_to_vector(multiply_matrix(A, vector_to_column_matrix(v)));
}

/**
 * @brief Multiplies a vector by a scalar constant in place.
 * 
 * @param v the vector to be multiplied
 * @param k the scaling factor
 * @return Vector: the scaled vector
 */
Vector scale_vector(Vector v, elem_t k)
{
    check_null(v);
    int d = v->dim;
    Vector vcpy = copy_vector(v);
    while (d-- > 0)
    {
        vcpy->elements[d] *= k;
    }
    return vcpy;
}

/**
 * @brief Computes the length (norm) of a vector
 * 
 * @param v the vector
 * @return the length
 */
elem_t norm(Vector v)
{
    check_null(v);
    int d = v->dim;

    elem_t len_squared = 0;
    for (int i = 0; i < d; i++)
    {
        len_squared += v->elements[i] * v->elements[i];
    }
    return sqrt(len_squared);
}

/**
 * @brief Normalizes a vector in place.
 * 
 * @param v the vector
 * @return elem_t: the length
 */
elem_t normalize(Vector v)
{
    check_null(v);
    int d = v->dim;
    elem_t length = norm(v);

    for (int i = 0; i < d; i++)
    {
        v->elements[i] /= length;
    }
    return length;
}

/**
 * @brief Normalizes the columns of a matrix in place
 * 
 * @param A the matrix
*/
void normalize_columns(Matrix A)
{
    int c=A->cols;
    Vector col = NULL;
    for (int i = 0; i < c; i++)
    {
        col = get_column_vector(A, 0);
        normalize(col);
        set_matrix_column(A, i, col->elements);
    }
}

/**
 * @brief Computes the dot product between two vectors.
 * 
 * @param v the first vector
 * @param w the second vector
 * @return the dot product
 */
elem_t dot_product(Vector v, Vector w)
{
    check_equal_dimension(v, w);
    int d = v->dim;
    elem_t result = 0;

    for (int i = 0; i < d; i++)
    {
        result += v->elements[i] * w->elements[i];
    }
    return result;
}

/**
 * @brief Computes the cross product of two vectors of either 2 or 3 elements.
 * 
 * @param v the first vector, will be promoted to 3D if passed in 2D vector
 * @param w the second vector, will be promoted to 3D if passed in 2D vector
 * @return the crossed vector
 */
Vector cross_product(Vector v, Vector w)
{
    check_null(v);
    check_null(w);
    check_equal_dimension(v, w);
    Vector result = new_vector(3);
    elem_t *ve=v->elements, *we=w->elements;
    
    switch(v->dim)
    {
        case 2:
            result->elements[2] = ve[0]*we[1] - ve[1]*we[0];
            return result;
        case 3:
            result->elements[0] = ve[1]*we[2] - ve[2]*we[1];
            result->elements[1] = ve[2]*we[0] - ve[0]*we[2];
            result->elements[2] = ve[0]*we[1] - ve[1]*we[0];
            return result;
        default:
            fprintf(stderr, "Cross product only defined for 2D and 3D.\n");
            exit(1);
    }
}

/**
 * @brief Converts a vector to a row matrix.
 * 
 * @param v the vector to be converted
 * @return the corresponding matrix
 */
Matrix vector_to_row_matrix(Vector v)
{
    int d = v->dim;
    Matrix row = new_matrix(1, d);

    for (int i = 0; i < d; i++)
    {
        row->elements[0][i] = v->elements[i];
    }
    return row;
}

/**
 * @brief Converts a vector to a column matrix.
 * 
 * @param v the vector to be converted
 * @return the corresponding matrix
 */
Matrix vector_to_column_matrix(Vector v)
{
    int d = v->dim;
    Matrix column = new_matrix(d, 1);

    for (int i = 0; i < d; i++)
    {
        column->elements[i][0] = v->elements[i];
    }
    return column;
}

/**
 * @brief Converts a matrix to a vector.
 * 
 * @param A the matrix to be converted, can be either 1xn or nx1
 * @return the corresponding vector
 */
Vector matrix_to_vector(Matrix A)
{
    check_null(A);
    int r=A->rows, c=A->cols;
    if (c == 1)
    {
        Vector v = new_vector(r);
        for (int i = 0; i < r; i++)
        {
            v->elements[i] = A->elements[i][0];
        }
        return v;
    }
    else
    {
        Vector v = new_vector(c);
        for (int i = 0; i < c; i++)
        {
            v->elements[i] = A->elements[0][c];
        }
        return v;
    }
}

/**
 * @brief Extracts a row vector from a matrix
 * 
 * @param A the matrix
 * @param r the row number (index from 0)
 * @return the row vector
 */
Vector get_row_vector(Matrix A, int r)
{
    check_null(A);
    int c = A->cols;
    Vector row = new_vector(c);
    
    for (int i = 0; i < c; i++)
    {
        row->elements[i] = A->elements[r][i];
    }
    return row;
}

/**
 * @brief Extracts a column vector from a matrix.
 * 
 * @param A the matrix
 * @param c the column number (index from 0)
 * @return the column vector
 */
Vector get_column_vector(Matrix A, int c)
{
    check_null(A);
    int r = A->rows;
    Vector col = new_vector(r);
    
    for (int i = 0; i < r; i++)
    {
        col->elements[i] = A->elements[i][c];
    }
    return col;
}

Matrix combine_row_vectors(Vector *vecs, int nvecs)
{
    Matrix combined = new_matrix(nvecs, vecs[0]->dim);
    for (int i = 0; i < nvecs; i++)
    {
        set_matrix_row(combined, i, vecs[i]->elements);
    }
    return combined;
}

Matrix combine_column_vectors(Vector *vecs, int nvecs)
{
    Matrix combined = new_matrix(vecs[0]->dim, nvecs);
    for (int i = 0; i < nvecs; i++)
    {
        set_matrix_column(combined, i, vecs[i]->elements);
    }
    return combined;
}

Vector *get_row_vectors(Matrix A)
{
    int r = A->rows;
    Vector *vecs = (Vector *) calloc(r, sizeof(Vector));
    for (int i = 0; i < r; i++)
    {
        vecs[i] = get_row_vector(A, i);
    }
    return vecs;
}

Vector *get_column_vectors(Matrix A)
{
    int c = A->cols;
    Vector *vecs = (Vector *) calloc(c, sizeof(Vector));
    for (int i = 0; i < c; i++)
    {
        vecs[i] = get_row_vector(A, i);
    }
    return vecs;
}

/**
 * @brief Generates an identity matrix.
 * 
 * @param d the size of the identity
 * @return the identity
 */
Matrix identity(int d)
{
    Matrix iden = new_matrix(d, d);
    for (int i = 0; i < d; i++)
    {
        iden->elements[i][i] = 1;
    }
    return iden;
}

/**
 * @brief Generates a standard unit vector in n-D space. 
 * 
 * @param dim the dimension of the unit vector
 * @param n the index (from 0) of the unit vector, ie. which column of the identity it is
 * @return the unit vector
 */
Vector std_unit_vector(int dim, int n)
{
    if (n > dim)
    {
        fprintf(stderr, "The specified unit vector does not exist (n > dim).\n");
        exit(1);
    }
    Vector b = new_vector(dim);
    b->elements[n] = 1;
    return b;
}

/**
 * @brief Computes the projection matrix of the space onto a vector. 
 * 
 * @param v the vector to be projected on to
 * @return the projection matrix
 */
Matrix vector_projection_matrix(const Vector v)
{
    check_null(v);
    elem_t length = norm(v);
    if (length == 0)
    {
        return identity(v->dim);
    }
    Matrix proj_v = multiply_matrix(vector_to_column_matrix(v), vector_to_row_matrix(v));
    proj_v = scale_matrix(proj_v, 1 / dot_product(v, v));
    return proj_v;
}

/**
 * @brief Computes the Householder reflection matrix of the n-D space about a
 * vector, ie. rotating the space 180 degrees around that vector.
 * 
 * @param v the vector about which the space is reflected
 * @return the reflection matrix
 */
Matrix householder_reflection(Vector v)
{
    check_null(v);
    Matrix mir = vector_projection_matrix(v);
    mir = scale_matrix(mir, -2);
    Matrix I = identity(v->dim);
    return add_matrix(mir, I);
}