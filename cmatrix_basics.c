#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cmatrix_basics.h"

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

int compare_matrix(Matrix A, Matrix B)
{
    check_null(A);
    check_null(B);
    check_equal_size(A, B);
    int r=A->rows, c=A->cols;
    
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            if (!nearly_zero(A->elements[i][j] - B->elements[i][j]))
            {
                return 1;
            }
        }
    }
    return 0;
}

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

void set_matrix_row(Matrix A, int r, elem_t *rowvals)
{
    int c = A->cols;
    for (int i = 0; i < c; i++)
    {
        A->elements[r][i] = *rowvals++;
    }
}

void set_matrix_column(Matrix A, int c, elem_t *colvals)
{
    int r = A->rows;
    for (int i = 0; i < r; i++)
    {
        A->elements[i][r] = *colvals++;
    }
}


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

void clear_matrix(Matrix A)
{
    set_matrix_by_value(A, 0);
}

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


void print_matrix(Matrix A)
{
    check_null(A);
    int r=A->rows, c=A->cols;
    if (r == 1)
    {
        printf("( ");
        print_row(A, 0);
        printf(" )\n");
        return;
    }

    printf("/ ");
    print_row(A, 0);
    printf(" \\\n");
    for (int i = 1; i < r-1; i++)
    {
        printf("| ");
        print_row(A, i);
        printf(" |\n");
    }
    printf("\\ ");
    print_row(A, r-1);
    printf(" /\n\n");
}

void print_row(Matrix A, int r)
{
    check_null(A);
    int c = A->cols;

    for (int i = 0; i < c; i++)
    {
        printf("%*.*f", ELEMENTWIDTH, ELEMENTPREC, (float) A->elements[r][i]);
        printf((i == c-1) ? "" : " ");
    }
}

Matrix add_matrix(Matrix A, Matrix B)
{
    check_null(A);
    check_null(B);
    int r1=A->rows, c1=A->cols, r2=B->rows, c2=B->cols;
    
    if (r1!=r2 || c1!=c2)
    {
        printf("ERROR (A_add): Arices must be same size.\n");
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
 * @brief Scales a matrix with a scalar in place.
 * 
 * @param A the matrix to be scaled
 * @param n the scaling factor
 */
void scale_matrix(Matrix A, elem_t n)
{
    check_null(A);
    int r=A->rows, c=A->cols;

    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            A->elements[i][j] *= n;
        }
    }
}

/**
 * @brief Raises a square matrix to an arbitrary integer pattern.
 * 
 * @param A the base of the power
 * @param pow the power for A to be raised. If 0, then the identity is returned.
 * If <0, then we take the power of inv(A) instead.
 * @return the resultant matrix
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
 * @brief Transposes a matrix.
 * 
 * @param A the matrix to be trasnposed
 * @return A^T
 */
Matrix transpose(const Matrix A)
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
 * @brief Inverts a square matrix.
 * 
 * @param A the matrix to be inverted
 * @return the inverted matrix
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

Matrix append_horizontal(const Matrix A, const Matrix B)
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

Vector new_vector(int d)
{
    Vector v = (Vector) malloc(sizeof(struct VecTor));
    check_null(v);

    v->dim = d;
    v->elements = (elem_t *) calloc(d, sizeof(elem_t));
    return v;
}

Vector copy_vector(Vector v)
{
    int d = v->dim;
    Vector cpy = new_vector(d);
    
    for (int i = 0; i < d; i++)
    {
        cpy->elements[i] = v->elements[i];
    }
}

int compare_vector(Vector v, Vector w)
{
    check_equal_dimension(v, w);
    int d = v->dim;
    
    for (int i = 0; i < d; i++)
    {
        if (!nearly_zero(v->elements[d] - w->elements[d]))
        {
            return 1;
        }
    }
    return 0;
}

void set_vector_by_value(Vector v, const elem_t e)
{
    int d = v->dim;
    for (int i = 0; i < d; i++)
    {
        v->elements[i] = e;
    }
}

void set_vector_by_function(Vector v, elem_t (*gen)(int))
{
    check_null(v);
    int d = v->dim;

    for (int i = 0; i < d; i++)
    {
        v->elements[i] = gen(i);
    }
}

void set_vector(Vector v, elem_t *elements)
{
    check_null(v);
    int d = v->dim;

    for (int i = 0; i < d; i++)
    {
        v->elements[i] = *elements++;
    }
}

void clear_vector(Vector v)
{
    set_vector_by_value(v, 0);
}

void print_vector(Vector v)
{
    check_null(v);
    int d = v->dim;

    if (d == 1)
    {
        printf("(%*.*f)\n", ELEMENTWIDTH, ELEMENTPREC, v->elements[0]);
        return;
    }
    printf("/%*.*f\\\n", ELEMENTWIDTH, ELEMENTPREC, v->elements[0]);
    for (int i = 1; i < d-1; i++)
    {
        printf("|%*.*f|\n", ELEMENTWIDTH, ELEMENTPREC, v->elements[i]);
    }
    printf("\\%*.*f/\n", ELEMENTWIDTH, ELEMENTPREC, v->elements[d-1]);
}

Vector add_vector(Vector v, Vector w)
{
    check_equal_dimension(v, w);
    int d = v->dim;
    Vector result = new_vector(d);

    while (d-- > 0)
    {
        result->elements[d] = v->elements[d] + w->elements[d];
    }
}

Vector multiply_matrix_vector(Matrix A, Vector v)
{
    check_null(A);
    check_null(v);
    return matrix_to_vector(multiply_matrix(A, vector_to_column_matrix(v)));
}

void scale_vector(Vector v, elem_t k)
{
    check_null(v);
    int d = v->dim;

    while (d-- > 0)
    {
        v->elements[d] *= k;
    }
}

elem_t vector_length(const Vector v)
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
 * @param v the vector to be normalized
 */
void normalize(Vector v)
{
    check_null(v);
    int d = v->dim;
    elem_t length = vector_length(v);

    for (int i = 0; i < d; i++)
    {
        v->elements[i] /= length;
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

Vector matrix_to_vector(Matrix A)
{
    check_null(A);
    int r=A->rows;
    Vector v = new_vector(r);

    for (int i = 0; i < r; i++)
    {
        v->elements[i] = A->elements[i][0];
    }
    return v;
}

Vector get_row_vector(Matrix A, int r)
{
    check_null(A);
    int c = A->cols;
    Vector col = new_vector(c);
    
    for (int i = 0; i < c; i++)
    {
        col->elements[i] = A->elements[i][r];
    }
    return col;
}

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
 * @param n the index of the unit vector, ie. which column of the identity it is. 
 * @return the unit vector.
 */
Vector std_unit_vector(int dim, int n)
{
    if (n > dim)
    {
        fprintf(stderr, "The specified unit vector does not exist (n > dim).\n");
        exit(1);
    }
    Vector b = new_vector(dim);
    b->elements[n-1] = 1;
    return b;
}