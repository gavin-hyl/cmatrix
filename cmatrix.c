#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#include "cmatrix.h"

#define nearly_zero(e) (e<EPSILON && e>(-1)*EPSILON)

#define check_null(A) \
({  \
    if (A == NULL)  \
    {   \
        fprintf(stderr, "Cannot operate on NULL matrix.\n");    \
        exit(1);    \
    }   \
})

#define check_equal_size(A, B)  \
({  \
    if (A->rows!=B->rows || A->cols!=B->cols)   \
    {   \
        fprintf(stderr, "Cannot operate on two matrices of different size.\n"); \
    }   \
})

#define check_equal_dimension(v, w) \
({  \
    if (v->dim != w->dim)   \
    {   \
        fprintf(stderr, "Cannot operate on two vectors of different dimension.\n"); \
        exit(1);    \
    }   \
})


#define check_multipliable(A, B)    \
({  \
    if (A->rows != B->cols) \
    {   \
        fprintf(stderr, "Cannot multiply when r1 is not equal to c2.\n");   \
        exit(1);    \
    }   \
})

#define check_square(A) \
({  \
    if (A->rows != A->cols) \
    {   \
        fprintf(stderr, "Cannot operate on non-square matrix.\n");  \
        exit(1);    \
    }   \
})

Matrix identity(int d)
{
    Matrix iden = new_matrix(d, d);
    for (int i = 0; i < d; i++)
    {
        iden->elements[i][i] = 1;
    }
    return iden;
}

Matrix row_combine_matrix(int d, int fr, int to)
{
    Matrix combrow = identity(d);
    combrow->elements[to][fr] = 1;
    return combrow;
}

Matrix row_swap_matrix(int d, int r1, int r2)
{
    Matrix swaprow = identity(d);
    swaprow->elements[r1][r1] = 0;
    swaprow->elements[r2][r2] = 0;
    swaprow->elements[r1][r2] = 1;
    swaprow->elements[r2][r1] = 1;
    return swaprow;
}

Matrix row_multiply_matrix(int d, int r, elem_t mult)
{
    Matrix multrow = identity(d);
    multrow->elements[r][r] = mult;
    return multrow;
}

/**
 * A 2D counterclockwise rotation centered at the origin.
*/
Matrix M2_rotation(elem_t angle)
{
    Matrix mat = new_matrix(2, 2);
    
    mat->elements[0][0] = cos(angle);
    mat->elements[1][0] = sin(angle);
    mat->elements[0][1] = -sin(angle);
    mat->elements[1][1] = cos(angle);
    return mat;
}

Matrix M3_rotation(elem_t angle1, elem_t angle2)
{

}

/**
 * A projection matrix of any vector onto v, based on the formula (v * vT) / ||v||^2.
*/
Matrix vector_projection_matrix(const Vector v)
{
    Matrix proj_v = multiply_matrix(vector_to_column_matrix(v), vector_to_row_matrix(v));
    scalar_multiply_matrix(proj_v, 1 / dot_product(v, v));
    return proj_v;
}


void set_vector_by_value(Vector v, const elem_t e)
{
    int d = v->dim;
    for (int i = 0; i < d; i++)
    {
        v->elements[i] = e;
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

Vector new_vector(int d)
{
    Vector v = (Vector) malloc(sizeof(struct VecTor));
    check_null(v);

    v->dim = d;
    v->elements = (elem_t *) calloc(d, sizeof(elem_t));
    return v;
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

Vector copy_vector(Vector v)
{
    int d = v->dim;
    Vector cpy = new_vector(d);
    
    for (int i = 0; i < d; i++)
    {
        cpy->elements[i] = v->elements[i];
    }
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
        case BY_ROW:
            for (int i = 0; i < r; i++)
            {
                for (int j = 0; j < c; j++)
                {
                    A->elements[i][j] = *elements++;
                }
            }
            break;
        case BY_COLUMN:
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
            if (nearly_zero(A->elements[i][j] - B->elements[i][j]))
            {
                return 1;
            }
        }
    }
    return 0;
}

int compare_vector(Vector v, Vector w)
{
    check_null(v);
    check_null(w);
    check_equal_dimension(v, w);
    int d = v->dim;
    
    for (int i = 0; i < d; i++)
    {
        if (nearly_zero(v->elements[d] - w->elements[d]))
        {
            return 1;
        }
    }
    return 0;
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

void clear(Matrix A)
{
    check_null(A);
    set_matrix_by_value(A, 0);
}

void set_matrix_row(Matrix A, int r, elem_t *rowvals)
{
    int c = A->cols;
    for (int i = 0; i < c; i++)
    {
        A->elements[r][i] = *rowvals++;
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

void set_matrix_column(Matrix A, int c, elem_t *colvals)
{
    int r = A->rows;
    for (int i = 0; i < r; i++)
    {
        A->elements[i][r] = *colvals++;
    }
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

elem_t *flatten(Matrix A, char mode)
{
    int r=A->rows, c=A->cols;
    elem_t *flat = (elem_t *) calloc(r*c, sizeof(elem_t));

    switch(mode)
    {
        case BY_ROW:
            for (int i = 0; i < r; i++)
            {
                for (int j = 0; j < c; j++)
                {
                    flat[i*c + j] = A->elements[i][j];
                }
            }
            break;
        case BY_COLUMN:
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

void swap_row(Matrix A, int r1, int r2)
{
    if (r1 != r2)
    {
        elem_t *temp = get_row_matrix(A, r1)->elements[0];
        set_matrix_row(A, r1, get_row_matrix(A, r2)->elements[0]);
        set_matrix_row(A, r2, temp);
    }
}

void multiply_row(Matrix A, int r, elem_t k)
{
    if (!nearly_zero(k-1))
    {
        Matrix row = get_row_matrix(A, r);
        scalar_multiply_matrix(row, k);
        set_matrix_row(A, r, row->elements[0]);
    }
}

void combine_row(Matrix A, int fr, int to)
{
    int c = A->cols;
    for (int i = 0; i < c; i++)
    {
        A->elements[to][i] += A->elements[fr][i];
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

void scalar_multiply_matrix(Matrix A, elem_t n)
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
 * 
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
 * Reduces A to row-echelon form in place.
 * 
 * @param   A   the matrix to be reduced
 * @returns     the new determinant divided by the old. 
*/
elem_t row_echelon_form(Matrix A)
{
    register int r=A->rows, c=A->cols, i, j;
    elem_t result=1, factor;
    for (i = 0; i < c; i++)
    {
        for (j = i; j<r && nearly_zero(A->elements[j][i]); j++) ;  // j<r must be put in front
        if (j == r)
        {   // entire column is zero
            continue;
        }
        swap_row(A, i, j);
        result *= (i==j) ? 1 : -1;
        for (j = i+1; j < r; j++)
        {   // clear the rest of the column
            if (!nearly_zero(A->elements[j][i]))
            {   // if it needs clearing
                factor = - (A->elements[j][i]) / A->elements[i][i];
                multiply_row(A, i, factor);
                combine_row(A, i, j);
                result *= factor;
            }
        }
    }
    return result;
}

elem_t determinant(Matrix A)
{
    check_square(A);
    int r = A->rows;
    elem_t determinant = 1;
    Matrix Acpy = copy_matrix(A);

    if (r == 2)
    {   // ad - bc
        determinant = A->elements[0][0] * A->elements[1][1] - A->elements[0][1] * A->elements[1][0];
        return determinant;
    }
    
    determinant *= row_echelon_form(Acpy);
    for (int i = 0; i < r; i++)
    {
        determinant *= Acpy->elements[i][i];
    }
    return !nearly_zero(determinant) * determinant;
}

Matrix solve_linear_system(const Matrix A, const Vector b)
{
    register int d=A->cols, i, j;
    elem_t factor, pivot;
    Matrix Aug = append(A, vector_to_column_matrix(b), RIGHT);

    for (i = 0; i < d; i++)
    {
        for (j = i; j<d && nearly_zero(Aug->elements[j][i]); j++) ;
        if (j == d)
        {
            return NULL;
        }
        (i == j) ? : swap_row(Aug, i, j);
        for (j = i+1; j < d; j++)
        {
            if (!nearly_zero(Aug->elements[j][i]))
            {
                factor = - (Aug->elements[j][i]) / Aug->elements[i][i];
                multiply_row(Aug, i, factor);
                combine_row(Aug, i, j);
            }
        }
    }

    while (d-- > 1)
    {
        pivot = Aug->elements[d][d];
        for (int i = 0; i < d; i++)
        {
            if (!nearly_zero(Aug->elements[i][d]))
            {
                factor = - (Aug->elements[i][d]) / pivot;
                multiply_row(Aug, d, factor);
                combine_row(Aug, d, i);
            }
        }
    }

    return get_column_matrix(Aug, A->cols);
}

Matrix inverse(Matrix A)
{
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

Matrix append(const Matrix A, const Matrix B, const char mode)
{
    int Acols=A->cols, Arows=A->rows, Bcols=B->cols, Brows=B->rows;
    int Mrows, Mcols;
    Matrix M;
    switch(mode)
    {
        case UP:
        case DOWN:
        case LEFT:
        case RIGHT:
            Mrows = Arows;
            Mcols = Acols + Bcols;
            M = new_matrix(Mrows, Mcols);
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
            break;
        default:
            fprintf(stderr, "Unrecognized mode.\n");
            M = NULL;
            break;
    }

    return M;
}

Matrix basis_transform(Matrix A, Matrix T)
{
    return multiply_matrix(inverse(T), multiply_matrix(A, T));
}

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

Vector standard_unit_vector(int dim, int n)
{
    Vector b = new_vector(dim);
    b->elements[n] = 1;
    return b;
}

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