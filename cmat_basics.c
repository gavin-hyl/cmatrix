#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "cmat_basics.h"
#include "cmat_stdmat.h"

void m_checknull(Matrix A)
{
    if (A == NULL)
    {
        perror("ERROR: cannot manipulate null matrix.\n");
        exit(EXIT_FAILURE);
    }
}

Matrix M_new(int r, int c)
{
    Matrix A = (Matrix) malloc(sizeof(struct m));

    A->nrows = r;
    A->ncols = c;
    A->elements = (elem_t **) calloc(r, sizeof(elem_t *));
    for (int i = 0; i < r; i++)
    {
        A->elements[i] = (elem_t *) calloc(c, sizeof(elem_t));
    }
    return A;
}

Matrix M_copy(Matrix A)
{
    int r=A->nrows, c=A->ncols;
    Matrix cpy = M_new(r, c);

    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            cpy->elements[i][j] = A->elements[i][j];
        }
    }
    return cpy;
}

void set_by_value(Matrix A, elem_t val)
{
    m_checknull(A);
    int r=A->nrows, c=A->ncols;

    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            A->elements[i][j] = val;
        }
    }
}

void set_by_function(Matrix A, elem_t (*gen)(int, int))
{
    m_checknull(A);
    int r=A->nrows, c=A->ncols;
    
    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            A->elements[i][j] = (*gen)(i, j);
        }
    }
}

void set_by_line(Matrix A, elem_t *elements, char mode)
{
    m_checknull(A);
    int r=A->nrows, c=A->ncols;

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
            perror("[ERROR] Unknown mode (set_by_line): must be BY_ROW or BY_COLUMN.\n");
            break;
    }
}

void clear(Matrix A)
{
    m_checknull(A);
    set_by_value(A, 0);
}

void set_row(Matrix A, int r, elem_t *rowvals)
{
    int c = A->ncols;
    for (int i = 0; i < c; i++)
    {
        A->elements[r][i] = *rowvals++;
    }
}

Matrix M_get_row(Matrix A, int r)
{
    int c = A->ncols;
    Matrix row = M_new(1, c);

    for (int i = 0; i < c; i++)
    {
        row->elements[0][i] = A->elements[r][i];
    }
    return row;
}

void set_column(Matrix A, int c, elem_t *colvals)
{
    int r = A->nrows;
    for (int i = 0; i < r; i++)
    {
        A->elements[i][r] = *colvals++;
    }
}

Matrix M_get_column(Matrix A, int c)
{
    int r = A->nrows;
    Matrix col = M_new(r, 1);
    
    for (int i = 0; i < r; i++)
    {
        col->elements[i][0] = A->elements[i][c];
    }
    return col;
}

elem_t *flatten(Matrix A, char mode)
{
    int r=A->nrows, c=A->ncols;
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


void print_matrix(Matrix A, int width)
{
    m_checknull(A);
    int wid = (width == 0) ? ELEMENTWIDTH : width;
    int r=A->nrows, c=A->ncols;
    if (r == 1)
    {
        printf("( ");
        print_row(A, 0, wid);
        printf(" )");
        return;
    }

    for (int i = 0; i < r; i++)
    {
        if (i == 0)
        {
            printf("/ ");
            print_row(A, i, wid);
            printf(" \\");
        }
        else if (i == r-1)
        {
            printf("\\ ");
            print_row(A, i, wid);
            printf(" /");
        }
        else
        {
            printf("| ");
            print_row(A, i, wid);
            printf(" |");
        }
        printf("\n");
    }
    printf("\n");
}

void print_row(Matrix A, int r, int width)
{
    m_checknull(A);
    int c = A->ncols;

    for (int i = 0; i < c; i++)
    {
        printf("%*.*f", (width == 0) ? ELEMENTWIDTH : width, ELEMENTPREC, (float) A->elements[r][i]);
        printf((i == c-1) ? "" : " ");
    }
}

void swap_row(Matrix A, int r1, int r2)
{
    if (r1 != r2)
    {
        elem_t *temp = M_get_row(A, r1)->elements[0];
        set_row(A, r1, M_get_row(A, r2)->elements[0]);
        set_row(A, r2, temp);
    }
}

void multiply_row(Matrix A, int r, elem_t k)
{
    if (!nearlyzero(k-1))
    {
        Matrix row = M_get_row(A, r);
        scalar_multiply(row, k);
        set_row(A, r, row->elements[0]);
    }
}

void combine_row(Matrix A, int fr, int to)
{
    int c = A->ncols;
    for (int i = 0; i < c; i++)
    {
        A->elements[to][i] += A->elements[fr][i];
    }
}

Matrix M_add(Matrix A, Matrix B)
{
    m_checknull(A);
    m_checknull(B);
    int r1=A->nrows, c1=A->ncols, r2=B->nrows, c2=B->ncols;
    
    if (r1!=r2 || c1!=c2)
    {
        printf("ERROR (A_add): Arices must be same size.\n");
    }
    Matrix result = M_new(r1, c1);
    for (int i = 0; i < r1; i++)
    {
        for (int j = 0; j < c1; j++)
        {
            result->elements[i][j] = A->elements[i][j] + B->elements[i][j];
        }
    }
    return result;
}

Matrix M_multiply(Matrix A, Matrix B)
{
    m_checknull(A);
    m_checknull(B);
    int r1=A->nrows, c1=A->ncols, r2=B->nrows, c2=B->ncols;

    if (c1 != r2)
    {
        printf("ERROR (A_multiply): nrow of M1 must be equal to ncol of M2.\n");
        return NULL;
    }

    Matrix result = M_new(r1, c2);
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

void scalar_multiply(Matrix A, elem_t n)
{
    m_checknull(A);
    int r=A->nrows, c=A->ncols;

    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            A->elements[i][j] *= n;
        }
    }
}

Matrix M_transpose(const Matrix A)
{
    int r=A->nrows, c=A->ncols;
    Matrix result = M_new(c, r);

    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < c; j++)
        {
            result->elements[j][i] = A->elements[i][j];
        }
    }
    return result;
}

elem_t row_echelon_form(Matrix A)
{
    register int r=A->nrows, c=A->ncols, i, j;
    elem_t result=1, factor;
    for (i = 0; i < c; i++)
    {
        for (j = i; j<r && nearlyzero(A->elements[j][i]); j++) ;  // j<r must be put in front
        if (j == r)
        {   // entire column is zero
            continue;
        }
        swap_row(A, i, j);
        result *= (i==j) ? 1 : -1;
        for (j = i+1; j < r; j++)
        {   // clear the rest of the column
            if (!nearlyzero(A->elements[j][i]))
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
    if (!issquare(A))
    {
        return 0;
    }
    int r = A->nrows;
    elem_t determinant = 1;
    Matrix Acpy = M_copy(A);

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
    return !nearlyzero(determinant) * determinant;
}

/**
 * Solves an equation in the form Ax = b, where A is square and non-singular
 * @param A
 * @param b
 * @returns the solution vector, or NULL if A is singular.
*/
Matrix M_solve(const Matrix A, const Matrix b)
{
    register int d=A->ncols, i, j;
    elem_t factor, pivot;
    Matrix Aug = M_append(A, b, RIGHT);

    for (i = 0; i < d; i++)
    {   // reduce to row-echelon form
        for (j = i; j<d && nearlyzero(Aug->elements[j][i]); j++) ;  // j<r must be put in front
        if (j == d)
        {   // entire column is zero
            return NULL;
        }
        (i == j) ? : swap_row(Aug, i, j);
        for (j = i+1; j < d; j++)
        {   // clear the rest of the column
            if (!nearlyzero(Aug->elements[j][i]))
            {   // if it needs clearing
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
            if (!nearlyzero(Aug->elements[i][d]))
            {
                factor = - (Aug->elements[i][d]) / pivot;
                multiply_row(Aug, d, factor);
                combine_row(Aug, d, i);
            }
        }
    }

    return M_get_column(Aug, A->ncols);
}

/**
 * ? why does this work
*/
Matrix M_inverse(Matrix A)
{
    int d = A->nrows;
    Matrix Acpy = M_copy(A);
    Matrix I = M_identity(d);

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

Matrix M_append(const Matrix A, const Matrix B, const char mode)
{
    int Acols=A->ncols, Arows=A->nrows, Bcols=B->ncols, Brows=B->nrows;
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
            M = M_new(Mrows, Mcols);
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
            perror("unrecognized mode.");
            M = NULL;
            break;
    }

    return M;
}

Matrix M_basis_transform(Matrix A, Matrix T)
{
    return M_multiply(M_inverse(T), M_multiply(A, T));
}
