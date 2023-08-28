#include <math.h>
#include "cmatrix_algebra.h"
#include "cmatrix_basics.h"

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
        scale_matrix(row, k);
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

// UNUSED - probably useless
// Matrix row_combine_matrix(int d, int fr, int to)
// {
//     Matrix combrow = identity(d);
//     combrow->elements[to][fr] = 1;
//     return combrow;
// }

// Matrix row_swap_matrix(int d, int r1, int r2)
// {
//     Matrix swaprow = identity(d);
//     swaprow->elements[r1][r1] = 0;
//     swaprow->elements[r2][r2] = 0;
//     swaprow->elements[r1][r2] = 1;
//     swaprow->elements[r2][r1] = 1;
//     return swaprow;
// }

// Matrix row_multiply_matrix(int d, int r, elem_t mult)
// {
//     Matrix multrow = identity(d);
//     multrow->elements[r][r] = mult;
//     return multrow;
// }

/**
 * @brief Reduces A to row-echelon form in place.
 * 
 * @param A the matrix to be reduced
 * @return the new determinant divided by the old. 
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

/**
 * @brief Computes the determinant of a square matrix.
 * 
 * We perform Gaussian elimination, and multiply the elements on the diagonal.
 * 
 * @param A the matrix of which determinant is to be calculated
 * @return the determinant
 */
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

/**
 * @brief Solves a linear system of the form Ax = b.
 * 
 * The algorithm used is Gaussian elimination. 
 * 
 * @param A the coefficient matrix
 * @param b the constants
 * @return the solution vector
 */
Vector solve_linear_system(const Matrix A, const Vector b)
{
    register int d=A->cols, i, j;
    elem_t factor, pivot;
    Matrix Aug = append_horizontal(A, vector_to_column_matrix(b));

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

    return get_column_vector(Aug, A->cols);
}

/**
 * @brief Transforms a matrix between basis representations.
 * 
 * @param A the matrix in E representation (E is any representation, that does
 * not need to be specified).
 * @param T the representation of the basis vectors of E in the desired basis.
 * @return the matrix after basis transformation 
 */
Matrix basis_transform(const Matrix A, const Matrix T)
{
    return multiply_matrix(inverse(T), multiply_matrix(A, T));
}

Matrix *qr_decomposition(Matrix A)
{

}

Matrix eigenvectors(Matrix A)
{

}

elem_t *eigenvalues(Matrix A)
{

}