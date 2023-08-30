#include <math.h>
#include "cmatrix_algebra.h"
#include "cmatrix_basics.h"
#include "cmatrix_geometry.h"

/**
 * @brief Swaps two rows of a matrix in place.
 * 
 * @param A the matrix
 * @param r1 one row to be swapped
 * @param r2 the other row to be swapped
 */
void swap_row(Matrix A, int r1, int r2)
{
    if (r1 != r2)
    {
        int c = A->cols;
        elem_t temp;
        for (int i = 0; i < c; i++)
        {
            temp = A->elements[r1][i];
            A->elements[r1][i] = A->elements[r2][i];
            A->elements[r2][i] = temp;
        }
    }
}

/**
 * @brief Multiplies a row of a matrix by a scalar constant in place.
 * 
 * @param A the matrix
 * @param r the row (index from 0)
 * @param k the scaling constant
 */
void multiply_row(Matrix A, int r, elem_t k)
{
    if (!nearly_zero(k-1))
    {
        int c = A->cols;
        for (int i = 0; i < c; i++)
        {
            A->elements[r][i] *= k;
        }
    }
}

/**
 * @brief Combines two rows of a matrix in place.
 * 
 * @param A the matrix
 * @param fr the row to add
 * @param to the row to which fr is added
 */
void combine_row(Matrix A, int fr, int to)
{
    if (fr != to)
    {
        int c = A->cols;
        for (int i = 0; i < c; i++)
        {
            A->elements[to][i] += A->elements[fr][i];
        }
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
 * @brief Reduces A to row-echelon form in place using Gaussian elimination.
 * 
 * @param A the matrix to be reduced
 * @return elem_t: the new determinant divided by the old. 
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
 * @brief Computes the determinant of a square matrix using Gaussian elimination 
 * and multipling the elements on the diagonal of the row-echelon form matrix.
 * 
 * @param A the matrix of which determinant is to be calculated
 * @return elem_t: the determinant
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
 * @brief Solves a linear system of the form Ax = b using Gaussian elimination. 
 * 
 * @param A the coefficient matrix
 * @param b the constants
 * @return Vector: the solution vector
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
 * @param A the matrix in E representation (E is any representation, does
 * not need to be specified).
 * @param T the representation of the basis vectors of E in the desired basis.
 * @return Matrix: the matrix after basis transformation
 */
Matrix basis_transform(const Matrix A, const Matrix T)
{
    return multiply_matrix(inverse(T), multiply_matrix(A, T));
}

/**
 * @brief Tests for linear independence of a set of vectors with Gaussian elimination.
 * 
 * @param vecs the vectors to be tested
 * @param nvecs the number of vectors to be tested
 * @return int: 1 if they are linearly independent, 0 if not
 */
int linear_independence(Vector *vecs, int nvecs)
{
    register int d = vecs[0]->dim;
    if(nvecs > d)
    {
        return 0;
    }
    Matrix vec_matrix = new_matrix(nvecs, d);
    for (int i = 0; i < nvecs; i++)
    {
        set_matrix_row(vec_matrix, i, vecs[i]->elements);
    }
    row_echelon_form(vec_matrix);
    return !vector_is_equal(new_vector(d), get_row_vector(vec_matrix, nvecs-1));
}

int row_linear_independence(Matrix vecs)
{
    row_echelon_form(vecs);
    return !vector_is_equal(new_vector(vecs->cols), get_row_vector(vecs, vecs->rows-1));
}


/**
 * @brief Computes an orthognal basis for a set of vectors using the Gram-Schmidt process.
 * 
 * Assumes the set to be linearly independent.
 * 
 * @param vecs the set of vectors
 * @param nvecs the number of vectors in the set
 * @return Matrix: the vectors consisting the orthogonal basis
 */
Matrix gram_schmidt_orthogonal_basis(Vector *vecs, int nvecs)
{
    int d = vecs[0]->dim;
    Matrix basis = vector_to_column_matrix(vecs[0]);
    Matrix vec_to_perp = identity(d);

    for (int i = 1; i < nvecs; i++)
    {
        Vector vi_perp = vecs[i];
        Matrix v_prev_proj = vector_projection_matrix(vecs[i-1]);
        scale_matrix(v_prev_proj, -1);
        vec_to_perp = add_matrix(vec_to_perp, v_prev_proj);
        vi_perp = multiply_matrix_vector(vec_to_perp, vi_perp);
        basis = append_horizontal(basis, vector_to_column_matrix(vi_perp));
    }

    return basis;
}

/**
 * @brief Computes the QR decomposition for a matrix.
 * 
 * @param A the matrix to be decomposed
 * @return Matrix*: [0] is Q, and [1] is R.
 */
Matrix *qr_decomposition(Matrix A)
{

}

Matrix eigenvectors(Matrix A)
{

}

elem_t *eigenvalues(Matrix A)
{

}