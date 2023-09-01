/**
 * @file cmatrix_defs.h
 * @author Gavin Hua (139950129+GavinHYL@users.noreply.github.com)
 * @brief Type definitions, common macros, and configs for the cmatrix library.
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef _CMATRIX_DEFS_H_
#define _CMATRIX_DEFS_H_

#include <stdio.h>
#include <stdlib.h>

#define ELEMENT_PRINT_WIDTH 7       // width of elements in print_matrix and print_vector
#define ELEMENT_PRINT_PRECISION 2   // precision of elements in print_matrix and print_vector
#define EPSILON 1e-5                // the upper limit of what two floating points are considered equal
#define QR_ITER 20                  // maximum iterations of the QR algorithm.
#define QR_CHECK_PERIOD 10          // determines how often the matrix is checked to be upper-triangular in the QR algorithm.


// The floating-point type that is used in calculations. 
// Must either be double or float, NOT LONG DOUBLE, for maximum portability.
typedef double flt_t;  
#define nearly_zero(e) ((e)<=EPSILON && (e)>=EPSILON*-1)

typedef struct MaTrix {
    int rows;
    int cols;
    flt_t **elements;
} *Matrix;

enum {ROW_MAJOR, COLUMN_MAJOR};

typedef struct VecTor {
    int dim;
    flt_t *elements;
} *Vector;

/**
 * A series of error-checking macros are defined here as well.
 * They will exit the program immediately if the condition is not satisfied.
*/
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
    check_null(A);  \
    check_null(B);  \
    if (A->rows!=B->rows || A->cols!=B->cols)   \
    {   \
        fprintf(stderr, "Cannot operate on two matrices of different size.\n"); \
    }   \
})

#define check_equal_dimension(v, w) \
({  \
    check_null(v);  \
    check_null(w);  \
    if (v->dim != w->dim)   \
    {   \
        fprintf(stderr, "Cannot operate on two vectors of different dimension.\n"); \
        exit(1);    \
    }   \
})


#define check_multipliable(A, B)    \
({  \
    check_null(A);  \
    check_null(B);  \
    if (A->cols != B->rows) \
    {   \
        fprintf(stderr, "Cannot multiply when r1 is not equal to c2.\n");   \
        exit(1);    \
    }   \
})

#define check_square(A) \
({  \
    check_null(A);  \
    if (A->rows != A->cols) \
    {   \
        fprintf(stderr, "Cannot operate on non-square matrix.\n");  \
        exit(1);    \
    }   \
})

#endif