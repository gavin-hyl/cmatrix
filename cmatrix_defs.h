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

#define ELEMENT_PRINT_WIDTH 7
#define ELEMENT_PRINT_PRECISION 2
#define EPSILON 1e-5
#define QR_ITER 20
#define QR_CHECK_PERIOD 10


// The floating-point type that is used in calculations. 
// Must either be double or float, NOT LONG DOUBLE, for maximum portability.
typedef double flt_t;  
#define nearly_zero(e) ((e)==0)

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

#define ex2 std_unit_vector(2, 0)
#define ey2 std_unit_vector(2, 1)
#define ex3 std_unit_vector(3, 0)
#define ey3 std_unit_vector(3, 1)
#define ez3 std_unit_vector(3, 2)

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