#ifndef _CMATRIX_DEFS_H_
#define _CMATRIX_DEFS_H_

#define EPSILON 0.01
#define ELEMENTWIDTH 7  // for printing the matrix, has no influence on calculations.
#define ELEMENTPREC 2

typedef double elem_t;  // must be either float or double, NOT long double
#define nearly_zero(e) (e<EPSILON && e>(-1)*EPSILON)

typedef struct MaTrix {
    int rows;
    int cols;
    elem_t **elements;
} *Matrix;

enum {ROW_MAJOR, COLUMN_MAJOR};

typedef struct VecTor {
    int dim;
    elem_t *elements;
} *Vector;

// const struct Vector e2x = E2X;


// A series of error-checking macros are defined here as well.

#include <stdio.h>
#include <stdlib.h>

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