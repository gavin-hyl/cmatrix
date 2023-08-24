#ifndef _CMAT_DEFS_H_
#define _CMAT_DEFS_H_

#define EPSILON 0.01
#define nearlyzero(e) ((e)<EPSILON && (e)>-EPSILON)
#define issquare(mat) ((mat)->nrows == (mat)->ncols)
#define ELEMENTWIDTH 7
#define ELEMENTPREC 2

typedef double elem_t;

typedef struct m {
    int nrows;
    int ncols;
    elem_t **elements;
} *Matrix;

#endif