#ifndef _CMAT_VEC_H_
#define _CMAT_VEC_H_

#include "cmat_defs.h"

/**
 * This 
*/

int is_column(Matrix);
int is_row(Matrix);
void normalize(Matrix);
elem_t length(Matrix);
elem_t dot_product(Matrix, Matrix);
Matrix cross_product(Matrix, Matrix);

#endif