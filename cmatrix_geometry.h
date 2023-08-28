#ifndef _CMATRIX_GEOMETRY_H_
#define _CMATRIX_GEOMETRY_H_

#include "cmatrix_defs.h"

Matrix vector_projection_matrix(Vector);
Matrix rotation_matrix_2d(elem_t);
Matrix reflect_about_vector(Vector);
Matrix rotation_matrix_3d(elem_t, Vector);
Matrix reflect_about_plane(Vector, Vector);
Matrix reflect_about_normal(Vector);
Vector unit_normal(Vector, Vector);
Vector polar_to_cartesian(elem_t, elem_t);
Vector spherical_to_cartesian(elem_t, elem_t, elem_t);
Vector cylindrical_to_cartesian(elem_t, elem_t, elem_t);
Matrix attitude_description(char *, elem_t, elem_t, elem_t);

#endif