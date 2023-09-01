/**
 * @file cmatrix_geometry.h
 * @author Gavin Hua (139950129+GavinHYL@users.noreply.github.com)
 * @brief More advanced geometric operations not included in cmatrix_basics.h
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef _CMATRIX_GEOMETRY_H_
#define _CMATRIX_GEOMETRY_H_

#include "cmatrix_defs.h"

Matrix rotation_matrix_2d(flt_t);
Matrix rotation_matrix_3d(Vector, flt_t);
Matrix reflect_about_plane(Vector, Vector);
Matrix reflect_about_normal(Vector);
Vector unit_normal(Vector, Vector);
Vector polar_to_cartesian(flt_t, flt_t);
Vector spherical_to_cartesian(flt_t, flt_t, flt_t);
Vector cylindrical_to_cartesian(flt_t, flt_t, flt_t);
Matrix attitude_description(char *, flt_t *);

#endif