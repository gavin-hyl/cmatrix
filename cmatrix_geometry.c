/**
 * @file cmatrix_geometry.c
 * @author Gavin Hua (139950129+GavinHYL@users.noreply.github.com)
 * @brief Implementation of cmatrix_geometry.h
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#include <math.h>
#include <string.h>
#include "cmatrix_basics.h"
#include "cmatrix_geometry.h"

/**
 * @brief A 2D counterclockwise rotation centered at the origin.
 * 
 * @param angle the angle to rotate 
 * @return the rotation matrix
*/
Matrix rotation_matrix_2d(const flt_t angle)
{
    Matrix mat = new_matrix(2, 2);
    
    mat->elements[0][0] = cos(angle);
    mat->elements[1][0] = sin(angle);
    mat->elements[0][1] = -sin(angle);
    mat->elements[1][1] = cos(angle);
    return mat;
}

/**
 * @brief Computes the rotation matrix of the 3D space around a vector.
 * 
 * @param axis the axis around which the space is rotated
 * @param theta the rotation angle
 * @return the rotation matrix
 */
Matrix rotation_matrix_3d(const Vector axis, const flt_t theta)
{
    normalize(axis);
    flt_t ux=axis->elements[0], uy=axis->elements[1], uz=axis->elements[2];
    flt_t c=cos(theta), s=sin(theta);
    Matrix rot_3d = new_matrix(3, 3);
    set_matrix_row(rot_3d, 0, (flt_t []) {c + ux*ux*(1-c), ux*uy*(1-c) - uz*s, ux*uz*(1-c) + uy*s});
    set_matrix_row(rot_3d, 1, (flt_t []) {ux*uy*(1-c) + uz*s, c + uy*uy*(1-c), uy*uz*(1-c) - ux*s});
    set_matrix_row(rot_3d, 2, (flt_t []) {ux*uz*(1-c) - uz*s, uy*uz*(1-c) + ux*s, c + uz*uz*(1-c)});
    return rot_3d;
}

/**
 * @brief Computes the reflection matrix of the 3D space about a plane.
 * 
 * @param v1 the first vector to specify the plane
 * @param v2 the second vector to specify the plane
 * @return the reflection matrix
 */
Matrix reflect_about_plane(const Vector v1, const Vector v2)
{
    check_equal_dimension(v1, v2);
    Vector normal = cross_product(v1, v2);
    Matrix normal_project = vector_projection_matrix(normal);
    normal_project = scale_matrix(normal_project, -2);
    return add_matrix(normal_project, identity(v1->dim));
}

/**
 * @brief Computes the reflection matrix of the 3D space about a plane's normal vector.
 * 
 * @param v the normal vector
 * @return Matrix: the reflection matrix
 */
Matrix reflect_about_normal(const Vector v)
{   
    Matrix normal_project = vector_projection_matrix(v);
    normal_project = scale_matrix(normal_project, -2);
    return add_matrix(normal_project, identity(v->dim));
}

/**
 * @brief Computes the unit normal vector of a plane.
 * 
 * @param v the first vector to specify a plane
 * @param w the second vector to specify a plane
 * @return the unit normal vector of the plane
 */
Vector unit_normal(const Vector v, const Vector w)
{
    Vector unit_norm = cross_product(v, w);
    unit_norm = scale_vector(unit_norm, 1 / (norm(v) * norm(w)));
    return unit_norm;
}

/**
 * @brief Converts from polar coordinates to cartesian coordinates in 2D.
 * 
 * @param r the length of the vector
 * @param phi the angle between the vector and the x-axis
 * @return the resultant vector in the standard orthonormal basis
 */
Vector polar_to_cartesian(const flt_t r, const flt_t phi)
{
    Vector v = new_vector(2);
    v->elements[0] = r * cos(phi);
    v->elements[1] = r * sin(phi);
    return v;
}

/**
 * @brief Converts from spherical coordinates to cartesian coordinates in 3D.
 * 
 * @param r the length of the vector
 * @param phi the angle on the xOy plane (as in polar coordinates) 
 * @param theta the angle between the vector and the z-axis
 * @return the resultant vector in the standard orthonormal basis
 */
Vector spherical_to_cartesian(const flt_t r, const flt_t phi, const flt_t theta)
{
    Vector v = new_vector(3);
    v->elements[0] = r * sin(theta) * cos(phi);
    v->elements[1] = r * sin(theta) * sin(phi);
    v->elements[2] = r * cos(theta);
    return v;
}

/**
 * @brief Converts from cylindrical coordinates to cartesian coordinates in 3D.
 * 
 * @param r the length of the vector
 * @param phi the angle on the xOy plane (as in polar coordinates)
 * @param z the z-component of the point
 * @return the resultant vector in the standard orthonormal basis
 */
Vector cylindrical_to_cartesian(const flt_t r, const flt_t phi, const flt_t z)
{
    Vector v = new_vector(3);
    v->elements[0] = r * cos(phi);
    v->elements[1] = r * sin(phi);
    v->elements[2] = z;
    return v;
}

/**
 * @brief Computes the rotation matrix corresponding to an attitude description.
 * 
 * @param mode a string specifying the the order of rotations, eg. "313" for 
 * rotation in the order z-x-z
 * @param angles the rotation angles, must be in the same order specified by mode
 * @return the rotation matrix
 */
Matrix attitude_description(char *mode, flt_t *angles)
{
    int axis;
    Matrix attitude = identity(3);
    while (axis = *mode++)
    {
        switch(axis)
        {
            case '1':
                attitude = multiply_matrix(rotation_matrix_3d(get_column_vector(attitude, 0), *angles++), attitude);
                break;
            case '2':
                attitude = multiply_matrix(rotation_matrix_3d(get_column_vector(attitude, 1), *angles++), attitude);
                break;
            case '3':
                attitude = multiply_matrix(rotation_matrix_3d(get_column_vector(attitude, 2), *angles++), attitude);
                break;
            default:
                fprintf(stderr, "Unknown axis in attitude description.\n");
                exit(1);
        }
    }
    return attitude;
}