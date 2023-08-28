#include <math.h>
#include "cmatrix_basics.h"
#include "cmatrix_geometry.h"

/**
 * @brief Computes the projection matrix of the space onto a vector. 
 * 
 * Formula is (v * vT) / ||v||^2.
 * 
 * @param v the vector to be projected on to.
 * @return the projection matrix. 
 */
Matrix vector_projection_matrix(const Vector v)
{
    check_null(v);
    elem_t length = vector_length(v);
    if (length == 0)
    {
        return identity(v->dim);
    }
    Matrix proj_v = multiply_matrix(vector_to_column_matrix(v), vector_to_row_matrix(v));
    scale_matrix(proj_v, 1 / dot_product(v, v));
    return proj_v;
}

/**
 * @brief A 2D counterclockwise rotation centered at the origin.
 * 
 * @param angle the angle to rotate 
 * @return the rotation matrix
*/
Matrix rotation_matrix_2d(const elem_t angle)
{
    Matrix mat = new_matrix(2, 2);
    
    mat->elements[0][0] = cos(angle);
    mat->elements[1][0] = sin(angle);
    mat->elements[0][1] = -sin(angle);
    mat->elements[1][1] = cos(angle);
    return mat;
}

/**
 * @brief Computes the reflection matrix of the n-D space about a vector.
 * 
 * The formula used is (2P-I), where P is the projection matrix and I is the identity.
 * 
 * @param v the vector about which the space is reflected
 * @return the reflection matrix. 
 */
Matrix reflect_about_vector(const Vector v)
{
    check_null(v);
    int d = v->dim;

    Matrix mir = vector_projection_matrix(v);
    scale_matrix(mir, 2);
    Matrix I = identity(d);
    scale_matrix(I, -1);
    return add_matrix(mir, I);
}

/**
 * @brief Computes the rotation matrix of the 3D space around a vector.
 * 
 * @param angle the rotation angle
 * @param axis the axis around which the space is rotated
 * @return the rotation matrix.
 */
Matrix rotation_matrix_3d(elem_t angle, Vector axis)
{

}

/**
 * @brief Computes the reflection matrix of the 3D space about a plane.
 * 
 * We find the normal vector of the plane and create its projection matrix.
 * Then, we employ an algorithm similar to vector reflection, subtracting the 
 * projection twice from the identity. 
 * 
 * @param v1 the first vector that defines the plane.
 * @param v2 the second vector that defines the plane.
 * @return the reflection matrix.
 */
Matrix reflect_about_plane(Vector v1, Vector v2)
{
    check_equal_dimension(v1, v2);
    Vector normal = cross_product(v1, v2);
    Matrix normal_project = vector_projection_matrix(normal);
    scale_matrix(normal_project, -2);
    return add_matrix(normal_project, identity(v1->dim));
}

Matrix reflect_about_normal(Vector v)
{   
    Matrix normal_project = vector_projection_matrix(v);
    scale_matrix(normal_project, -2);
    return add_matrix(normal_project, identity(v->dim));
}

Vector unit_normal(Vector v, Vector w)
{
    Vector unit_norm = cross_product(v, w);
    scale_vector(unit_norm, 1 / (vector_length(v) * vector_length(w)));
    return unit_norm;
}

Vector polar_to_cartesian(elem_t r, elem_t phi)
{

}

Vector spherical_to_cartesian(elem_t r, elem_t phi, elem_t theta)
{

}

Vector cylindrical_to_cartesian(elem_t r, elem_t phi, elem_t z)
{

}

Matrix attitude_description(char *mode, elem_t a1, elem_t a2, elem_t a3)
{

}