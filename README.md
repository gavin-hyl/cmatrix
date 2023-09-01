# cmatrix - C-based Linear Algebra
A lightweight C-based linear algebra library. Work in progress, please fork/raise an issue if you'd like to contribute!

# File Structure (How to Use)
## cmatrix.h 
Includes everything that the library has to offer, essentially a one-stop approach for using this library. 
## cmatrix_defs.h 
Includes type definitions for vectors/matrices and basic error checking macros.

## cmatrix_basics.h 
Includes all the basic opertions needed for matrix/vector creation and manipulation.

## cmatrix_algebra.h
Includes more advanced algebraic opertions, such as constructing Gram-Schmidt orthogonal bases and performing QR decompositions.

## cmatrix_geometry.h 
Includes generators for various rotation, reflection, and projection matrices in 2D and 3D. Also includes functions for coordinate system transformations.
