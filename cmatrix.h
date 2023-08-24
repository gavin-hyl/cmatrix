/**
 * CMatrix - a lightweight linear algebra library for embedded systems.
 * 
 * Includes functions and types for creating and manipulating matrices, as well as various
 * linear algebra operations. Also includes a set of generators for common
 * matrices.
 * 
 * Author: Gavin Hua (ghua@caltech.edu)
 * 
 * On function argument modification:
 * 
 * 1. If a function requires two or more matrix arguments, they will be const.
 * 2. If the function requires one matrix argument but returns either a matrix 
 * of different size or doesn't return a matrix, the arguments will be const.
 * 3. In all other cases, the arguments may be modified.
 * 
 * 
 * On function names:
 * 
 * 1. Any function returning a matrix has a type "M_".
 * 2. Any other function has no one-letter prefix. 
 * 
*/

#include "cmat_defs.h"
#include "cmat_basics.h"
#include "cmat_stdmat.h"
#include "cmat_vec.h"