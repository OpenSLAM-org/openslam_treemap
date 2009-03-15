/*
Copyright (c) 2009, Universitaet Bremen
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Universitaet Bremen nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

//!\author Udo Frese

/*!\file xymOperations.h Contains functions providing operations on vectors and
   matrices. (i.e. Add, Multiply, Cholesky). The function name have global scope
   but all use an xy matrix package class type in its signature.
*/

#ifndef XYMOPERATIONS_H
#define XYMOPERATIONS_H

#include "xymMatrixVC.h"
#include "xymMatrixC.h"
#include "xymVectorV.h"
#include "xymVector.h"

#include <vector>

//************************************************************
/*!\name Vector operations
 */
//@{


//! Returns the scalar product between 'v1' and 'v2'
double dot (const XymVectorV& v1, const XymVectorV& v2);

//! Returns the 2 norm of v
double norm (const XymVectorV& v);

//! Returns the infinity norm of v
double norminf (const XymVectorV& v);

//! Returns the E norm of A (root of sum of square entries)
double normE (const XymMatrixVC& A);


//! Sets 'res=v1+v2', operands can be identical
void add (XymVectorV& res, const XymVectorV& v1, const XymVectorV& v2);

//! Sets 'res=v1+v2', operands can be identical. 
/*! The result is created if necessary */
void add (XymVector& res, const XymVectorV& v1, const XymVectorV& v2);

//! Sets 'res=v1-v2', operands can be indentical
void sub (XymVectorV& res, const XymVectorV& v1, const XymVectorV& v2);

//! Sets 'res=v1-v2', operands can be indentical
/*! The result is created if necessary. */
void sub (XymVector& res, const XymVectorV& v1, const XymVectorV& v2);

//! Sets \c res=u*factor. Both may be identical
void scale (XymVectorV& res, const XymVectorV& u, double scaleFactor);

//! Sets \c res=u*factor. Both may be identical
/*! The result is created if necessary */
void scale (XymVector& res, const XymVectorV& u, double scaleFactor);

//! Sets \c res+=u*factor.
void scaleAdd (XymVectorV& res, const XymVectorV& u, double scale);

//! Computes \c res = lambda*u + mu*v. 
/*!
  \c res, \c u and \c v can be the same.
*/
void linear (XymVectorV& res, const XymVectorV& u, double lambda, const XymVectorV& v, double mu);

//! Computes \c res = lambda*u + mu*v. 
/*!
  \c res, \c u and \c v can be the same. \c res can be uninitialised.
*/
void linear (XymVector& res, const XymVectorV& u, double lambda, const XymVectorV& v, double mu);

//@}

//************************************************************
/*!\name Matrix operations
 */
//@{

//! Sets \c res=a*b*factor. 
/*! 'res' must either have the correct dimension or be uninitialised
 */
void multiply (XymMatrixVC& res, const XymMatrixVC& a, const XymMatrixVC& b, double factor=1);

//! Sets \c res=res+a*b*factor.
void multiplyAdd (XymMatrixVC& res, const XymMatrixVC& a, const XymMatrixVC& b, double factor=1);

//! Sets c\ a=lambda*a. 
void multiply (XymMatrixVC& a, double lambda);

//! Sets c\ res=lambda*a. 
void multiply (XymMatrixVC& res, const XymMatrixVC& a, double lambda);

//! Sets \c res=a*v
void multiply (XymVectorV& res, const XymMatrixVC& a, const XymVectorV& v, double factor=1);

//! Sets \c res=a*v
void multiply (XymVector& res, const XymMatrixVC& a, const XymVectorV& v, double factor=1);

//! Sets \c res+=a*v
void multiplyAdd (XymVectorV& res, const XymMatrixVC& a, const XymVectorV& v, double factor=1);

//! Sets 'res = at*v'
void multiplyT (XymVectorV& res, const XymMatrixVC& a, const XymVectorV& v, double factor=1);

//! Sets \c res=at*v, \c res may be unallocated
void multiplyT (XymVector& res, const XymMatrixVC& a, const XymVectorV& v, double factor=1);

//! Sets 'res += at*v'
void multiplyTAdd (XymVectorV& res, const XymMatrixVC& a, const XymVectorV& v, double factor=1);


//! Sets res=res+a
void add (XymMatrixVC& res, const XymMatrixVC& a);

//! Sets res=a+b
void add (XymMatrixVC& res, const XymMatrixVC& a, const XymMatrixVC& b);

//! Sets res=a+b
/*! If \c res is empty it is created. */
void add (XymMatrixC& res, const XymMatrixVC& a, const XymMatrixVC& b);

//! Sets res=res-a
void sub (XymMatrixVC& res, const XymMatrixVC& a);

//! Sets res=a-b
void sub (XymMatrixVC& res, const XymMatrixVC& a, const XymMatrixVC& b);


//! Sets res=a-b
/*! If \c res is empty it is created. */
void sub (XymMatrixC& res, const XymMatrixVC& a, const XymMatrixVC& b);

//! Sets res+=lambda*a
void linear (XymMatrixVC& res, const XymMatrixVC& a, double lambda);

//! Sets \c res+=a*lambda. Alias for linear()
inline void multiplyAdd(XymMatrixVC &res, const XymMatrixVC& a, double lambda){
    linear(res, a, lambda);
}

//! Sets \c res=A*lambda. Both can be identical
void scale (XymMatrixVC& res, const XymMatrixVC& A, double lambda);

//! Computes \c M+=scale*vv^t, that is performs a rank one update of
//! a symmetric matrix.
void addvvt (XymMatrixVC& res, const XymVectorV& v, double scale=1.0);

//! Computes \c M+=scale*uv^t, that is performs a rank one update of
//! a matrix.
void adduvt (XymMatrixVC& res, const XymVectorV& u, const XymVectorV& v, double scale=1.0);

//! Convenience function that returns the outer product \f$ scale*uv^T \f$ of two vectors
/*! Note that calling \c adduvt is faster. */
XymMatrixC outer (const XymVectorV& u, const XymVectorV& v, double scale=1.0);

//! Computes \c M += scale*JCJ^t, that is performs a rank \c J.cols() update of
//! a symmetric matrix. \c C must be symmetric too.
void addJCJt (XymMatrixVC& res, const XymMatrixVC& J, const XymMatrixVC& C, double scale=1.0);

//! Computes \c u^TMu of a symmetric matrix \c v
/*! This is about a factor of two faster than using \c VMV.
 */
double symVMV (const XymVectorV& u, const XymMatrixVC& M);

//! Computes \c u^TMv of a symmetric matrix \c v
double VMV (const XymVectorV& u, const XymMatrixVC& M, const XymVectorV& v);

//! Internal function for a dot product
double xymDot (double* a, double* b, int n);

//! Internal  function for a dot product
double xymDot (double* a, int aOfs, double* b, int bOfs, int n);

//! trace of matrix A (sum of diagonal entries)
double trace (const XymMatrixVC& A);

//! Multiply two matrics.
/*! 'res' can be uninitialised.
*/
void multiply (XymMatrixC& res, const XymMatrixVC& a, const XymMatrixVC& b, double factor=1);

//! Transposes matrix A
void transpose (XymMatrixC& result, const XymMatrixVC& A);

//! Transposes matrix A
void transpose (XymMatrixVC& result, const XymMatrixVC& A);

//! Transposes square matrix A
void transpose (XymMatrixVC& A);

//! Replaces \c A by \c \f$ \frac{A+A^T}{2} \f$ to make it symmetric again
void symmetrize (XymMatrixVC& A);

//! Computes the maximal difference between 'A(i,j)' and 'A(J,i)'
/*! Can be used to check, whether a symmetric matrix is really symmetric.
 */
double symmetry (const XymMatrixVC& A);

//! Copies the lower half of \c A into the upper half
void copyLowerToUpper (XymMatrixVC& A);

//! Copies the upper half of \c A into the lower half
void copyUpperToLower (XymMatrixVC& A);

//! Sets the strict lower half of \c A to 0 
/*! Works fdr general rectangular matrices. */
void zeroLower (XymMatrixVC& A);

//! Sets the strict upper half of \c A to 0
/*! Works fdr general rectangular matrices. */
void zeroUpper (XymMatrixVC& A);



//! Adds in an efficient way \c factor*Q(i,j) to \c A(idxI[i], idxJ[j])
/*! Indices -1 are ommited. \c idxI and \c idxJ can be identical. */
void addIndexed (XymMatrixVC& A, const XymMatrixVC& Q, const vector<int>& idxI, const vector<int>& idxJ, double factor=1);

//! Adds in an efficient way \c factor*Q(i,j) to \c A(i, idxJ[j])
void addIndexedCol (XymMatrixVC& A, const XymMatrixVC& Q, const vector<int>& idxJ, double factor=1);

//! Adds in an efficient way v(i) to u(idx[i])
/*! Indices -1 are ommited. */
void addIndexed (XymVectorV& u, const XymVectorV& v, const vector<int>& idx, double facto=1);

//! Rotates a matrix representing alternately x and y coordinates by \c angle
/*! \c A must have even dimensions (not necessarily square). Let R be a block diagonal
    matrix where each block is \c ((cos angle, -sin angle), (sin angle, cos angle))
    Then \c result is \c R A R^T. Result and A may be identical.

    This form of rotation means: all eigenvalues rotate by the same angle. For instance,
    if a is one-dimensional A=uu^T, than Result = RAR^T = Ruu^TR^T = (Ru)(Ru)^T is one
    dimensional being generated by a rotated vector.

    This is equivalent to calling \c rotate2DRows followed by \c rotate2DCols.
 */
void rotate2D (XymMatrixVC& result, const XymMatrixVC& A, double angle);

//! See \c rotate2D 
/*! If \c result is uninitialised, it is initialised with the right format. */
void rotate2D (XymMatrixC& result, const XymMatrixVC& A, double angle);

//! Rotates the rowspace of a matrix alternately x and y coordinates by \c angle
/*! Compare \c rotate2D. This is equivalent to computing \c R A. */
void rotate2DRows (XymMatrixVC& result, const XymMatrixVC& A, double angle);

//! See \c rotate2DRows
/*! If \c result is uninitialised, it is initialised with the right format. */
void rotate2DRows (XymMatrixC& result, const XymMatrixVC& A, double angle);


//! Rotates the colspace of a matrix alternately x and y coordinates by \c angle
/*! Compare \c rotate2D. This is equivalent to computing \c A R^T. */
void rotate2DCols (XymMatrixVC& result, const XymMatrixVC& A, double angle);

//! See \c rotate2DCols
/*! If \c result is uninitialised, it is initialised with the right format. */
void rotate2DCols (XymMatrixC& result, const XymMatrixVC& A, double angle);




//! Transforms (translation, rotation) a vector representing alternately x and y coordinates
/*! \c v must have even dimension. Each block of two values (x, y) is
    rotated by angle and then \c (dx, dy) is added.
 */
void transform2D (XymVectorV& result, const XymVectorV& v, 
                  double dx, double dy, double angle);

//! See \c transform2D
/*! If \c result is uninitialised, it is initialised with the right format. */
void transform2D (XymVector& result, const XymVectorV& v, 
                  double dx, double dy, double angle);


//! Cholesky decomposition of a symmetric positive definite matrix A
/*! \c U is a upper triangular matrix and \f$ U^TU \f$ equals \c A.
    \c U is just the transpose of the more familiar matrix \c L of
    the \f$ LL^T=A \f$ decomposition. This way, memory access in the
    inner loop is more efficient.

    If the cholesky decomposition fails (matrix not SPD), than U will
    be uninitialised (\c U.isValid()==false). If \c minDiagonal>0, the
    routine will return a failure (\c U uninitialised) if at any stage
    of the algorithm the diagonal entry is \c <minDiagonal. However if
    \c setToMin is \c true, the diagonal entry is set to \c
    minDiagonal and the decoposition is continued. NOTE: This feature
    is a not very good way to deal with singular matrices. It is mch
    better to use SVD (LAPACK) instead.
  
    This code is only provided to allow very simple linear algebra to be
    more self contained. If you need more sophisticated operations
    use LAPACK. (This class is specifically designed to cooperate
    with LAPACK)

    \sa 
    Numerical Recipes Chapter 2.9
*/
void cholesky (XymMatrixC& U, const XymMatrixVC& A, double minDiagonal=0, bool setToMin=false);

//! Solve the equation \f$ Ax=b \f$, where \c U is the cholesky decomposition of \c A.
/*! \c A must be symmetric positive definite.

    If \c U is not valid (for instance passing \c XymMatrixVC()) \c cholesky is called.
    if \c refine is \c true, one step of iterative refinement is performed. The effort
    is much smaller than of the cholesky decomposition, so almost always worth it.  
    (Refinement is currently not implemented)

    \c x and \c b may be identical.

    \sa 
    Numerical Recipes Chapter 2.9
*/
void choleskySolve (XymVector& x, const XymMatrixVC& A, const XymVectorV& b, 
                    const XymMatrixVC& U=XymMatrixVC(), bool refine=true);

//! Inverts a symmetric positive definite matrix \c A
/*!    
    Refinement is currently not implemented.

    If \c U is not valid (for instance passing \c XymMatrixVC()) \c cholesky is called.

    \sa 
    Numerical Recipes Chapter 2.9

    Caution: This function inverts a SPD matrix using Cholesky decomposition, as opposed
    to \c inverseCholeskyFactor that computes the inverse of the Cholesky factor U^T.
 */
void choleskyInvert (XymMatrixC& result, const XymMatrixVC& A, 
                     const XymMatrixVC& U=XymMatrixVC(), bool refine=false);

//! Computes the inverse \f$ L^{-1} \f$ of the Cholesky factor \f$ L \f$, with \f$ LL^T=A \f$
/*! 
    Refinement is currently not implemented.

    A popular use of this function is to scale the rows of a least square minimization system
    by \c result to incorporate that the uncertainty is the covariance A.

    Note, that this function does no invert \c A. (cf. \c choleskyInvert)
 */
void inverseCholeskyFactor (XymMatrixC& result, const XymMatrixVC& A, bool refine=false);



//! Eliminates column \c col of \c A (SPD) by a one rank update \f$ uu^T \f$
/*! This is like one step of the cholesky decomposition. The vector \c u
    is determined and returned in \c u. The eliminated row / column is set
    to exactly 0 (without rounding errors). If the diagonal entry is smaller
    than \c minDiagonal the column is assumed to be 0 with numerical roundoff
    errors. In this case it is checked, that the off-diagonal entries are
    small too.
*/
void choleskyStep (XymMatrixVC& A, int col, XymVectorV& u, double minDiagonal=0);


//! Computes \c result to be \f$ UU^T \f$ for an upper triangular matrix \c U
void multiplyLTL (XymMatrixC& result, const XymMatrixVC& U);

//! Computes \c result to be \f$ J^TJ \f$ for any matrix \c J
void multiplyJTJ (XymMatrixVC& result, const XymMatrixVC& J);


//! Inverts and transposes an upper triangular matrix \c U.
/*! The result is upper triangular again. 
    \c refine is currently not implemented. The inverse is tranposed, since this
    makes access both for this routine and for using this routine to invert a
    symmetric matrix easier.
 */
void invertUT (XymMatrixC& result, const XymMatrixVC& U, bool refine=false);


//! Modify the inverse of a symmetric matrix after a rank one update
//! (sherman-morison formula)
/*! Let \c AI be the inverse of \c A. Than \c AI is modified to be the
    \f$ (A + factor uu^T)^{-1} \f$. This is done by application of the
    Sherman - Morrison Formula.

    \sa 
    Numerical Recipes Chapter 2.7
 */
void shermanMorisonUpdate (XymMatrixVC& AI, XymVectorV& u, double factor=1);

//@}



#endif
