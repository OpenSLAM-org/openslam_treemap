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

#ifndef XYMBLAS_H
#define XYMBLAS_H

//! Wrapper routines for BLAS
/*!\file xymBLAS.h 
   \author Udo Frese 

   Contains wrapper routines for BLAS (Basic Linear Algebra
   Subroutines) Only the functions containing dense double vectors and
   matrices are provided. Eventually this file must be compiled with
   the flag \c UNDER and / or \c CAPS set depending on the linker
   convention of the FORTRAN compiler used for the actual BLAS
   implementation. Since FORTRAN / BLAS use 1 as starting index and
   the classes in this package use 0, indices are converted.

   The wrapper functions are documented rather poorly. Use the BLAS
   documentation. The names correspond to a prefix \c xym and the BLAS
   name with capital letters. The leading 'D' is ommited, since C++
   functions can be overloaded. Only those functions working on double
   vectors and matrices are wrapped.

   Normally the program must be linked with (order is important)

  \code
      -lBLAS -lF77 -lm
  \endcode

  \sa
      http://www.netlib.org/blas/
*/

#include <xymatrix/xymVectorV.h>
#include <xymatrix/xymMatrixVC.h>




//************************************************************
/*!\name Level 1 BLAS
  Routines working on vectors. (\f$ O(n) \f$)
*/

//@{

//! Apply plane rotation (Givens) to x and y
void xymROT (XymVectorV& x, XymVectorV& y, double c, double s);

//! Apply modified plane rotation to x and y
/*! \warning Not implemented in current BLAS version. */
void xymROTM (XymVectorV& x, XymVectorV& y, double param[5]);

//! Swaps x and y
void xymSWAP (XymVectorV& x, XymVectorV& y);

//! Scales \c x by a factor of \c alpha
void xymSCAL (double alpha, XymVectorV& x);

//! Copies x to y
void xymCOPY (const XymVectorV& x, XymVectorV& y);

//! Copies x to y
/*! If y is uninitialised it is initialised to the right dimension
 */
void xymCOPY (const XymVectorV& x, XymVector& y);

//! \c y=alpha*x+y */
void xymAXPY (double alpha, const XymVectorV& x, XymVectorV& y);


//! returns \c x^ty
double xymDOT (const XymVectorV& x, const XymVectorV& y);

//! returns the 2 norm of \c v
double xymNRM2 (const XymVectorV& v);

//! returns the 1 norm of \c v
double xymASUM (const XymVectorV& v);

//! returns which entry has the maximum norm in \c v
/*! The index is converted to "start with 0" convention. 
*/
int xymAMAX (const XymVectorV& v);

//@}


//************************************************************
/*!\name Level 2 BLAS
  Routines working on vectors and matrices. (\f$ O(n^2) \f$)
*/
//@{

//! Whether to tranpose a matrix or not
enum XymTranspose {NOTRANSPOSE, TRANSPOSE};

//! Whether to use the upper or lower part (or both) of a symmetric matrix
enum XymUpperLower {UPPER, LOWER, BOTH};

//! Whether the diagonal can assumed to be 1
enum XymDiagonal {UNITDIAGONAL, NONUNITDIAGONAL};

//! Multiply a matrix with a vector
/*! Sets \c y=alpha*A*x+beta*y if (\c transpose==NOTRANSPOSE) or \c
    y=alpha*A^T*x+beta*y if (\c transpose==TRANSPOSE) 
*/
void xymGEMV (XymTranspose transpose, double alpha, 
              const XymMatrixVC& A, const XymVectorV& x, double beta, XymVectorV& y);

//! Multiply a symmetric matrix with a vector
/*! Sets \c y=alpha*A*x+beta*y if (\c transpose==false) or \c
    y=alpha*A^T*x+beta*y if (\c transpose==true) for a symmetric
    matrix A.  \c upLo specifies, whether the upper or lower part of
    \c A shall actually be used. (So it can sometimes be avoided to
    copy each entry of a symmetric matrix. 
*/
void xymSYMV (XymUpperLower upLo, double alpha, const XymMatrixVC& A, const XymVectorV& x, double beta, XymVectorV& y);


//! Multiply a triangular matrix with a vector
/*! Sets \c v=A*v if \c transpose==NOTRANSPOSE or \c v=A^T*v if \c
    transpose==TRANSPOSE for triangular matrix A.  \c upLo specifies,
    whether the matrix is upper or lower trianguler. \c diag
    specifies, whether the matrix diagonal shall be assumed to be 1
    (\c UNITDIAGONAL) or taken as in the matrix.
 */
void xymTRMV (XymUpperLower upLo, XymTranspose transpose, XymDiagonal diag, const XymMatrixVC& A, const XymVectorV& v);


//! Solve a triangular equation system.
/*! Sets \c v=A^-1*v if \c transpose==NOTRANSPOSE or \c v=A^-1T*v if
    \c transpose==TRANSPOSE for triangular matrix A.  \c upLo
    specifies, whether the matrix is upper or lower trianguler. \c
    diag specifies, whether the matrix diagonal shall be assumed to be
    1 (\c UNITDIAGONAL) or taken as in the matrix.
 */
void xymTRSV (XymUpperLower upLo, XymTranspose transpose, XymDiagonal diag, const XymMatrixVC& A, const XymVectorV& v);


//! Add a rank 1 update to a matrix
/*! Add \c alpha*xy^T to \c A.
 */
void xymGER (double alpha, const XymVectorV& x, const XymVectorV& y, XymMatrixVC& A);

//! Add a rank 1 update to a symmetric matrix
/*! Add \c alpha*xx^T to \c A. \c upLo specifies, whether to update the lower, upper or
    both parts of the matrix. The latter is and xy matrix extension. In this case the
    result is copied using a non BLAS routine. This is always slower and may be much
    slower if an optimized BLAS is present.
 */
void xymSYR (XymUpperLower upLo, double alpha, const XymVectorV& x, XymMatrixVC& A);

//! Add a rank 2 update to a symmetric matrix
/*! Add \c alpha*x*y^T to \c A. \c upLo specifies, whether to update the lower, upper or
    both parts of the matrix. The latter is and xy matrix extension. In this case the
    result is copied using a non BLAS routine. This is always slower and may be much
    slower if an optimized BLAS is present.
 */
void xymSYR2 (XymUpperLower upLo, double alpha, const XymVectorV& x, const XymVectorV& y, XymMatrixVC& A);

//@}


//************************************************************
/*!\name Level 3 BLAS
  Routines working on matrices. (\f$ O(n^3) \f$)
*/
//@{

//! Matrix matrix multiply
/*! Sets \c C=alpha*A*B+beta*C, where \c A and \c B are transposed before if this is
    specified by \c transA and \c transB.
*/
void xymGEMM (XymTranspose transA, XymTranspose transB, double alpha, 
              const XymMatrixVC& A, const XymMatrixVC& B, double beta, XymMatrixVC& C);


//! Matrix matrix multiply for a symmetric matrix
/*! Sets \c C=alpha*A*B+beta*C if \c bAInsteadAB==false or \c
    C=alpha*B*A=beta*C if \c bAInsteadAB==true. \c A must be
    symmetric. \c upLo specifies, which part of \c A actually to use.
 */
void xymSYMM (bool bAInsteadAB, XymUpperLower upLo, double alpha, 
              const XymMatrixVC& A, const XymMatrixVC& B, double beta, XymMatrixVC& C);


//! Symmetric rank k update
/*! Sets \c C=alpha*AA^T+beta*C if \c trans==NOTRANSPOSE or \c
    C=alpha*A^TA+beta*C if \c trans==TRANSPOSE. \c C must be
    symmetric.  \c upLo specifies, whether to update the upper or
    lower part of \c C. If \c upLo==BOTH, one part is copied not using
    BLAS, which may take some time.
*/
void xymSYRK (XymUpperLower upLo, XymTranspose trans, double alpha, 
              const XymMatrixVC& A, double beta, XymMatrixVC& C);


//! Symmetric rank 2*k update
/*! Sets \c C=alpha*AB^T+alpha*B*A^T+beta*C if \c trans==NOTRANSPOSE
    or \c C=alpha*A^TB+alpha*B^TA+beta*C if \c trans==TRANSPOSE. \c C
    must be symmetric.  \c upLo specifies, whether to update the upper
    or lower part of \c C. If \c upLo==BOTH, one part is copied not
    using BLAS, which may take some time.
*/
void xymSYRK2 (XymUpperLower upLo, XymTranspose trans, double alpha, 
              const XymMatrixVC& A, const XymMatrixVC& B, double beta, XymMatrixVC& C);


//! Matrix matrix multiply for a triangular matrix
/*! Sets \c C=alpha*A*B if \c bAInsteadAB==false or \c C=alpha*B*A if
    \c bAInsteadAB==true. If \c transA==TRANSPOSE \c A is transposed
    before. \c upLo specifies, whether \c A is an upper or lower
    triangular matrix.  \c diag specifies, whether the diagonal is to
    be assumed to be 1 (\c UNITDIAGONAL) or taken as in \c A (\c
    NONUNITDIAGONAL).
*/
void xymTRMM (bool bAInsteadAB, XymUpperLower upLo, XymTranspose transA, XymDiagonal diag,
              double alpha, const XymMatrixVC& A, XymMatrixVC& B);

//! Equation solving for a triangular matrix and a matrix
/*! Sets \c C=alpha*A^-1*B if \c bAInsteadAB==false or \c
    C=alpha*B*A^-1 if \c bAInsteadAB==true. If \c transA==TRANSPOSE \c
    A is transposed before. \c upLo specifies, whether \c A is an
    upper or lower triangular matrix.  \c diag specifies, whether the
    diagonal is to be assumed to be 1 (\c UNITDIAGONAL) or taken as in
    \c A (\c NONUNITDIAGONAL).
*/
void xymTRSM (bool bAInsteadAB, XymUpperLower upLo, XymTranspose transA, XymDiagonal diag,
              double alpha, const XymMatrixVC& A, XymMatrixVC& B);



//@}


//! Copies the upper half of \c A into the lower half.
void upperToLower (XymMatrixVC& A);

#endif
