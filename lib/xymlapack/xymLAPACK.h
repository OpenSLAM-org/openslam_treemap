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

#ifndef XYMLAPACK_H
#define XYMLAPACK_H

//! \file xymLAPACK.h Wrapper routines for LAPACK
/*! \file xymLAPACK.h
    \author Udo Frese

   Contains wrapper routines for LAPACK. Only the functions containing
   dense double vectors and matrices are provided. Eventually this
   file must be compiled with the flag \c UNDER and / or \c CAPS set
   depending on the linker convention of the FORTRAN compiler used for
   the actual BLAS implementation. Since FORTRAN / BLAS / LAPACK use 1
   as starting index and the classes in this package use 0, indices
   are converted.

   The wrapper functions are documented rather poorly. Use the LAPACK
   documentation. The names correspond to a prefix \c xym and the LAPACK
   name with capital letters. The leading 'D' is ommited, since C++
   functions can be overloaded. Only those functions working on double
   vectors and matrices are wrapped.

   Currently only very few functions are implemented. Add further
   functions as necessary.

   Normally the program must be linked with (order is important)

  \code
     -llapack -lblas -lF77
  \endcode
*/

#include <xymatrix/xymMatrixC.h>
#include <xymatrix/xymVector.h>

#include "clapack.h"

#include <stdexcept>

//! Exception thrown by LAPACK wrappers
class XymLAPACKException : public runtime_error
{
public:
  //! Std constructor taking a message
    XymLAPACKException (const string msg)
            :runtime_error (msg)
        {};
};



//! Computes an extended eigenvalue problem for \c A*x = lambda*B*x
/*! \c A and \c B must be symmetric and \c B must be SPD. The
    eigenvalues are returned in \c w in ascending order. The
    eigenvectors are returned in the columns of \c Z. It holds: \c
    Z^TBZ=I and \c Z^TAZ=diag{w}.

    Some copying is performed.
*/
void xymSYGV (const XymMatrixC& A, const XymMatrixC& B, XymMatrixC& Z, XymVector& w)
    throw (XymLAPACKException);


//! Computes an QR factorization of a rectangular matrix R
/*! \c A can be a general matrix, both rectangular and rank deficient.
    The returned matrix R is triangular with dimension min(A.cols(), A.rows()).
    It may be trapezoidal if R is rank deficient. The matrix Q is implicitly
    computed by the LAPACK routine and currently thrown away.

    A and R can be identical.
    Some copying is performed.
    The routine needs a workspace \c work. If \c work is capable of holding enough
    space it is used otherwise reallocated. This way memory allocation can be avoided.
 */
void xymGEQR2 (const XymMatrixC& A, XymMatrixC& R, XymVector& work);

//! Overloaded function
/*! The workspace is allocated and deallocated. */
void xymGEQR2 (const XymMatrixC& A, XymMatrixC& R);



#endif  /* XYMLAPACK_H */
