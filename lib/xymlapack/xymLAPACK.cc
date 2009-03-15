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

//! \file xymLAPACK.cc implementation of xymLAPACK.h

#include "xymLAPACK.h"
//#include "xymBLAS.h"

#include "clapack.h"

#include <xymatrix/xymOperations.h>

void xymSYGV (const XymMatrixC& A, const XymMatrixC& B, XymMatrixC& Z, XymVector& w)
    throw (XymLAPACKException)
{
    checkSquare (A);
    checkSameFormat (A, B);
    
    int n = A.rows();
    copyFrom (Z, A);
    w.create (n);

    int nb = 20; // Assume LAPACK block size
    int nWork = (nb+2)*n;
    double* work = new double[nWork];

    int info;
    
    XymMatrixC BCopy(B);
    char ulC = 'U';
    char jobZC = 'V';
    int typeOfProblem = 1;
    dsygv_ (&typeOfProblem, &jobZC, &ulC, // Type Ax = lambda *Bx, compute eigenvectors, use upper triangle
            &n, Z.base(), (int*) Z.ld(), BCopy.base(), (int*) BCopy.ld(), 
            w.base(), work, &nWork, &info); 
    // 1, 1 are the length of ulC and jobZC

    delete[] work;

    if (info<0) throw XymLAPACKException ("Illegal parameter calling dsygv (FORTRAN)");
    if (info>0) throw XymLAPACKException ("B is not SPD");
}


void xymGEQR2 (const XymMatrixC& A, XymMatrixC& R, XymVector& work)
{
  int m = A.rows();
  int n = A.cols();  
  int rDim;
  if (m<n) rDim = m;
  else rDim = n;
  work.resize (rDim+n, false);  
  int info;

  if (&A!=&R) R = A;  
  if (R.rows()>0 && R.cols()>0) dgeqr2_ (&m, &n, R.base(), (int*) R.ld(), work.base(), work.base()+rDim, &info);
  if (rDim<R.rows()) R.deleteLastRow (R.rows()-rDim);
  zeroLower (R);

  if (info<0) throw XymLAPACKException ("Illegal parameter calling dgeqr2 (FORTRAN)");
}


void xymGEQR2 (const XymMatrixC& A, XymMatrixC& R)
{
  XymVector work;
  xymGEQR2 (A, R, work);
}

