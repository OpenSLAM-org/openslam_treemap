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


//! \file xymBLAS.cc implementation of \c xymBLAS.h
#include "xymBLAS.h"
#include <xymatrix/xymGeneral.h>

extern "C" 
{
#include "cblas.h"
}


// Level 1 BLAS

void xymROT (XymVectorV& x, XymVectorV& y, double c, double s)
{
    checkEqualLength (x, y);
    int sz  = x.size();
    drot_ (&sz, x.base(), x.st(), y.base(), y.st(), &c, &s);
}


void xymROTM (XymVectorV& x, XymVectorV& y, double param[5])
{
    checkEqualLength (x, y);
    int sz  = x.size();
//    drotm_ (&sz, x.base(), x.st(), y.base(), y.st(), param);
}


void xymSWAP (XymVectorV& x, XymVectorV& y)
{
    checkEqualLength (x, y);
    int sz  = x.size();
    dswap_ (&sz, x.base(), x.st(), y.base(), y.st());    
}


void xymSCAL (double alpha, XymVectorV& v)
{
    int sz  = v.size();
    dscal_ (&sz, &alpha, v.base(), v.st());
}


void xymCOPY (const XymVectorV& x, XymVectorV& y)
{
    checkEqualLength (x, y);
    int sz  = x.size();
    dcopy_ (&sz, x.base(), x.st(), y.base(), y.st());
}

void xymCOPY (const XymVectorV& x, XymVector& y)
{
    if (!y.isValid()) y.create (x.size());
    else checkEqualLength (x, y);
    int sz  = x.size();
    dcopy_ (&sz, x.base(), x.st(), y.base(), y.st());
}


void xymAXPY (double alpha, const XymVectorV& x, XymVectorV& y)
{
    checkEqualLength (x, y);
    int sz  = x.size();
    daxpy_ (&sz, &alpha, x.base(), x.st(), y.base(), y.st());
}


double xymDOT (const XymVectorV& x, const XymVectorV& y)
{
    checkEqualLength (x, y);
    int sz  = x.size();
    return ddot_ (&sz, x.base(), x.st(), y.base(), y.st());
}


double xymNRM2 (const XymVectorV& v)
{
    int sz  = v.size();
    return dnrm2_ (&sz, v.base(), v.st());
}


double xymASUM (const XymVectorV& v)
{
    int sz  = v.size();
    return dasum_ (&sz, v.base(), v.st());
}
 

int xymAMAX (const XymVectorV& v)
{
    int sz  = v.size();
    return idamax_ (&sz, v.base(), v.st())-1;
}


// Level 2 BLAS
void xymGEMV (XymTranspose transpose, double alpha, const XymMatrixVC& A, const XymVectorV& x, double beta, XymVectorV& y)
{
    if (transpose==NOTRANSPOSE) checkForMultiply (y, A, x);
    else checkForMultiply (x, A, y);
    int m = A.rows();
    int n = A.cols();
    char tC = (transpose==TRANSPOSE)?'T':'N';
    dgemv_ (&tC, &m, &n, &alpha, A.base(), A.ld(), x.base(), x.st(), &beta, y.base(), y.st(), 1); 
    // 1 is length of 'tC'
}


void xymSYMV (XymUpperLower part, double alpha, const XymMatrixVC& A, const XymVectorV& x, double beta, XymVectorV& y)
{
    checkForMultiply (y, A, x);
    checkSquare (A);
    int n = A.cols();
    char ulC = (part==LOWER)?'L':'U';
    dsymv_ (&ulC, &n, &alpha, A.base(), A.ld(), x.base(), x.st(), &beta, y.base(), y.st(), 1); 
    // 1 is length of 'ulC'
}


void xymTRMV (XymUpperLower upLo, XymTranspose transpose, XymDiagonal diag, const XymMatrixVC& A, const XymVectorV& v)
{
    checkForMultiply (v, A, v);
    checkSquare (A);
    int n = A.cols();
    char ulC = (upLo==LOWER)?'L':'U';
    char tpC = (transpose==TRANSPOSE)?'T':'N';
    char dC = (diag==UNITDIAGONAL)?'U':'N';
    dtrmv_ (&ulC, &tpC, &dC, &n, A.base(), A.ld(), v.base(), v.st(), 1, 1, 1); 
    // 1, 1, 1 is length of 'ulC', 'tpC', 'dC'
}


void xymTRSV (XymUpperLower upLo, XymTranspose transpose, XymDiagonal diag, const XymMatrixVC& A, const XymVectorV& v)
{
    checkForMultiply (v, A, v);
    checkSquare (A);
    int n = A.cols();
    char ulC = (upLo==LOWER)?'L':'U';
    char tpC = (transpose==TRANSPOSE)?'T':'N';
    char dC = (diag==UNITDIAGONAL)?'U':'N';
    dtrsv_ (&ulC, &tpC, &dC, &n, A.base(), A.ld(), v.base(), v.st(), 1, 1, 1);     
     // 1, 1, 1 is length of 'ulC', 'tpC', 'dC'
}


void xymGER (double alpha, const XymVectorV& x, const XymVectorV& y, XymMatrixVC& A)
{
    int m, n;
    m = A.rows();
    n = A.cols();
    XYMFORMAT (m==x.size());
    XYMFORMAT (n==y.size());
    dger_ (&m, &n, &alpha, x.base(), x.st(), y.base(), y.st(), A.base(), A.ld());
}


void xymSYR (XymUpperLower upLo, double alpha, const XymVectorV& x, XymMatrixVC& A)
{
    checkSquare (A);
    int n;
    n = A.cols();
    XYMFORMAT (n==x.size());
    char ulC = (upLo==LOWER)?'L':'U';
    dsyr_ (&ulC, &n, &alpha, x.base(), x.st(), A.base(), A.ld(), 1); 
    // 1 is length for 'ulC'
    if (upLo==BOTH) upperToLower (A);
}


void xymSYR2 (XymUpperLower upLo, double alpha, const XymVectorV& x, const XymVectorV& y, XymMatrixVC& A)
{
    checkSquare (A);
    int n;
    n = A.cols();
    XYMFORMAT (n==x.size());
    XYMFORMAT (n==y.size());
    char ulC = (upLo==LOWER)?'L':'U';
    dsyr2_ (&ulC, &n, &alpha, x.base(), x.st(), y.base(), y.st(), A.base(), A.ld(), 1); 
    // 1 is length for 'ulC'
    if (upLo==BOTH) upperToLower (A);
}


// Level 3 BLAS

void xymGEMM (XymTranspose transA, XymTranspose transB, double alpha, 
              const XymMatrixVC& A, const XymMatrixVC& B, double beta, XymMatrixVC& C)
{
    int m = C.rows();
    int n = C.cols();
    int k; 
    if (transA==NOTRANSPOSE) {
        XYMFORMAT (A.rows()==m);
        k = A.cols();
    }
    else {
        XYMFORMAT (A.cols()==m);
        k = A.rows();
    }
    if (transB==NOTRANSPOSE) {
        XYMFORMAT (B.rows()==k);
        XYMFORMAT (B.cols()==n);
    }
    else {
        XYMFORMAT (B.cols()==k);
        XYMFORMAT (B.rows()==n);
    }
    
    char trAC = (transA==TRANSPOSE)?'T':'N';
    char trBC = (transB==TRANSPOSE)?'T':'N';

    dgemm_ (&trAC, &trBC, &m, &n, &k, &alpha, A.base(), A.ld(), B.base(), B.ld(), &beta, C.base(), C.ld(), 1, 1);
    // 1 is length for 'trAC' and 'trBC'
}


void xymSYMM (bool bAInsteadAB, XymUpperLower upLo, double alpha, 
              const XymMatrixVC& A, const XymMatrixVC& B, double beta, XymMatrixVC& C)
{
    checkSquare (A);
    int m = C.rows();
    int n = C.cols();
    if (bAInsteadAB) checkForMultiply (C, A, B);
    else checkForMultiply (C, B, A);
    
    char ulC = (upLo==LOWER)?'L':'U';
    char sdC = (bAInsteadAB)?'R':'L';
    
    dsymm_ (&sdC, &ulC, &m, &n, &alpha, A.base(), A.ld(), B.base(), B.ld(), &beta, C.base(), C.ld(), 1, 1);
    // 1 is length for 'ulC' and 'sdC'
}


void xymSYRK (XymUpperLower upLo, XymTranspose trans, double alpha, 
              const XymMatrixVC& A, double beta, XymMatrixVC& C)
{
    checkSquare(C);
    int  n = C.cols();
    int  k;
    if (trans==NOTRANSPOSE) { // AA^T
        XYMFORMAT (n==A.rows());
        k = A.cols();
    }
    else { // A^TA
        XYMFORMAT (n==A.cols());
        k = A.rows();
    }
    char ulC = (upLo==LOWER)?'L':'U';
    char trC = (trans==TRANSPOSE)?'T':'N';
    dsyrk_ (&ulC, &trC, &n, &k, &alpha, A.base(), A.ld(), &beta, C.base(), C.ld(), 1, 1);
    // 1 is length for 'ulC' and 'trC'
    if (upLo==BOTH) upperToLower(C);
}


void xymSYRK2 (XymUpperLower upLo, XymTranspose trans, double alpha, 
               const XymMatrixVC& A, const XymMatrixVC& B, double beta, XymMatrixVC& C)
{
    checkSquare(C);
    int  n = C.cols();
    int  k;
    if (trans==NOTRANSPOSE) { // AB^T+B^TA
        XYMFORMAT (n==A.rows());
        XYMFORMAT (n==B.rows());
        XYMFORMAT (A.cols()==B.cols());
        k = A.cols();
    }
    else { // A^TB+AB^T
        XYMFORMAT (n==A.cols());
        XYMFORMAT (n==B.cols());
        XYMFORMAT (A.rows()==B.rows());
        k = A.rows();
    }
    char ulC = (upLo==LOWER)?'L':'U';
    char trC = (trans==TRANSPOSE)?'T':'N';
    dsyr2k_ (&ulC, &trC, &n, &k, &alpha, A.base(), A.ld(), B.base(), B.ld(), &beta, C.base(), C.ld(), 1, 1);
    // 1 is length for 'ulC' and 'trC'
    if (upLo==BOTH) upperToLower(C);
}


void xymTRMM (bool bAInsteadAB, XymUpperLower upLo, XymTranspose transA, XymDiagonal diag,
              double alpha, const XymMatrixVC& A, XymMatrixVC& B)
{
    checkSquare(A);
    int m = B.rows();
    int n = B.cols();
    if (bAInsteadAB) XYMFORMAT(A.cols()==n);
    else XYMFORMAT(A.rows()==m);
    char ulC = (upLo==LOWER)?'L':'U';
    char trC = (transA==TRANSPOSE)?'T':'N';
    char sdC = (bAInsteadAB)?'R':'L';
    char dC = (diag==UNITDIAGONAL)?'U':'N';
    dtrmm_ (&sdC, &ulC, &trC, &dC, &m, &n, &alpha, A.base(), A.ld(), B.base(), B.ld(), 1, 1, 1, 1);
    // 1, 1, 1, 1 is length for 'sdC', 'ulC', 'trC', 'dC'
}


void xymTRSM (bool bAInsteadAB, XymUpperLower upLo, XymTranspose transA, XymDiagonal diag,
              double alpha, const XymMatrixVC& A, XymMatrixVC& B)
{
    checkSquare(A);
    int m = B.rows();
    int n = B.cols();
    if (bAInsteadAB) XYMFORMAT(A.cols()==n);
    else XYMFORMAT(A.rows()==m);
    char ulC = (upLo==LOWER)?'L':'U';
    char trC = (transA==TRANSPOSE)?'T':'N';
    char sdC = (bAInsteadAB)?'R':'L';
    char dC = (diag==UNITDIAGONAL)?'U':'N';
    dtrsm_ (&sdC, &ulC, &trC, &dC, &m, &n, &alpha, A.base(), A.ld(), B.base(), B.ld(), 1, 1, 1, 1);
    // 1, 1, 1, 1 is length for 'sdC', 'ulC', 'trC', 'dC'
}


// Auxiliary functions

void upperToLower (XymMatrixVC& A)
{
    int n      = A.rows();
    int srcOfs = A.colOfs();
    int dstOfs = A.rowOfs();
    for (int i=0; i<n-1; i++) {
        // Copy row 'i' to column 'i' starting with entry 'i+1'
        double* src = &(A(i,i+1));
        double* dst = &(A(i+1,i));
        double* end = src+(n-i-1);
        while (src!=end) {
            *dst = *src;
            src += srcOfs;
            dst += dstOfs;
        }       
    }
}


