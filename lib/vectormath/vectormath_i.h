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
/*! \file vectormath_i.h
    \brief Implementation of \c vectormath.h
    \author Udo Frese

    Implementation of the inline functions for small matrix computations.
*/
#include<limits>
#include <math.h>

inline void vmSet (VmVector2& v, VmReal v0, VmReal v1)
{
    v[0] = v0;
    v[1] = v1;
}



inline void vmSet (VmVector3& v, VmReal v0, VmReal v1, VmReal v2)
{
    v[0] = v0;
    v[1] = v1;
    v[2] = v2;
}



inline void vmSet (VmMatrix2x2& A, VmReal a00, VmReal a01, VmReal a10, VmReal a11)
{
    A[0][0] = a00;
    A[0][1] = a01;
    A[1][0] = a10;
    A[1][1] = a11;
}




inline void vmSet (VmMatrix2x3& A, VmReal a00, VmReal a01, VmReal a02, VmReal a10, VmReal a11, VmReal a12)
{
    A[0][0] = a00;
    A[0][1] = a01;
    A[0][2] = a02;
    A[1][0] = a10;
    A[1][1] = a11;
    A[1][2] = a12;
}



inline void vmSet (VmMatrix3x2& A, VmReal a00, VmReal a01, VmReal a10, VmReal a11, VmReal a20, VmReal a21)
{
    A[0][0] = a00;
    A[0][1] = a01;
    A[1][0] = a10;
    A[1][1] = a11;
    A[2][0] = a20;
    A[2][1] = a21;
}



inline void vmSet (VmMatrix3x3& A, VmReal a00, VmReal a01, VmReal a02, VmReal a10, VmReal a11, VmReal a12, VmReal a20, VmReal a21, VmReal a22)
{
    A[0][0] = a00;
    A[0][1] = a01;
    A[0][2] = a02;
    A[1][0] = a10;
    A[1][1] = a11;
    A[1][2] = a12;
    A[2][0] = a20;
    A[2][1] = a21;
    A[2][2] = a22;
}


inline void vmZero (VmVector2& v)
{
    v[0] = v[1] = 0;
}


inline void vmZero (VmVector3& v)
{
    v[0] = v[1] = v[2] = 0;
}



inline void vmZero (VmMatrix2x2& A)
{
    A[0][0] = A[0][1] = 
      A[1][0] = A[1][1] = 0;
}



inline void vmZero (VmMatrix2x3& A)
{
    A[0][0] = A[0][1] = A[0][2] = 
      A[1][0] = A[1][1] = A[1][2] = 0;
}



inline void vmZero (VmMatrix3x2& A)
{
    A[0][0] = A[0][1] =
      A[1][0] = A[1][1] =
      A[2][0] = A[2][1] = 0;    
}



inline void vmZero (VmMatrix3x3& A)
{
    A[0][0] = A[0][1] = A[0][2] = 
      A[1][0] = A[1][1] = A[1][2] = 
      A[2][0] = A[2][1] = A[2][2] = 0;    
}


inline void vmZero (VmMatrix4x4& A)
{
    A[0][0] = A[0][1] = A[0][2] = A[0][3] = 
      A[1][0] = A[1][1] = A[1][2] = A[1][3] = 
      A[2][0] = A[2][1] = A[2][2] = A[2][3] =     
      A[3][0] = A[3][1] = A[3][2] = A[3][3] = 0;      
}


inline void vmNaN (VmVector2& v)
{
  v[0] = v[1] = vmNaN();
}


inline void vmNaN (VmVector3& v)
{
  v[0] = v[1] = v[2] = vmNaN();
}



inline void vmNaN (VmMatrix2x2& A)
{
  A[0][0] = A[0][1] = 
    A[1][0] = A[1][1] = vmNaN();
}



inline void vmNaN (VmMatrix2x3& A)
{
    A[0][0] = A[0][1] = A[0][2] =
      A[1][0] = A[1][1] = A[1][2] = vmNaN();
}



inline void vmNaN (VmMatrix3x2& A)
{
  A[0][0] = A[0][1] = 
    A[1][0] = A[1][1] = 
    A[2][0] = A[2][1] = vmNaN();    
}



inline void vmNaN (VmMatrix3x3& A)
{
    A[0][0] = A[0][1] = A[0][2] = 
    A[1][0] = A[1][1] = A[1][2] = 
      A[2][0] = A[2][1] = A[2][2] = vmNaN();    
}


inline void vmNaN (VmMatrix4x4& A)
{
    A[0][0] = A[0][1] = A[0][2] = A[0][3] = 
    A[1][0] = A[1][1] = A[1][2] = A[1][3] = 
    A[2][0] = A[2][1] = A[2][2] = A[2][3] =     
      A[3][0] = A[3][1] = A[3][2] = A[3][3] = vmNaN ();      
}





inline void vmOne (VmMatrix2x2& A)
{
  A[0][0] = A[1][1] = 1;
  A[0][1] = A[1][0] = 0;
}


inline void vmOne (VmMatrix3x3& A)
{
  A[0][0] = A[1][1] = A[2][2] = 1;
  A[0][1] = A[0][2] = A[1][2] = 0;
  A[1][0] = A[2][0] = A[2][1] = 0;
}


inline void vmOne (VmMatrix4x4& A)
{
    A[0][1] = A[0][2] = A[0][3] = 0;
    A[1][0] = A[1][2] = A[1][3] = 0;
    A[2][0] = A[2][1] = A[2][3] = 0;    
    A[3][0] = A[3][1] = A[3][2] = 0;      
    A[0][0] = A[1][1] = A[2][2] = A[3][3] = 1;    
}


inline bool isFinite (const VmVector3& x)
{
  return isfinite(x[0]) && isfinite(x[1]) && isfinite(x[2]);
}


inline bool isFinite (const VmVector2& x)
{
  return isfinite(x[0]) && isfinite(x[1]) && isfinite(x[2]);
}


inline bool isFinite (const VmMatrix2x2& A)
{
  return isfinite(A[0][0]) && isfinite(A[0][1]) &&
    isfinite(A[1][0]) && isfinite(A[1][1]);
}

inline bool isFinite (const VmMatrix2x3& A)
{
  return isfinite(A[0][0]) && isfinite(A[0][1]) && isfinite(A[0][2]) &&
    isfinite(A[1][0]) && isfinite(A[1][1]) && isfinite(A[1][2]);
}


inline bool isFinite (const VmMatrix3x2& A)
{
  return isfinite(A[0][0]) && isfinite(A[0][1]) &&
    isfinite(A[1][0]) && isfinite(A[1][1]) &&
    isfinite(A[2][0]) && isfinite(A[2][1]);
}


inline bool isFinite (const VmMatrix3x3& A)
{  
  return isfinite(A[0][0]) && isfinite(A[0][1]) && isfinite(A[0][2]) &&
    isfinite(A[1][0]) && isfinite(A[1][1]) && isfinite(A[1][2]) && 
  isfinite(A[2][0]) && isfinite(A[2][1]) && isfinite(A[2][2]);
}



inline void vmCopy (const VmMatrix2x2& A, VmMatrix2x2& B)
{
    B[0][0] = A[0][0];
    B[0][1] = A[0][1];
    B[1][0] = A[1][0];
    B[1][1] = A[1][1];
}




inline void vmCopy (const VmMatrix3x3& A, VmMatrix3x3& B)
{
    B[0][0] = A[0][0];
    B[0][1] = A[0][1];
    B[0][2] = A[0][2];
    B[1][0] = A[1][0];
    B[1][1] = A[1][1];
    B[1][2] = A[1][2];
    B[2][0] = A[2][0];
    B[2][1] = A[2][1];
    B[2][2] = A[2][2];
}


inline void vmCopy (const VmMatrix2x3& A, VmMatrix2x3& B)
{
    B[0][0] = A[0][0];
    B[0][1] = A[0][1];
    B[0][2] = A[0][2];
    B[1][0] = A[1][0];
    B[1][1] = A[1][1];
    B[1][2] = A[1][2];
}


inline void vmCopy (const VmMatrix4x4& A, VmMatrix4x4& B)
{
  B[0][0] = A[0][0];
  B[0][1] = A[0][1];
  B[0][2] = A[0][2];
  B[0][3] = A[0][3];
  B[1][0] = A[1][0];
  B[1][1] = A[1][1];
  B[1][2] = A[1][2];
  B[1][3] = A[1][3];
  B[2][0] = A[2][0];
  B[2][1] = A[2][1];
  B[2][2] = A[2][2];
  B[2][3] = A[2][3];
  B[3][0] = A[3][0];
  B[3][1] = A[3][1];
  B[3][2] = A[3][2];
  B[3][3] = A[3][3];
}




inline void vmCopy (const VmVector3& a, VmVector3& b)
{
    b[0] = a[0];
    b[1] = a[1];
    b[2] = a[2];
}



inline void vmCopy (const VmVector2& a, VmVector2& b)
{
    b[0] = a[0];
    b[1] = a[1];
}



inline void vmTranspose (VmMatrix2x2& A)
{
    vmSwap(A[0][1], A[1][0]);
}



inline void vmTranspose (VmMatrix3x3& A)
{
    vmSwap(A[0][1], A[1][0]);
    vmSwap(A[0][2], A[2][0]);
    vmSwap(A[1][2], A[2][1]);
}



inline void vmTranspose (VmMatrix2x2& result, const VmMatrix2x2& A)
{
    result[0][0] = A[0][0];
    result[0][1] = A[1][0];
    result[1][0] = A[0][1];
    result[1][1] = A[1][1];
}

    
inline void vmTranspose (VmMatrix2x3& result, const VmMatrix3x2& A)
{
    result[0][0] = A[0][0];
    result[0][1] = A[1][0];
    result[0][2] = A[2][0];
    result[1][0] = A[0][1];
    result[1][1] = A[1][1];
    result[1][2] = A[2][1];
}



inline void vmTranspose (VmMatrix3x2& result, const VmMatrix2x3& A)
{
    result[0][0] = A[0][0];
    result[0][1] = A[1][0];
    result[1][0] = A[0][1];
    result[1][1] = A[1][1];
    result[2][0] = A[0][2];
    result[2][1] = A[1][2];
}

    
inline void vmTranspose (VmMatrix3x3& result, const VmMatrix3x3& A)
{
    result[0][0] = A[0][0];
    result[0][1] = A[1][0];
    result[0][2] = A[2][0];
    result[1][0] = A[0][1];
    result[1][1] = A[1][1];
    result[1][2] = A[2][1];
    result[2][0] = A[0][2];
    result[2][1] = A[1][2];
    result[2][2] = A[2][2];
}



inline void vmMultiply (VmVector2& result, const VmMatrix2x2& A, const VmVector2& v)
{
    result[0] = A[0][0]*v[0] + A[0][1]*v[1];
    result[1] = A[1][0]*v[0] + A[1][1]*v[1];
}



inline void vmMultiply (VmVector2& result, const VmMatrix2x3& A, const VmVector3& v)
{
    result[0] = A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2];
    result[1] = A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2];
}



inline void vmMultiply (VmVector3& result, const VmMatrix3x2& A, const VmVector2& v)
{
    result[0] = A[0][0]*v[0] + A[0][1]*v[1];
    result[1] = A[1][0]*v[0] + A[1][1]*v[1];
    result[2] = A[2][0]*v[0] + A[2][1]*v[1];
}



inline void vmMultiply (VmVector3& result, const VmMatrix3x3& A, const VmVector3& v)
{
    result[0] = A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2];
    result[1] = A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2];
    result[2] = A[2][0]*v[0] + A[2][1]*v[1] + A[2][2]*v[2];
}




inline void vmMultiplyAdd (VmVector2& result, const VmMatrix2x2& A, const VmVector2& v)
{
    result[0] += A[0][0]*v[0] + A[0][1]*v[1];
    result[1] += A[1][0]*v[0] + A[1][1]*v[1];
}



inline void vmMultiplyAdd (VmVector2& result, const VmMatrix2x3& A, const VmVector3& v)
{
    result[0] += A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2];
    result[1] += A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2];
}



inline void vmMultiplyAdd (VmVector3& result, const VmMatrix3x2& A, const VmVector2& v)
{
    result[0] += A[0][0]*v[0] + A[0][1]*v[1];
    result[1] += A[1][0]*v[0] + A[1][1]*v[1];
    result[2] += A[2][0]*v[0] + A[2][1]*v[1];
}



inline void vmMultiplyAdd (VmVector3& result, const VmMatrix3x3& A, const VmVector3& v)
{
    result[0] += A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2];
    result[1] += A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2];
    result[2] += A[2][0]*v[0] + A[2][1]*v[1] + A[2][2]*v[2];
}



inline void vmMultiplyAddScale (VmVector2& result, const VmMatrix2x2& A, const VmVector2& v, double factor)
{
    result[0] += factor*(A[0][0]*v[0] + A[0][1]*v[1]);
    result[1] += factor*(A[1][0]*v[0] + A[1][1]*v[1]);
}



inline void vmMultiplyAddScale (VmVector2& result, const VmMatrix2x3& A, const VmVector3& v, double factor)
{
    result[0] += factor*(A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2]);
    result[1] += factor*(A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2]);
}



inline void vmMultiplyAddScale (VmVector3& result, const VmMatrix3x2& A, const VmVector2& v, double factor)
{
    result[0] += factor*(A[0][0]*v[0] + A[0][1]*v[1]);
    result[1] += factor*(A[1][0]*v[0] + A[1][1]*v[1]);
    result[2] += factor*(A[2][0]*v[0] + A[2][1]*v[1]);
}



inline void vmMultiplyAddScale (VmVector3& result, const VmMatrix3x3& A, const VmVector3& v, double factor)
{
    result[0] += factor*(A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2]);
    result[1] += factor*(A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2]);
    result[2] += factor*(A[2][0]*v[0] + A[2][1]*v[1] + A[2][2]*v[2]);
}



inline void vmMultiplySub (VmVector2& result, const VmMatrix2x2& A, const VmVector2& v)
{
    result[0] -= A[0][0]*v[0] + A[0][1]*v[1];
    result[1] -= A[1][0]*v[0] + A[1][1]*v[1];
}



inline void vmMultiplySub (VmVector2& result, const VmMatrix2x3& A, const VmVector3& v)
{
    result[0] -= A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2];
    result[1] -= A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2];
}



inline void vmMultiplySub (VmVector3& result, const VmMatrix3x2& A, const VmVector2& v)
{
    result[0] -= A[0][0]*v[0] + A[0][1]*v[1];
    result[1] -= A[1][0]*v[0] + A[1][1]*v[1];
    result[2] -= A[2][0]*v[0] + A[2][1]*v[1];
}



inline void vmMultiplySub (VmVector3& result, const VmMatrix3x3& A, const VmVector3& v)
{
    result[0] -= A[0][0]*v[0] + A[0][1]*v[1] + A[0][2]*v[2];
    result[1] -= A[1][0]*v[0] + A[1][1]*v[1] + A[1][2]*v[2];
    result[2] -= A[2][0]*v[0] + A[2][1]*v[1] + A[2][2]*v[2];
}



inline void vmApplyTransformToPoint (VmVector3& result, const VmMatrix4x4& T, const VmVector3& p)
{
  result[0] = T[0][0]*p[0] + T[0][1]*p[1] + T[0][2]*p[2] + T[0][3];
  result[1] = T[1][0]*p[0] + T[1][1]*p[1] + T[1][2]*p[2] + T[1][3];
  result[2] = T[2][0]*p[0] + T[2][1]*p[1] + T[2][2]*p[2] + T[2][3];  
}


inline void vmApplyInverseTransformToPoint (VmVector3& result, const VmMatrix4x4& T, const VmVector3& p)
{
  VmVector3 v;
  v[0] = p[0]-T[0][3];
  v[1] = p[1]-T[1][3];
  v[2] = p[2]-T[2][3];
  result[0] = T[0][0]*v[0] + T[1][0]*v[1] + T[2][0]*v[2];
  result[1] = T[0][1]*v[0] + T[1][1]*v[1] + T[2][1]*v[2];
  result[2] = T[0][2]*v[0] + T[1][2]*v[1] + T[2][2]*v[2];  
}


inline void vmApplyTransformToDirection (VmVector3& result, const VmMatrix4x4& T, const VmVector3& d)
{
  result[0] = T[0][0]*d[0] + T[0][1]*d[1] + T[0][2]*d[2];
  result[1] = T[1][0]*d[0] + T[1][1]*d[1] + T[1][2]*d[2];
  result[2] = T[2][0]*d[0] + T[2][1]*d[1] + T[2][2]*d[2];  
}


inline void vmApplyInverseTransformToDirection (VmVector3& result, const VmMatrix4x4& T, const VmVector3& d)
{
  result[0] = T[0][0]*d[0] + T[1][0]*d[1] + T[2][0]*d[2];
  result[1] = T[0][1]*d[0] + T[1][1]*d[1] + T[2][1]*d[2];
  result[2] = T[0][2]*d[0] + T[1][2]*d[1] + T[2][2]*d[2];  
}


inline void vmInverseTransform (VmMatrix4x4& result, const VmMatrix4x4& A)
{
  result[0][0] = A[0][0];
  result[1][0] = A[0][1];
  result[2][0] = A[0][2];
  result[3][0] = 0;

  result[0][1] = A[1][0];
  result[1][1] = A[1][1];
  result[2][1] = A[1][2];
  result[3][1] = 0;

  result[0][2] = A[2][0];
  result[1][2] = A[2][1];
  result[2][2] = A[2][2];
  result[3][2] = 0;

  result[0][3] = -A[0][0]*A[0][3] -A[1][0]*A[1][3] -A[2][0]*A[2][3];
  result[1][3] = -A[0][1]*A[0][3] -A[1][1]*A[1][3] -A[2][1]*A[2][3];
  result[2][3] = -A[0][2]*A[0][3] -A[1][2]*A[1][3] -A[2][2]*A[2][3];
  result[3][3] = 1;
}


inline void vmMultiply (VmMatrix4x4& result, const VmMatrix4x4& A, const VmMatrix4x4& B)
{
  for (int i=0; i<3; i++)
    for (int j=0; j<4; j++)
      result[i][j] = A[i][0]*B[0][j] + A[i][1]*B[1][j] + A[i][2]*B[2][j] + A[i][3]*B[3][j];  
  result[3][0] = result[3][1] = result[3][2] = 0;
  result[3][3] = 1;  
}




inline void vmAdd (VmMatrix2x2& result, const VmMatrix2x2& A, const VmMatrix2x2& B)
{
    result[0][0] = A[0][0] + B[0][0];
    result[0][1] = A[0][1] + B[0][1];
    result[1][0] = A[1][0] + B[1][0];
    result[1][1] = A[1][1] + B[1][1];
}


inline void vmAdd (VmVector3& result, VmVector3& v)
{
    result[0] += v[0];
    result[1] += v[1];
    result[2] += v[2];
}



inline void vmAdd (VmVector2& result, VmVector2& v)
{
    result[0] += v[0];
    result[1] += v[1];
}



inline void vmAdd (VmMatrix3x3& result, const VmMatrix3x3& A, const VmMatrix3x3& B)
{
    result[0][0] = A[0][0] + B[0][0];
    result[0][1] = A[0][1] + B[0][1];
    result[0][2] = A[0][2] + B[0][2];
    result[1][0] = A[1][0] + B[1][0];
    result[1][1] = A[1][1] + B[1][1];
    result[1][2] = A[1][2] + B[1][2];
    result[2][0] = A[2][0] + B[2][0];
    result[2][1] = A[2][1] + B[2][1];
    result[2][2] = A[2][2] + B[2][2];
}



inline void vmAdd (VmMatrix2x2& result, const VmMatrix2x2& A)
{
    result[0][0] += A[0][0];
    result[0][1] += A[0][1];
    result[1][0] += A[1][0];
    result[1][1] += A[1][1];
}




inline void vmAdd (VmMatrix3x3& result, const VmMatrix3x3& A)
{
    result[0][0] += A[0][0];
    result[0][1] += A[0][1];
    result[0][2] += A[0][2];
    result[1][0] += A[1][0];
    result[1][1] += A[1][1];
    result[1][2] += A[1][2];
    result[2][0] += A[2][0];
    result[2][1] += A[2][1];
    result[2][2] += A[2][2];
}



inline void vmAdd (VmVector2& result, const VmVector2& a, const VmVector2& b)
{
    result[0] = a[0] + b[0];
    result[1] = a[1] + b[1];
}



inline void vmAdd (VmVector3& result, const VmVector3& a, const VmVector3& b)
{
    result[0] = a[0] + b[0];
    result[1] = a[1] + b[1];
    result[2] = a[2] + b[2];
}



inline void vmScale (VmMatrix2x2& result, VmReal scale)
{
    result[0][0] *= scale;
    result[0][1] *= scale;
    result[1][0] *= scale;
    result[1][1] *= scale;
}



inline void vmScale (VmMatrix3x3& result, VmReal scale)
{
    result[0][0] *= scale;
    result[0][1] *= scale;
    result[0][2] *= scale;
    result[1][0] *= scale;
    result[1][1] *= scale;
    result[1][2] *= scale;
    result[2][0] *= scale;
    result[2][1] *= scale;
    result[2][2] *= scale;
}


inline void vmScale (VmVector2& result, VmReal scale)
{
    result[0] *= scale;
    result[1] *= scale;
}



inline void vmScale (VmVector3& result, VmReal scale)
{
    result[0] *= scale;
    result[1] *= scale;
    result[2] *= scale;
}



inline void vmScale (VmMatrix2x2& result, const VmMatrix2x2& A, VmReal scale)
{
    result[0][0] = A[0][0] * scale;
    result[0][1] = A[0][1] * scale;
    result[1][0] = A[1][0] * scale;
    result[1][1] = A[1][1] * scale;
}



inline void vmScale (VmMatrix3x3& result, const VmMatrix3x3& A, VmReal scale)
{
    result[0][0] = A[0][0] * scale;
    result[0][1] = A[0][1] * scale;
    result[0][2] = A[0][2] * scale;
    result[1][0] = A[1][0] * scale;
    result[1][1] = A[1][1] * scale;
    result[1][2] = A[1][2] * scale;
    result[2][0] = A[2][0] * scale;
    result[2][1] = A[2][1] * scale;
    result[2][2] = A[2][2] * scale;
}



inline void vmScale (VmVector2& result, const VmVector2& v, VmReal scale)
{
    result[0] = v[0] * scale;
    result[1] = v[1] * scale;
}



inline void vmScale (VmVector3& result, const VmVector3& v, VmReal scale)
{
    result[0] = v[0] * scale;
    result[1] = v[1] * scale;
    result[2] = v[2] * scale;
}



inline void vmAddScale (VmMatrix2x2& result, const VmMatrix2x2& A, VmReal  scale)
{
    result[0][0] += A[0][0]*scale;
    result[0][1] += A[0][1]*scale;
    result[1][0] += A[1][0]*scale;
    result[1][1] += A[1][1]*scale;
}




inline void vmAddScale (VmMatrix3x3& result, const VmMatrix3x3& A, VmReal  scale)
{
    result[0][0] += A[0][0]*scale;
    result[0][1] += A[0][1]*scale;
    result[0][2] += A[0][2]*scale;
    result[1][0] += A[1][0]*scale;
    result[1][1] += A[1][1]*scale;
    result[1][2] += A[1][2]*scale;
    result[2][0] += A[2][0]*scale;
    result[2][1] += A[2][1]*scale;
    result[2][2] += A[2][2]*scale;
}



inline void vmAddScale (VmVector2& result, const VmVector2& v, VmReal scale)
{
  result[0] += v[0]*scale;
  result[1] += v[1]*scale;  
}


inline void vmAddScale (VmVector3& result, const VmVector3& v, VmReal scale)
{
  result[0] += v[0]*scale;
  result[1] += v[1]*scale;  
  result[2] += v[2]*scale;  
}


inline void vmSub (VmMatrix2x2& result, const VmMatrix2x2& A)
{
    result[0][0] -= A[0][0];
    result[0][1] -= A[0][1];
    result[1][0] -= A[1][0];
    result[1][1] -= A[1][1];
}



inline void vmSub (VmMatrix3x3& result, const VmMatrix3x3& A)
{
    result[0][0] -= A[0][0];
    result[0][1] -= A[0][1];
    result[0][2] -= A[0][2];
    result[1][0] -= A[1][0];
    result[1][1] -= A[1][1];
    result[1][2] -= A[1][2];
    result[2][0] -= A[2][0];
    result[2][1] -= A[2][1];
    result[2][2] -= A[2][2];
}



inline void vmSub (VmMatrix2x2& result, const VmMatrix2x2& A, const VmMatrix2x2& B)
{
    result[0][0] = A[0][0] - B[0][0];
    result[0][1] = A[0][1] - B[0][1];
    result[1][0] = A[1][0] - B[1][0];
    result[1][1] = A[1][1] - B[1][1];
}
    


inline void vmSub (VmMatrix3x3& result, const VmMatrix3x3& A, const VmMatrix3x3& B)
{
    result[0][0] = A[0][0] - B[0][0];
    result[0][1] = A[0][1] - B[0][1];
    result[0][2] = A[0][2] - B[0][2];
    result[1][0] = A[1][0] - B[1][0];
    result[1][1] = A[1][1] - B[1][1];
    result[1][2] = A[1][2] - B[1][2];
    result[2][0] = A[2][0] - B[2][0];
    result[2][1] = A[2][1] - B[2][1];
    result[2][2] = A[2][2] - B[2][2];
}



inline void vmSub (VmVector2& result, const VmVector2& a, const VmVector2& b)
{
    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];
}



inline void vmSub (VmVector3& result, const VmVector3& a, const VmVector3& b)
{
    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];
    result[2] = a[2] - b[2];
}



inline void vmMultiply (VmMatrix2x2& result, const VmMatrix2x2& A, const VmMatrix2x2& B)
{
    // result[i][j] = A[i][k]*B[k][j]
    result[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0];
    result[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1];
    result[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0];
    result[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1];
}



inline void vmMultiply (VmMatrix2x3& result, const VmMatrix2x2& A, const VmMatrix2x3& B)
{
    // result[i][j] = A[i][k]*B[k][j]
    result[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0];
    result[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1];
    result[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2];
    result[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0];
    result[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1];
    result[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2];
}



inline void vmMultiply (VmMatrix2x2& result, const VmMatrix2x3& A, const VmMatrix3x2& B)
{
    // result[i][j] = A[i][k]*B[k][j]
    result[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
    result[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1];
    result[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0];
    result[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1];
}



inline void vmMultiply (VmMatrix2x3& result, const VmMatrix2x3& A, const VmMatrix3x3& B)
{
    // result[i][j] = A[i][k]*B[k][j]
    result[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
    result[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1];
    result[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2] + A[0][2]*B[2][2];
    result[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0];
    result[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1];
    result[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2] + A[1][2]*B[2][2];
}



inline void vmMultiply (VmMatrix3x2& result, const VmMatrix3x2& A, const VmMatrix2x2& B)
{
    // result[i][j] = A[i][k]*B[k][j]
    result[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0];
    result[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1];
    result[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0];
    result[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1];
    result[2][0] = A[2][0]*B[0][0] + A[2][1]*B[1][0];
    result[2][1] = A[2][0]*B[0][1] + A[2][1]*B[1][1];
}



inline void vmMultiply (VmMatrix3x3& result, const VmMatrix3x2& A, const VmMatrix2x3& B)
{
    // result[i][j] = A[i][k]*B[k][j]
    result[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0];
    result[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1];
    result[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2];
    result[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0];
    result[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1];
    result[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2];
    result[2][0] = A[2][0]*B[0][0] + A[2][1]*B[1][0];
    result[2][1] = A[2][0]*B[0][1] + A[2][1]*B[1][1];
    result[2][2] = A[2][0]*B[0][2] + A[2][1]*B[1][2];
}



inline void vmMultiply (VmMatrix3x2& result, const VmMatrix3x3& A, const VmMatrix3x2& B)
{
    // result[i][j] = A[i][k]*B[k][j]
    result[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
    result[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1];
    result[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0];
    result[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1];
    result[2][0] = A[2][0]*B[0][0] + A[2][1]*B[1][0] + A[2][2]*B[2][0];
    result[2][1] = A[2][0]*B[0][1] + A[2][1]*B[1][1] + A[2][2]*B[2][1];
}



inline void vmMultiply (VmMatrix3x3& result, const VmMatrix3x3& A, const VmMatrix3x3& B)
{
    // result[i][j] = A[i][k]*B[k][j]
    result[0][0] = A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
    result[0][1] = A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1];
    result[0][2] = A[0][0]*B[0][2] + A[0][1]*B[1][2] + A[0][2]*B[2][2];
    result[1][0] = A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0];
    result[1][1] = A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1];
    result[1][2] = A[1][0]*B[0][2] + A[1][1]*B[1][2] + A[1][2]*B[2][2];
    result[2][0] = A[2][0]*B[0][0] + A[2][1]*B[1][0] + A[2][2]*B[2][0];
    result[2][1] = A[2][0]*B[0][1] + A[2][1]*B[1][1] + A[2][2]*B[2][1];
    result[2][2] = A[2][0]*B[0][2] + A[2][1]*B[1][2] + A[2][2]*B[2][2];
}



inline void vmMultiplyAdd (VmMatrix2x2& result, const VmMatrix2x2& A, const VmMatrix2x2& B)
{
    // result[i][j] += A[i][k]*B[k][j]
    result[0][0] += A[0][0]*B[0][0] + A[0][1]*B[1][0];
    result[0][1] += A[0][0]*B[0][1] + A[0][1]*B[1][1];
    result[1][0] += A[1][0]*B[0][0] + A[1][1]*B[1][0];
    result[1][1] += A[1][0]*B[0][1] + A[1][1]*B[1][1];
}



inline void vmMultiplyAdd (VmMatrix3x3& result, const VmMatrix3x3& A, const VmMatrix3x3& B)
{
    // result[i][j] += A[i][k]*B[k][j]
    result[0][0] += A[0][0]*B[0][0] + A[0][1]*B[1][0] + A[0][2]*B[2][0];
    result[0][1] += A[0][0]*B[0][1] + A[0][1]*B[1][1] + A[0][2]*B[2][1];
    result[0][2] += A[0][0]*B[0][2] + A[0][1]*B[1][2] + A[0][2]*B[2][2];
    result[1][0] += A[1][0]*B[0][0] + A[1][1]*B[1][0] + A[1][2]*B[2][0];
    result[1][1] += A[1][0]*B[0][1] + A[1][1]*B[1][1] + A[1][2]*B[2][1];
    result[1][2] += A[1][0]*B[0][2] + A[1][1]*B[1][2] + A[1][2]*B[2][2];
    result[2][0] += A[2][0]*B[0][0] + A[2][1]*B[1][0] + A[2][2]*B[2][0];
    result[2][1] += A[2][0]*B[0][1] + A[2][1]*B[1][1] + A[2][2]*B[2][1];
    result[2][2] += A[2][0]*B[0][2] + A[2][1]*B[1][2] + A[2][2]*B[2][2];
}



inline void vmSymmetrize (VmMatrix2x2& A)
{
    A[0][1] = A[1][0] = (A[0][1] + A[1][0])/2;
}



inline void vmSymmetrize (VmMatrix3x3& A)
{
    A[0][1] = A[1][0] = (A[0][1] + A[1][0])/2;
    A[0][2] = A[2][0] = (A[0][2] + A[2][0])/2;
    A[1][2] = A[2][1] = (A[1][2] + A[2][1])/2;
}



inline void vmJCKt (double&    result, const   VmVector2& J, const VmMatrix2x2& C, const   VmVector2& K)
{
    // result[i][j] = J[i][k]*C[k][l]*K[j][l]
    // result[i][j] = J[i][0]*C[0][0] + J[i][1]*C[1][0])*K[j][0] + (J[i][0]*C[0][1] + J[i][1]*C[1][1])*K[j][1];
    result = (J[0]*C[0][0] + J[1]*C[1][0])*K[0] + (J[0]*C[0][1] + J[1]*C[1][1])*K[1];
}



inline void vmJCKt (double&    result, const   VmVector3& J, const VmMatrix3x3& C, const   VmVector3& K)
{
    // result[i][j] = J[i][k]*C[k][l]*J[j][l]
    // result[i][j] = (J[i][0]*C[0][0] + J[i][1]*C[1][0] + J[i][2]*C[2][0])*K[j][0] +
    //                (J[i][0]*C[0][1] + J[i][1]*C[1][1] + J[i][2]*C[2][1])*K[j][1] +
    //                (J[i][0]*C[0][2] + J[i][1]*C[1][2] + J[i][2]*C[2][2])*K[j][2];
    result = (J[0]*C[0][0] + J[1]*C[1][0] + J[2]*C[2][0])*K[0] +
        (J[0]*C[0][1] + J[1]*C[1][1] + J[2]*C[2][1])*K[1] +
        (J[0]*C[0][2] + J[1]*C[1][2] + J[2]*C[2][2])*K[2]; 
}




inline void vmJCKt (VmMatrix2x2& result, const   VmVector2& J, double C, const   VmVector2& K)
{
    // result[i][j] = J[i][k]*C[k][l]*K[j][l]
    // result[i][j] = J[i]*C*K[j];
    result[0][0] = J[0]*C*K[0];
    result[0][1] = J[0]*C*K[1];
    result[1][0] = J[1]*C*K[0];
    result[1][1] = J[1]*C*K[1];    
}



inline void vmJCKt (VmMatrix2x2& result, const VmMatrix2x2& J, const VmMatrix2x2& C, const VmMatrix2x2& K)
{
    // result[i][j] = J[i][k]*C[k][l]*K[j][l]
    // result[i][j] = (J[i][0]*C[0][0] + J[i][1]*C[1][0])*K[j][0] + (J[i][0]*C[0][1] + J[i][1]*C[1][1])*K[j][1];
    result[0][0] = (J[0][0]*C[0][0] + J[0][1]*C[1][0])*K[0][0] + (J[0][0]*C[0][1] + J[0][1]*C[1][1])*K[0][1];
    result[0][1] = (J[0][0]*C[0][0] + J[0][1]*C[1][0])*K[1][0] + (J[0][0]*C[0][1] + J[0][1]*C[1][1])*K[1][1];
    result[1][0] = (J[1][0]*C[0][0] + J[1][1]*C[1][0])*K[0][0] + (J[1][0]*C[0][1] + J[1][1]*C[1][1])*K[0][1];
    result[1][1] = (J[1][0]*C[0][0] + J[1][1]*C[1][0])*K[1][0] + (J[1][0]*C[0][1] + J[1][1]*C[1][1])*K[1][1];
}



inline void vmJCKt (VmMatrix2x2& result, const VmMatrix2x3& J, const VmMatrix3x3& C, const VmMatrix2x3& K)
{
    // result[i][j] = J[i][k]*C[k][l]*K[j][l]
    // result[i][j] = (J[i][0]*C[0][0] + J[i][1]*C[1][0] + J[i][2]*C[2][0])*K[j][0] +
    //                (J[i][0]*C[0][1] + J[i][1]*C[1][1] + J[i][2]*C[2][1])*K[j][1] +
    //                (J[i][0]*C[0][2] + J[i][1]*C[1][2] + J[i][2]*C[2][2])*K[j][2];
    result[0][0] = 
        (J[0][0]*C[0][0] + J[0][1]*C[1][0] + J[0][2]*C[2][0])*K[0][0] +
        (J[0][0]*C[0][1] + J[0][1]*C[1][1] + J[0][2]*C[2][1])*K[0][1] +
        (J[0][0]*C[0][2] + J[0][1]*C[1][2] + J[0][2]*C[2][2])*K[0][2];
    result[0][1] = 
        (J[0][0]*C[0][0] + J[0][1]*C[1][0] + J[0][2]*C[2][0])*K[1][0] +
        (J[0][0]*C[0][1] + J[0][1]*C[1][1] + J[0][2]*C[2][1])*K[1][1] +
        (J[0][0]*C[0][2] + J[0][1]*C[1][2] + J[0][2]*C[2][2])*K[1][2];
    result[1][0] = 
        (J[1][0]*C[0][0] + J[1][1]*C[1][0] + J[1][2]*C[2][0])*K[0][0] +
        (J[1][0]*C[0][1] + J[1][1]*C[1][1] + J[1][2]*C[2][1])*K[0][1] +
        (J[1][0]*C[0][2] + J[1][1]*C[1][2] + J[1][2]*C[2][2])*K[0][2];
    result[1][1] = 
        (J[1][0]*C[0][0] + J[1][1]*C[1][0] + J[1][2]*C[2][0])*K[1][0] +
        (J[1][0]*C[0][1] + J[1][1]*C[1][1] + J[1][2]*C[2][1])*K[1][1] +
        (J[1][0]*C[0][2] + J[1][1]*C[1][2] + J[1][2]*C[2][2])*K[1][2];
}



inline void vmJCKt (VmMatrix3x3& result, const   VmVector3& J, double C, const   VmVector3& K)
{
    // result[i][j] = J[i][k]*C[k][l]*K[j][l]
    // result[i][j] = J[i]*C*K[j];
    result[0][0] = J[0]*C*K[0];
    result[0][1] = J[0]*C*K[1];
    result[0][2] = J[0]*C*K[2];
    result[1][0] = J[1]*C*K[0];
    result[1][1] = J[1]*C*K[1];
    result[1][2] = J[1]*C*K[2];
    result[2][0] = J[2]*C*K[0];
    result[2][1] = J[2]*C*K[1];
    result[2][2] = J[2]*C*K[2];
}



inline void vmJCKt (VmMatrix3x3& result, const VmMatrix3x2& J, const VmMatrix2x2& C, const VmMatrix3x2& K)
{
    // result[i][j] = J[i][k]*C[k][l]*K[j][l]
    // result[i][j] = (J[i][0]*C[0][0] + J[i][1]*C[1][0])*J[j][0] + (J[i][0]*C[0][1] + J[i][1]*C[1][1])*K[j][1];
    result[0][0] = (J[0][0]*C[0][0] + J[0][1]*C[1][0])*K[0][0] + (J[0][0]*C[0][1] + J[0][1]*C[1][1])*K[0][1];
    result[0][1] = (J[0][0]*C[0][0] + J[0][1]*C[1][0])*K[1][0] + (J[0][0]*C[0][1] + J[0][1]*C[1][1])*K[1][1];
    result[0][2] = (J[0][0]*C[0][0] + J[0][1]*C[1][0])*K[2][0] + (J[0][0]*C[0][1] + J[0][1]*C[1][1])*K[2][1];
    result[1][0] = (J[1][0]*C[0][0] + J[1][1]*C[1][0])*K[0][0] + (J[1][0]*C[0][1] + J[1][1]*C[1][1])*K[0][1];
    result[1][1] = (J[1][0]*C[0][0] + J[1][1]*C[1][0])*K[1][0] + (J[1][0]*C[0][1] + J[1][1]*C[1][1])*K[1][1];
    result[1][2] = (J[1][0]*C[0][0] + J[1][1]*C[1][0])*K[2][0] + (J[1][0]*C[0][1] + J[1][1]*C[1][1])*K[2][1];
    result[2][0] = (J[2][0]*C[0][0] + J[2][1]*C[1][0])*K[0][0] + (J[2][0]*C[0][1] + J[2][1]*C[1][1])*K[0][1];
    result[2][1] = (J[2][0]*C[0][0] + J[2][1]*C[1][0])*K[1][0] + (J[2][0]*C[0][1] + J[2][1]*C[1][1])*K[1][1];
    result[2][2] = (J[2][0]*C[0][0] + J[2][1]*C[1][0])*K[2][0] + (J[2][0]*C[0][1] + J[2][1]*C[1][1])*K[2][1];
}



inline void vmJCKt (VmMatrix2x3& result, const VmMatrix2x2& J, const VmMatrix2x3& C, const VmMatrix3x3& K)
    // result[i][j] = J[i][k]*C[k][l]*K[j][l]
    // result[i][j] = (J[i][0]*C[0][0] + J[i][1]*C[1][0] + J[i][2]*C[2][0])*K[j][0] +
    //                (J[i][0]*C[0][1] + J[i][1]*C[1][1] + J[i][2]*C[2][1])*K[j][1] +
    //                (J[i][0]*C[0][2] + J[i][1]*C[1][2] + J[i][2]*C[2][2])*K[j][2];
{
    result[0][0] = 
        (J[0][0]*C[0][0] + J[0][1]*C[1][0] + J[0][2]*C[2][0])*K[0][0] +
        (J[0][0]*C[0][1] + J[0][1]*C[1][1] + J[0][2]*C[2][1])*K[0][1] +
        (J[0][0]*C[0][2] + J[0][1]*C[1][2] + J[0][2]*C[2][2])*K[0][2];
    result[0][1] = 
        (J[0][0]*C[0][0] + J[0][1]*C[1][0] + J[0][2]*C[2][0])*K[1][0] +
        (J[0][0]*C[0][1] + J[0][1]*C[1][1] + J[0][2]*C[2][1])*K[1][1] +
        (J[0][0]*C[0][2] + J[0][1]*C[1][2] + J[0][2]*C[2][2])*K[1][2];
    result[0][2] = 
        (J[0][0]*C[0][0] + J[0][1]*C[1][0] + J[0][2]*C[2][0])*K[2][0] +
        (J[0][0]*C[0][1] + J[0][1]*C[1][1] + J[0][2]*C[2][1])*K[2][1] +
        (J[0][0]*C[0][2] + J[0][1]*C[1][2] + J[0][2]*C[2][2])*K[2][2];
    result[1][0] = 
        (J[1][0]*C[0][0] + J[1][1]*C[1][0] + J[1][2]*C[2][0])*K[0][0] +
        (J[1][0]*C[0][1] + J[1][1]*C[1][1] + J[1][2]*C[2][1])*K[0][1] +
        (J[1][0]*C[0][2] + J[1][1]*C[1][2] + J[1][2]*C[2][2])*K[0][2];
    result[1][1] = 
        (J[1][0]*C[0][0] + J[1][1]*C[1][0] + J[1][2]*C[2][0])*K[1][0] +
        (J[1][0]*C[0][1] + J[1][1]*C[1][1] + J[1][2]*C[2][1])*K[1][1] +
        (J[1][0]*C[0][2] + J[1][1]*C[1][2] + J[1][2]*C[2][2])*K[1][2];
    result[1][2] = 
        (J[1][0]*C[0][0] + J[1][1]*C[1][0] + J[1][2]*C[2][0])*K[2][0] +
        (J[1][0]*C[0][1] + J[1][1]*C[1][1] + J[1][2]*C[2][1])*K[2][1] +
        (J[1][0]*C[0][2] + J[1][1]*C[1][2] + J[1][2]*C[2][2])*K[2][2];
}



inline void vmJCKt (VmMatrix3x3& result, const VmMatrix3x3& J, const VmMatrix3x3& C, const VmMatrix3x3& K)
{
    // result[i][j] = J[i][k]*C[k][l]*K[j][l]
    // result[i][j] = (J[i][0]*C[0][0] + J[i][1]*C[1][0] + J[i][2]*C[2][0])*K[j][0] +
    //                (J[i][0]*C[0][1] + J[i][1]*C[1][1] + J[i][2]*C[2][1])*K[j][1] +
    //                (J[i][0]*C[0][2] + J[i][1]*C[1][2] + J[i][2]*C[2][2])*K[j][2];
    result[0][0] = 
        (J[0][0]*C[0][0] + J[0][1]*C[1][0] + J[0][2]*C[2][0])*K[0][0] +
        (J[0][0]*C[0][1] + J[0][1]*C[1][1] + J[0][2]*C[2][1])*K[0][1] +
        (J[0][0]*C[0][2] + J[0][1]*C[1][2] + J[0][2]*C[2][2])*K[0][2];
    result[0][1] = 
        (J[0][0]*C[0][0] + J[0][1]*C[1][0] + J[0][2]*C[2][0])*K[1][0] +
        (J[0][0]*C[0][1] + J[0][1]*C[1][1] + J[0][2]*C[2][1])*K[1][1] +
        (J[0][0]*C[0][2] + J[0][1]*C[1][2] + J[0][2]*C[2][2])*K[1][2];
    result[0][2] = 
        (J[0][0]*C[0][0] + J[0][1]*C[1][0] + J[0][2]*C[2][0])*K[2][0] +
        (J[0][0]*C[0][1] + J[0][1]*C[1][1] + J[0][2]*C[2][1])*K[2][1] +
        (J[0][0]*C[0][2] + J[0][1]*C[1][2] + J[0][2]*C[2][2])*K[2][2];
    result[1][0] = 
        (J[1][0]*C[0][0] + J[1][1]*C[1][0] + J[1][2]*C[2][0])*K[0][0] +
        (J[1][0]*C[0][1] + J[1][1]*C[1][1] + J[1][2]*C[2][1])*K[0][1] +
        (J[1][0]*C[0][2] + J[1][1]*C[1][2] + J[1][2]*C[2][2])*K[0][2];
    result[1][1] = 
        (J[1][0]*C[0][0] + J[1][1]*C[1][0] + J[1][2]*C[2][0])*K[1][0] +
        (J[1][0]*C[0][1] + J[1][1]*C[1][1] + J[1][2]*C[2][1])*K[1][1] +
        (J[1][0]*C[0][2] + J[1][1]*C[1][2] + J[1][2]*C[2][2])*K[1][2];
    result[1][2] = 
        (J[1][0]*C[0][0] + J[1][1]*C[1][0] + J[1][2]*C[2][0])*K[2][0] +
        (J[1][0]*C[0][1] + J[1][1]*C[1][1] + J[1][2]*C[2][1])*K[2][1] +
        (J[1][0]*C[0][2] + J[1][1]*C[1][2] + J[1][2]*C[2][2])*K[2][2];
    result[2][0] = 
        (J[2][0]*C[0][0] + J[2][1]*C[1][0] + J[2][2]*C[2][0])*K[0][0] +
        (J[2][0]*C[0][1] + J[2][1]*C[1][1] + J[2][2]*C[2][1])*K[0][1] +
        (J[2][0]*C[0][2] + J[2][1]*C[1][2] + J[2][2]*C[2][2])*K[0][2];
    result[2][1] = 
        (J[2][0]*C[0][0] + J[2][1]*C[1][0] + J[2][2]*C[2][0])*K[1][0] +
        (J[2][0]*C[0][1] + J[2][1]*C[1][1] + J[2][2]*C[2][1])*K[1][1] +
        (J[2][0]*C[0][2] + J[2][1]*C[1][2] + J[2][2]*C[2][2])*K[1][2];
    result[2][2] = 
        (J[2][0]*C[0][0] + J[2][1]*C[1][0] + J[2][2]*C[2][0])*K[2][0] +
        (J[2][0]*C[0][1] + J[2][1]*C[1][1] + J[2][2]*C[2][1])*K[2][1] +
        (J[2][0]*C[0][2] + J[2][1]*C[1][2] + J[2][2]*C[2][2])*K[2][2];
}





inline void vmJCKtAdd (double&    result, const   VmVector2& J, const VmMatrix2x2& C, const   VmVector2& K)
{
    // result[i][j] += J[i][k]*C[k][l]*K[j][l]
    // result[i][j] += J[i][0]*C[0][0] + J[i][1]*C[1][0])*K[j][0] + (J[i][0]*C[0][1] + J[i][1]*C[1][1])*K[j][1];
    result += (J[0]*C[0][0] + J[1]*C[1][0])*K[0] + (J[0]*C[0][1] + J[1]*C[1][1])*K[1];
}




inline void vmJCKtAdd (double&    result, const   VmVector3& J, const VmMatrix3x3& C, const   VmVector3& K)
{
    // result[i][j] += J[i][k]*C[k][l]*K[j][l]
    // result[i][j] += (J[i][0]*C[0][0] + J[i][1]*C[1][0] + J[i][2]*C[2][0])*K[j][0] +
    //                (J[i][0]*C[0][1] + J[i][1]*C[1][1] + J[i][2]*C[2][1])*K[j][1] +
    //                (J[i][0]*C[0][2] + J[i][1]*C[1][2] + J[i][2]*C[2][2])*K[j][2];
    result += 
        (J[0]*C[0][0] + J[1]*C[1][0] + J[2]*C[2][0])*K[0] +
        (J[0]*C[0][1] + J[1]*C[1][1] + J[2]*C[2][1])*K[1] +
        (J[0]*C[0][2] + J[1]*C[1][2] + J[2]*C[2][2])*K[2]; 
}




inline void vmJCKtAdd (VmMatrix2x2& result, const   VmVector2& J, double C, const   VmVector2& K)
{
    // result[i][j] += J[i][k]*C[k][l]*K[j][l]
    // result[i][j] += J[i]*C*K[j];
    result[0][0] += J[0]*C*K[0];
    result[0][1] += J[0]*C*K[1];
    result[1][0] += J[1]*C*K[0];
    result[1][1] += J[1]*C*K[1];    
}



inline void vmJCKtAdd (VmMatrix2x2& result, const VmMatrix2x2& J, const VmMatrix2x2& C, const VmMatrix2x2& K)
{
    // result[i][j] += J[i][k]*C[k][l]*K[j][l]
    // result[i][j] += (J[i][0]*C[0][0] + J[i][1]*C[1][0])*K[j][0] + (J[i][0]*C[0][1] + J[i][1]*C[1][1])*K[j][1];
    result[0][0] += (J[0][0]*C[0][0] + J[0][1]*C[1][0])*K[0][0] + (J[0][0]*C[0][1] + J[0][1]*C[1][1])*K[0][1];
    result[0][1] += (J[0][0]*C[0][0] + J[0][1]*C[1][0])*K[1][0] + (J[0][0]*C[0][1] + J[0][1]*C[1][1])*K[1][1];
    result[1][0] += (J[1][0]*C[0][0] + J[1][1]*C[1][0])*K[0][0] + (J[1][0]*C[0][1] + J[1][1]*C[1][1])*K[0][1];
    result[1][1] += (J[1][0]*C[0][0] + J[1][1]*C[1][0])*K[1][0] + (J[1][0]*C[0][1] + J[1][1]*C[1][1])*K[1][1];
}
 


inline void vmJCKtAdd (VmMatrix2x2& result, const VmMatrix2x3& J, const VmMatrix3x3& C, const VmMatrix2x3& K)
{
    // result[i][j] += J[i][k]*C[k][l]*K[j][l]
    // result[i][j] += (J[i][0]*C[0][0] + J[i][1]*C[1][0] + J[i][2]*C[2][0])*K[j][0] +
    //                (J[i][0]*C[0][1] + J[i][1]*C[1][1] + J[i][2]*C[2][1])*K[j][1] +
    //                (J[i][0]*C[0][2] + J[i][1]*C[1][2] + J[i][2]*C[2][2])*K[j][2];
    result[0][0] += 
        (J[0][0]*C[0][0] + J[0][1]*C[1][0] + J[0][2]*C[2][0])*K[0][0] +
        (J[0][0]*C[0][1] + J[0][1]*C[1][1] + J[0][2]*C[2][1])*K[0][1] +
        (J[0][0]*C[0][2] + J[0][1]*C[1][2] + J[0][2]*C[2][2])*K[0][2];
    result[0][1] += 
        (J[0][0]*C[0][0] + J[0][1]*C[1][0] + J[0][2]*C[2][0])*K[1][0] +
        (J[0][0]*C[0][1] + J[0][1]*C[1][1] + J[0][2]*C[2][1])*K[1][1] +
        (J[0][0]*C[0][2] + J[0][1]*C[1][2] + J[0][2]*C[2][2])*K[1][2];
    result[1][0] += 
        (J[1][0]*C[0][0] + J[1][1]*C[1][0] + J[1][2]*C[2][0])*K[0][0] +
        (J[1][0]*C[0][1] + J[1][1]*C[1][1] + J[1][2]*C[2][1])*K[0][1] +
        (J[1][0]*C[0][2] + J[1][1]*C[1][2] + J[1][2]*C[2][2])*K[0][2];
    result[1][1] += 
        (J[1][0]*C[0][0] + J[1][1]*C[1][0] + J[1][2]*C[2][0])*K[1][0] +
        (J[1][0]*C[0][1] + J[1][1]*C[1][1] + J[1][2]*C[2][1])*K[1][1] +
        (J[1][0]*C[0][2] + J[1][1]*C[1][2] + J[1][2]*C[2][2])*K[1][2];
}



inline void vmJCKtAdd (VmMatrix3x3& result, const   VmVector3& J, double C, const   VmVector3& K)
{
    // result[i][j] += J[i][k]*C[k][l]*K[j][l]
    // result[i][j] += J[i]*C*K[j];
    result[0][0] += J[0]*C*K[0];
    result[0][1] += J[0]*C*K[1];
    result[0][2] += J[0]*C*K[2];
    result[1][0] += J[1]*C*K[0];
    result[1][1] += J[1]*C*K[1];
    result[1][2] += J[1]*C*K[2];
    result[2][0] += J[2]*C*K[0];
    result[2][1] += J[2]*C*K[1];
    result[2][2] += J[2]*C*K[2];
}



inline void vmJCKtAdd (VmMatrix3x3& result, const VmMatrix3x2& J, const VmMatrix2x2& C, const VmMatrix3x2& K)
{
    // result[i][j] += J[i][k]*C[k][l]*K[j][l]
    // result[i][j] += (J[i][0]*C[0][0] + J[i][1]*C[1][0])*J[j][0] + (J[i][0]*C[0][1] + J[i][1]*C[1][1])*K[j][1];
    result[0][0] += (J[0][0]*C[0][0] + J[0][1]*C[1][0])*K[0][0] + (J[0][0]*C[0][1] + J[0][1]*C[1][1])*K[0][1];
    result[0][1] += (J[0][0]*C[0][0] + J[0][1]*C[1][0])*K[1][0] + (J[0][0]*C[0][1] + J[0][1]*C[1][1])*K[1][1];
    result[0][2] += (J[0][0]*C[0][0] + J[0][1]*C[1][0])*K[2][0] + (J[0][0]*C[0][1] + J[0][1]*C[1][1])*K[2][1];
    result[1][0] += (J[1][0]*C[0][0] + J[1][1]*C[1][0])*K[0][0] + (J[1][0]*C[0][1] + J[1][1]*C[1][1])*K[0][1];
    result[1][1] += (J[1][0]*C[0][0] + J[1][1]*C[1][0])*K[1][0] + (J[1][0]*C[0][1] + J[1][1]*C[1][1])*K[1][1];
    result[1][2] += (J[1][0]*C[0][0] + J[1][1]*C[1][0])*K[2][0] + (J[1][0]*C[0][1] + J[1][1]*C[1][1])*K[2][1];
    result[2][0] += (J[2][0]*C[0][0] + J[2][1]*C[1][0])*K[0][0] + (J[2][0]*C[0][1] + J[2][1]*C[1][1])*K[0][1];
    result[2][1] += (J[2][0]*C[0][0] + J[2][1]*C[1][0])*K[1][0] + (J[2][0]*C[0][1] + J[2][1]*C[1][1])*K[1][1];
    result[2][2] += (J[2][0]*C[0][0] + J[2][1]*C[1][0])*K[2][0] + (J[2][0]*C[0][1] + J[2][1]*C[1][1])*K[2][1];
}



inline void vmJCKtAdd (VmMatrix2x3& result, const VmMatrix2x2& J, const VmMatrix2x3& C, const VmMatrix3x3& K)
{
    // result[i][j] = J[i][k]*C[k][l]*K[j][l]
    // result[i][j] = (J[i][0]*C[0][0] + J[i][1]*C[1][0] + J[i][2]*C[2][0])*K[j][0] +
    //                (J[i][0]*C[0][1] + J[i][1]*C[1][1] + J[i][2]*C[2][1])*K[j][1] +
    //                (J[i][0]*C[0][2] + J[i][1]*C[1][2] + J[i][2]*C[2][2])*K[j][2];
    result[0][0] += 
        (J[0][0]*C[0][0] + J[0][1]*C[1][0] + J[0][2]*C[2][0])*K[0][0] +
        (J[0][0]*C[0][1] + J[0][1]*C[1][1] + J[0][2]*C[2][1])*K[0][1] +
        (J[0][0]*C[0][2] + J[0][1]*C[1][2] + J[0][2]*C[2][2])*K[0][2];
    result[0][1] += 
        (J[0][0]*C[0][0] + J[0][1]*C[1][0] + J[0][2]*C[2][0])*K[1][0] +
        (J[0][0]*C[0][1] + J[0][1]*C[1][1] + J[0][2]*C[2][1])*K[1][1] +
        (J[0][0]*C[0][2] + J[0][1]*C[1][2] + J[0][2]*C[2][2])*K[1][2];
    result[0][2] += 
        (J[0][0]*C[0][0] + J[0][1]*C[1][0] + J[0][2]*C[2][0])*K[2][0] +
        (J[0][0]*C[0][1] + J[0][1]*C[1][1] + J[0][2]*C[2][1])*K[2][1] +
        (J[0][0]*C[0][2] + J[0][1]*C[1][2] + J[0][2]*C[2][2])*K[2][2];
    result[1][0] += 
        (J[1][0]*C[0][0] + J[1][1]*C[1][0] + J[1][2]*C[2][0])*K[0][0] +
        (J[1][0]*C[0][1] + J[1][1]*C[1][1] + J[1][2]*C[2][1])*K[0][1] +
        (J[1][0]*C[0][2] + J[1][1]*C[1][2] + J[1][2]*C[2][2])*K[0][2];
    result[1][1] += 
        (J[1][0]*C[0][0] + J[1][1]*C[1][0] + J[1][2]*C[2][0])*K[1][0] +
        (J[1][0]*C[0][1] + J[1][1]*C[1][1] + J[1][2]*C[2][1])*K[1][1] +
        (J[1][0]*C[0][2] + J[1][1]*C[1][2] + J[1][2]*C[2][2])*K[1][2];
    result[1][2] += 
        (J[1][0]*C[0][0] + J[1][1]*C[1][0] + J[1][2]*C[2][0])*K[2][0] +
        (J[1][0]*C[0][1] + J[1][1]*C[1][1] + J[1][2]*C[2][1])*K[2][1] +
        (J[1][0]*C[0][2] + J[1][1]*C[1][2] + J[1][2]*C[2][2])*K[2][2];
}



inline void vmJCKtAdd (VmMatrix3x3& result, const VmMatrix3x3& J, const VmMatrix3x3& C, const VmMatrix3x3& K)
{
    // result[i][j] = J[i][k]*C[k][l]*K[j][l]
    // result[i][j] = (J[i][0]*C[0][0] + J[i][1]*C[1][0] + J[i][2]*C[2][0])*K[j][0] +
    //                (J[i][0]*C[0][1] + J[i][1]*C[1][1] + J[i][2]*C[2][1])*K[j][1] +
    //                (J[i][0]*C[0][2] + J[i][1]*C[1][2] + J[i][2]*C[2][2])*K[j][2];
    result[0][0] += 
        (J[0][0]*C[0][0] + J[0][1]*C[1][0] + J[0][2]*C[2][0])*K[0][0] +
        (J[0][0]*C[0][1] + J[0][1]*C[1][1] + J[0][2]*C[2][1])*K[0][1] +
        (J[0][0]*C[0][2] + J[0][1]*C[1][2] + J[0][2]*C[2][2])*K[0][2];
    result[0][1] += 
        (J[0][0]*C[0][0] + J[0][1]*C[1][0] + J[0][2]*C[2][0])*K[1][0] +
        (J[0][0]*C[0][1] + J[0][1]*C[1][1] + J[0][2]*C[2][1])*K[1][1] +
        (J[0][0]*C[0][2] + J[0][1]*C[1][2] + J[0][2]*C[2][2])*K[1][2];
    result[0][2] += 
        (J[0][0]*C[0][0] + J[0][1]*C[1][0] + J[0][2]*C[2][0])*K[2][0] +
        (J[0][0]*C[0][1] + J[0][1]*C[1][1] + J[0][2]*C[2][1])*K[2][1] +
        (J[0][0]*C[0][2] + J[0][1]*C[1][2] + J[0][2]*C[2][2])*K[2][2];
    result[1][0] += 
        (J[1][0]*C[0][0] + J[1][1]*C[1][0] + J[1][2]*C[2][0])*K[0][0] +
        (J[1][0]*C[0][1] + J[1][1]*C[1][1] + J[1][2]*C[2][1])*K[0][1] +
        (J[1][0]*C[0][2] + J[1][1]*C[1][2] + J[1][2]*C[2][2])*K[0][2];
    result[1][1] += 
        (J[1][0]*C[0][0] + J[1][1]*C[1][0] + J[1][2]*C[2][0])*K[1][0] +
        (J[1][0]*C[0][1] + J[1][1]*C[1][1] + J[1][2]*C[2][1])*K[1][1] +
        (J[1][0]*C[0][2] + J[1][1]*C[1][2] + J[1][2]*C[2][2])*K[1][2];
    result[1][2] += 
        (J[1][0]*C[0][0] + J[1][1]*C[1][0] + J[1][2]*C[2][0])*K[2][0] +
        (J[1][0]*C[0][1] + J[1][1]*C[1][1] + J[1][2]*C[2][1])*K[2][1] +
        (J[1][0]*C[0][2] + J[1][1]*C[1][2] + J[1][2]*C[2][2])*K[2][2];
    result[2][0] += 
        (J[2][0]*C[0][0] + J[2][1]*C[1][0] + J[2][2]*C[2][0])*K[0][0] +
        (J[2][0]*C[0][1] + J[2][1]*C[1][1] + J[2][2]*C[2][1])*K[0][1] +
        (J[2][0]*C[0][2] + J[2][1]*C[1][2] + J[2][2]*C[2][2])*K[0][2];
    result[2][1] += 
        (J[2][0]*C[0][0] + J[2][1]*C[1][0] + J[2][2]*C[2][0])*K[1][0] +
        (J[2][0]*C[0][1] + J[2][1]*C[1][1] + J[2][2]*C[2][1])*K[1][1] +
        (J[2][0]*C[0][2] + J[2][1]*C[1][2] + J[2][2]*C[2][2])*K[1][2];
    result[2][2] += 
        (J[2][0]*C[0][0] + J[2][1]*C[1][0] + J[2][2]*C[2][0])*K[2][0] +
        (J[2][0]*C[0][1] + J[2][1]*C[1][1] + J[2][2]*C[2][1])*K[2][1] +
        (J[2][0]*C[0][2] + J[2][1]*C[1][2] + J[2][2]*C[2][2])*K[2][2];
}



inline VmReal vmDeterminant (const VmMatrix2x2& j)
{
    return j[0][0]*j[1][1]-j[1][0]*j[0][1];
}



inline VmReal vmDeterminant (const VmMatrix3x3& j)
{
    VmReal det =  j[0][0]*j[1][1]*j[2][2] + j[0][1]*j[1][2]*j[2][0] + j[0][2]*j[1][0]*j[2][1]
        -j[2][0]*j[1][1]*j[0][2] - j[2][1]*j[1][2]*j[0][0] - j[2][2]*j[1][0]*j[0][1];
    return det;
}



inline VmReal vmTrace (const VmMatrix2x2& A)
{
    return A[0][0] + A[1][1];
}



inline VmReal vmTrace (const VmMatrix3x3& A)
{
    return A[0][0] + A[1][1] + A[2][2];
}

inline VmReal vmDot (const VmVector2& v, const VmVector2& w)
{
  return v[0]*w[0] + v[1]*w[1];
}

inline VmReal vmDot (const VmVector3& v, const VmVector2& w)
{
  return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

inline VmReal vmLength2 (const VmVector2& v)
{
  return v[0]*v[0] + v[1]*v[1];
}

inline VmReal vmLength2 (const VmVector3& v)
{
  return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}


inline VmReal vmLength (const VmVector2& v)
{
  return sqrt(v[0]*v[0] + v[1]*v[1]);
}

inline VmReal vmLength (const VmVector3& v)
{
  return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}



inline VmReal vmInverse (VmMatrix2x2& result, const VmMatrix2x2& A, bool refine)
{
    VmReal det = vmDeterminant (A);
    VmReal dI  = 1/det;
    // Actual Inversion
    result[0][0] =  A[1][1] * dI;
    result[0][1] = -A[0][1] * dI;
    result[1][0] = -A[1][0] * dI;
    result[1][1] =  A[0][0] * dI;
    // Iterative improvement
    if (refine) {
        VmMatrix2x2 C;
        VmMatrix2x2 D;
        vmMultiply (C, A, result);
        vmMultiply (D, result, C);
        result[0][0] = 2*result[0][0] - D[0][0];
        result[1][0] = 2*result[1][0] - D[1][0];
        result[0][1] = 2*result[0][1] - D[0][1];
        result[1][1] = 2*result[1][1] - D[1][1];
    }
    return det;
}




inline VmReal  vmInverse (VmMatrix2x2& A, bool refine)
{
    VmMatrix2x2 inv;
    double det;
    det = vmInverse(inv, A, refine);
    vmCopy (inv, A);
    return det;
}




inline VmReal vmInverse (VmMatrix3x3& result, const VmMatrix3x3& a, bool refine)
{
    VmReal det = vmDeterminant (a);
    VmReal dI  = 1/det;
    // Actual Inversion
    result[0][0] = (-a[1][2]*a[2][1] + a[1][1]*a[2][2])*dI;
    result[0][1] = ( a[0][2]*a[2][1] - a[0][1]*a[2][2])*dI;
    result[0][2] = (-a[0][2]*a[1][1] + a[0][1]*a[1][2])*dI;
    result[1][0] = ( a[1][2]*a[2][0] - a[1][0]*a[2][2])*dI;
    result[1][1] = (-a[0][2]*a[2][0] + a[0][0]*a[2][2])*dI;
    result[1][2] = ( a[0][2]*a[1][0] - a[0][0]*a[1][2])*dI;
    result[2][0] = (-a[1][1]*a[2][0] + a[1][0]*a[2][1])*dI;
    result[2][1] = ( a[0][1]*a[2][0] - a[0][0]*a[2][1])*dI;
    result[2][2] = (-a[0][1]*a[1][0] + a[0][0]*a[1][1])*dI;
    // Iterative improvement
    if (refine) {
        VmMatrix3x3 C, D;
        vmMultiply (C, a, result);
        vmMultiply (D, result, C);
        result[0][0] = 2*result[0][0] - D[0][0];
        result[1][0] = 2*result[1][0] - D[1][0];
        result[2][0] = 2*result[2][0] - D[2][0];

        result[0][1] = 2*result[0][1] - D[0][1];
        result[1][1] = 2*result[1][1] - D[1][1];        
        result[2][1] = 2*result[2][1] - D[2][1];

        result[0][2] = 2*result[0][2] - D[0][2];
        result[1][2] = 2*result[1][2] - D[1][2];        
        result[2][2] = 2*result[2][2] - D[2][2];
    }
    return det;
}



inline VmReal vmInverse (VmMatrix3x3& A, bool refine)
{
    VmReal det;
    VmMatrix3x3 inv;
    det = vmInverse(inv, A, refine);
    vmCopy (inv, A);
    return det;
}



inline VmReal vmInverseSymmetric (VmMatrix2x2& result, const VmMatrix2x2& AOrig, bool refine)
{
    VmMatrix2x2 A;
    A[0][0] = AOrig[0][0];
    A[0][1] = A[1][0] = (AOrig[0][1]+AOrig[1][0])/2;
    A[1][1] = AOrig[1][1];

    VmReal det = vmDeterminant (A);
    VmReal dI  = 1/det;
    result[0][0] =  A[1][1] * dI;
    result[0][1] = -A[0][1] * dI;
    result[1][0] = result[0][1];
    result[1][1] =  A[0][0] * dI;    
    // Iterative improvement
    if (refine) {
        VmMatrix2x2 C, D;
        vmMultiply (C, A, result);
        vmMultiply (D, result, C);
        result[0][0] = 2*result[0][0] - D[0][0];
        result[1][0] = 2*result[1][0] - D[1][0];
        result[0][1] = result[1][0];
        result[1][1] = 2*result[1][1] - D[1][1];
    }
    return det;
}



inline VmReal  vmInverseSymmetric (VmMatrix2x2& A, bool refine)
{
    VmReal det;
    VmMatrix2x2 inv;
    det = vmInverseSymmetric(inv, A, refine);
    vmCopy (inv, A);
    return det;
}



inline VmReal vmInverseSymmetric (VmMatrix3x3& result, const VmMatrix3x3& AOrig, bool refine)
{
    VmMatrix3x3 A;
    A[0][0] = AOrig[0][0];
    A[0][1] = A[1][0] = (AOrig[0][1]+AOrig[1][0])/2;
    A[0][2] = A[2][0] = (AOrig[0][2]+AOrig[2][0])/2;
    A[1][1] = AOrig[1][1];
    A[1][2] = A[2][1] = (AOrig[1][2]+AOrig[2][1])/2;
    A[2][2] = AOrig[2][2];

    VmReal det = vmDeterminant (A);
    VmReal dI  = 1/det;
    result[0][0] = (-A[1][2]*A[2][1] + A[1][1]*A[2][2])*dI;
    result[0][1] = ( A[0][2]*A[2][1] - A[0][1]*A[2][2])*dI;
    result[0][2] = (-A[0][2]*A[1][1] + A[0][1]*A[1][2])*dI;
    result[1][0] = result[0][1];
    result[1][1] = (-A[0][2]*A[2][0] + A[0][0]*A[2][2])*dI;
    result[1][2] = ( A[0][2]*A[1][0] - A[0][0]*A[1][2])*dI;
    result[2][0] = result[0][2];
    result[2][1] = result[1][2];
    result[2][2] = (-A[0][1]*A[1][0] + A[0][0]*A[1][1])*dI;
    // Iterative improvement
    if (refine) {
        VmMatrix3x3 C, D;
        vmMultiply (C, A, result);
        vmMultiply (D, result, C);
        result[0][0] = 2*result[0][0] - D[0][0];
        result[1][0] = 2*result[1][0] - D[1][0];
        result[2][0] = 2*result[2][0] - D[2][0];

        result[0][1] = result[1][0];
        result[1][1] = 2*result[1][1] - D[1][1];        
        result[2][1] = 2*result[2][1] - D[2][1];

        result[0][2] = result[2][0];
        result[1][2] = result[2][1];
        result[2][2] = 2*result[2][2] - D[2][2];
    }
    return det;
}



inline VmReal vmInverseSymmetric (VmMatrix3x3& A, bool refine)
{
    VmReal det;
    VmMatrix3x3 inv;
    det = vmInverseSymmetric(inv, A, refine);
    vmCopy (inv, A);
    return det;
}



inline VmReal vmSolve (const VmMatrix2x2& A, VmVector2& x, const VmVector2& b)
{
    double det  = A[0][0]*A[1][1] - A[1][0]*A[0][1];
    double detX = b[0]   *A[1][1] - b[1]   *A[0][1];
    double detY = A[0][0]*b[1]    - A[1][0]*b[0];
    x[0] = detX/det;
    x[1] = detY/det;
    return det;
}



inline VmReal vmSolve (const VmMatrix3x3& A, VmVector3& x, const VmVector3& b)
{
    double det  =  A[0][0]*A[1][1]*A[2][2] + A[0][1]*A[1][2]*A[2][0] + A[0][2]*A[1][0]*A[2][1]
                  -A[2][0]*A[1][1]*A[0][2] - A[2][1]*A[1][2]*A[0][0] - A[2][2]*A[1][0]*A[0][1];
    double detX =  b[0]   *A[1][1]*A[2][2] + A[0][1]*A[1][2]*b[2]    + A[0][2]*b[1]   *A[2][1]
                  -b[2]   *A[1][1]*A[0][2] - A[2][1]*A[1][2]*b[0]    - A[2][2]*b[1]   *A[0][1];
    double detY =  A[0][0]*b[1]   *A[2][2] + b[0]   *A[1][2]*A[2][0] + A[0][2]*A[1][0]*b[2]   
                  -A[2][0]*b[1]   *A[0][2] - b[2]   *A[1][2]*A[0][0] - A[2][2]*A[1][0]*b[0];
    double detZ =  A[0][0]*A[1][1]*b[2]    + A[0][1]*b[1]   *A[2][0] + b[0]   *A[1][0]*A[2][1]
                  -A[2][0]*A[1][1]*b[0]    - A[2][1]*b[1]   *A[0][0] - b[2]   *A[1][0]*A[0][1];
    x[0] = detX/det;
    x[1] = detY/det;
    x[2] = detZ/det;
    return det;
}


inline void vmRotateSym (VmMatrix2x2& A, double cAngle, double sAngle)
{
    VmMatrix2x2 B;
    B[0][0] =  A[0][0]*cAngle + A[0][1]*sAngle;
    B[0][1] = -A[0][0]*sAngle + A[0][1]*cAngle;
    B[1][0] =  A[1][0]*cAngle + A[1][1]*sAngle;
    B[1][1] = -A[1][0]*sAngle + A[1][1]*cAngle;
    
    A[0][0] = B[0][0]*cAngle + B[1][0]*sAngle;
    A[0][1] = B[0][1]*cAngle + B[1][1]*sAngle;
    A[1][0] = B[1][0]*cAngle - B[0][0]*sAngle;
    A[1][1] = B[1][1]*cAngle - B[0][1]*sAngle;
}


inline void vmRotate (VmVector2& result, const VmVector2& v, double angle)
{
  double c = cos(angle), s = sin(angle);
  double v0 = v[0];  
  result[0] = c*v0 - s*v[1];
  result[1] = s*v0 + c*v[1];  
}
