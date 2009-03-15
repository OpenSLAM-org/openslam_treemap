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
/*! \file vectormath.h
    \brief matrix algebra for small matrices
    \author Udo Frese

    Contains matrix algebra routines for very small matrices and
    vectors of fixed size (2,3-vector, {2,3}x{2,3}-matrix). These
    matrices are all simple \c double[n][m] types and do not need
    any dynamic memory allocation
*/


#ifndef VECTORMATH_H
#define VECTORMATH_H


#include "vmMathMisc.h"

//! Small 2x2 matrix
typedef VmReal VmMatrix2x2[2][2];
//! Small 2x3 matrix
typedef VmReal VmMatrix2x3[2][3];
//! Small 3x2 matrix
typedef VmReal VmMatrix3x2[3][2];
//! Small 3x3 matrix
typedef VmReal VmMatrix3x3[3][3];
//! Small 4x4 matrix
typedef VmReal VmMatrix4x4[4][4];

//! Small 2 vector
typedef VmReal VmVector2  [2];
//! Small 3 vector
typedef VmReal VmVector3  [3];

//! Class wrapper for VmVector2
/*! This class is used to allow VmVector2 objects
    in an STL vector. 
 */
class VmVector2C 
{
 public:
  //! The vector data
  VmVector2 v;

  //! Zero vector
  VmVector2C (){v[0]=v[1]=0;}
  //! Constructor with given entries
  VmVector2C (VmReal v0, VmReal v1) {v[0]=v0; v[1]=v1;}  
  //! Copy constructor
  VmVector2C (const VmVector2& v2){v[0]=v2[0]; v[1]=v2[1]; }
  
  //! Entry access operator
  VmReal& operator [] (int i) {return v[i];}
  //! Entry access operator
  VmReal  operator [] (int i) const {return v[i];}
};

//! Class wrapper for VmVector3
/*! This class is used to allow VmVector3 objects
    in an STL vector. 
 */
class VmVector3C 
{
 public:
  //! The vector data
  VmVector3 v;

  //! Zero vector
  VmVector3C (){v[0]=v[1]=v[2]=0;}
  //! Constructor with given entries
  VmVector3C (VmReal v0, VmReal v1, VmReal v2) {v[0]=v0; v[1]=v1; v[2]=v2;}  
  //! Copy constructor
  VmVector3C (const VmVector3& v2) {v[0]=v2[0]; v[1]=v2[1]; v[2]=v2[2];}

  //! Entry access operator  
  VmReal& operator [] (int i) {return v[i];}
  //! Entry access operator
  VmReal  operator [] (int i) const {return v[i];}
};


//! \name Functions operating on vectors 
//@{  

//! Set components of \c v to zero
inline void vmZero (VmVector2& v);
//! Set components of \c v to zero
inline void vmZero (VmVector3& v);

//! Set components of \c v to \c NaN
inline void vmNaN (VmVector2& v);
//! Set components of \c v to \c NaN
inline void vmNaN (VmVector3& v);

//! Sets  components of \c v
inline void vmSet (VmVector2& v, VmReal v0, VmReal v1);
//! Sets  components of \c v
inline void vmSet (VmVector3& v, VmReal v0, VmReal v1, VmReal v2);

//! Copies \c a to \c b
inline void vmCopy (const VmVector3& a, VmVector3& b);
//! Copies \c a to \c b
inline void vmCopy (const VmVector2& a, VmVector2& b);

//! returns whether all components of \c x are finite
inline bool isFinite (const VmVector3& x);
//! returns whether all components of \c x are finite
inline bool isFinite (const VmVector2& x);

//! Returns the dot product of \c v and \c w
inline VmReal vmDot (const VmVector2& v, const VmVector2& w);
//! Returns the dot product of \c v and \c w
inline VmReal vmDot (const VmVector3& v, const VmVector2& w);

//! Returns the length of v
inline VmReal vmLength (const VmVector2& v);
//! Returns the length of v
inline VmReal vmLength (const VmVector3& v);

//! Returns the squared length of v
inline VmReal vmLength2 (const VmVector2& v);
//! Returns the squared length of v
inline VmReal vmLength2 (const VmVector3& v);

//! Computes  \c result=a+b
inline void vmAdd (VmVector2& result, const VmVector2& a, const VmVector2& b);
//! Computes  \c result=a+b
inline void vmAdd (VmVector3& result, const VmVector3& a, const VmVector3& b);

//! Computes \c result+=v
inline void vmAdd (VmVector3& result, VmVector3& v);
//! Computes \c result+=v
inline void vmAdd (VmVector2& result, VmVector2& v);

//! Computes \c result*=scale
inline void vmScale (VmVector2& result, VmReal scale);
//! Computes \c result*=scale
inline void vmScale (VmVector3& result, VmReal scale);

// Computes \c result=v*scale
inline void vmScale (VmVector2& result, const VmVector2& v, VmReal scale);
// Computes \c result=v*scale
inline void vmScale (VmVector3& result, const VmVector3& v, VmReal scale);

//! Computes \c result+=scale*v
inline void vmAddScale (VmVector2& result, const VmVector2& v, VmReal scale);
//! Computes \c result+=scale*v
inline void vmAddScale (VmVector3& result, const VmVector3& v, VmReal scale);

//! Computes \c result=a-b
inline void vmSub (VmVector2& result, const VmVector2& a, const VmVector2& b);
//! Computes \c result=a-b
inline void vmSub (VmVector3& result, const VmVector3& a, const VmVector3& b);

//! Prints \c v
void vmPrint (const VmVector2& v);
//! Prints \c v
void vmPrint (const VmVector3& v);
//@}


//! \name Functions operating on matrices (and vectors)
//@{  

//! Set components of \c A to zero
inline void vmZero (VmMatrix2x2& A);
//! Set components of \c A to zero
inline void vmZero (VmMatrix2x3& A);
//! Set components of \c A to zero
inline void vmZero (VmMatrix3x2& A);
//! Set components of \c A to zero
inline void vmZero (VmMatrix3x3& A);
//! Set components of \c A to zero
inline void vmZero (VmMatrix4x4& A);


//! Set components of \c A to NaN
inline void vmNaN (VmMatrix2x2& A);
//! Set components of \c A to NaN
inline void vmNaN (VmMatrix2x3& A);
//! Set components of \c A to NaN
inline void vmNaN (VmMatrix3x2& A);
//! Set components of \c A to NaN
inline void vmNaN (VmMatrix3x3& A);
//! Set components of \c A to NaN
inline void vmNaN (VmMatrix4x4& A);

//! Set \c A to identity matrix
inline void vmOne (VmMatrix2x2& A);
//! Set \c A to identity matrix
inline void vmOne (VmMatrix3x3& A);
//! Set \c A to identity matrix
inline void vmOne (VmMatrix4x4& A);

//! Set components of \c A
inline void vmSet (VmMatrix2x2& A, VmReal a00, VmReal a01, VmReal a10, VmReal a11);
//! Set components of \c A
inline void vmSet (VmMatrix2x3& A, VmReal a00, VmReal a01, VmReal a02, VmReal a10, VmReal a11, VmReal a12);
//! Set components of \c A
inline void vmSet (VmMatrix3x2& A, VmReal a00, VmReal a01, VmReal a10, VmReal a11, VmReal a20, VmReal a21);
//! Set components of \c A
inline void vmSet (VmMatrix3x3& A, VmReal a00, VmReal a01, VmReal a02, VmReal a10, VmReal a11, VmReal a12, VmReal a20, VmReal a21, VmReal a22);

//! Solves a 2x2 \c Ax=b system of linear equations 
//! Returns the determinant of \c A
inline VmReal vmDeterminant (const VmMatrix2x2& A);
//! Returns the determinant of \c A
inline VmReal vmDeterminant (const VmMatrix3x3& A);

//! Returns the trace of \c A
inline VmReal vmTrace (const VmMatrix2x2& A);
//! Returns the trace of \c A
inline VmReal vmTrace (const VmMatrix3x3& A);

//! Copies \c A to \c B
inline void vmCopy (const VmMatrix2x2& A, VmMatrix2x2& B);
//! Copies \c A to \c B
inline void vmCopy (const VmMatrix3x3& A, VmMatrix3x3& B);
//! Copies \c A to \c B
inline void vmCopy (const VmMatrix2x3& A, VmMatrix2x3& B);
//! Copies \c A to \c B
inline void vmCopy (const VmMatrix3x2& A, VmMatrix3x2& B);
//! Copies \c A to \c B
inline void vmCopy (const VmMatrix4x4& A, VmMatrix4x4& B);

//! Returns, whether all components of \c A are finite
inline bool isFinite (const VmMatrix2x2& A);
//! Returns, whether all components of \c A are finite
inline bool isFinite (const VmMatrix2x3& A);
//! Returns, whether all components of \c A are finite
inline bool isFinite (const VmMatrix3x2& A);
//! Returns, whether all components of \c A are finite
inline bool isFinite (const VmMatrix3x3& A);

//! Computes  \c result=A+B
inline void vmAdd (VmMatrix2x2& result, const VmMatrix2x2& A, const VmMatrix2x2& B);
//! Computes  \c result=A+B
inline void vmAdd (VmMatrix3x3& result, const VmMatrix3x3& A, const VmMatrix3x3& B);


//! Computes \c result+=A
inline void vmAdd (VmMatrix2x2& result, const VmMatrix2x2& A);
//! Computes \c result+=A
inline void vmAdd (VmMatrix3x3& result, const VmMatrix3x3& A);

//! Computes \c result*=scale
inline void vmScale (VmMatrix2x2& result, VmReal scale);
//! Computes \c result*=scale
inline void vmScale (VmMatrix3x3& result, VmReal scale);

//! Computes \c result=A*scale
inline void vmScale (VmMatrix2x2& result, const VmMatrix2x2& A, VmReal scale);
//! Computes \c result=A*scale
inline void vmScale (VmMatrix3x3& result, const VmMatrix3x3& A, VmReal scale);

//! Computes \c result+=scale*A
inline void vmAddScale (VmMatrix2x2& result, const VmMatrix2x2& A, VmReal  scale);
//! Computes \c result+=scale*A
inline void vmAddScale (VmMatrix3x3& result, const VmMatrix3x3& A, VmReal  scale);

//! Computes \c result-=A
inline void vmSub (VmMatrix2x2& result, const VmMatrix2x2& A);
//! Computes \c result-=A
inline void vmSub (VmMatrix3x3& result, const VmMatrix3x3& A);
//! Computes \c result=A-B
inline void vmSub (VmMatrix2x2& result, const VmMatrix2x2& A, const VmMatrix2x2& B);
//! Computes \c result=A-B
inline void vmSub (VmMatrix3x3& result, const VmMatrix3x3& A, const VmMatrix3x3& B);

//! Computes \c result=A^t
inline void vmTranspose (VmMatrix2x2& result, const VmMatrix2x2& A);
//! Computes \c result=A^t
inline void vmTranspose (VmMatrix2x3& result, const VmMatrix3x2& A);
//! Computes \c result=A^t
inline void vmTranspose (VmMatrix3x2& result, const VmMatrix2x3& A);
//! Computes \c result=A^t
inline void vmTranspose (VmMatrix3x3& result, const VmMatrix3x3& A);    

//! Computes \c A=A^t
inline void vmTranspose (VmMatrix2x2& A);
//! Computes \c A=A^t
inline void vmTranspose (VmMatrix3x3& A);    

//! Computes \c result=A*v
inline void vmMultiply (VmVector2& result, const VmMatrix2x2& A, const VmVector2& v);
//! Computes \c result=A*v
inline void vmMultiply (VmVector2& result, const VmMatrix2x3& A, const VmVector3& v);
//! Computes \c result=A*v
inline void vmMultiply (VmVector3& result, const VmMatrix3x2& A, const VmVector2& v);
//! Computes \c result=A*v
inline void vmMultiply (VmVector3& result, const VmMatrix3x3& A, const VmVector3& v);

//! Computes result+=A*v
inline void vmMultiplyAdd (VmVector2& result, const VmMatrix2x2& A, const VmVector2& v);
//! Computes result+=A*v
inline void vmMultiplyAdd (VmVector2& result, const VmMatrix2x3& A, const VmVector3& v);
//! Computes result+=A*v
inline void vmMultiplyAdd (VmVector3& result, const VmMatrix3x2& A, const VmVector2& v);
//! Computes result+=A*v
inline void vmMultiplyAdd (VmVector3& result, const VmMatrix3x3& A, const VmVector3& v);

//! Computes result-=A*v
inline void vmMultiplySub (VmVector2& result, const VmMatrix2x2& A, const VmVector2& v);
//! Computes result-=A*v
inline void vmMultiplySub (VmVector2& result, const VmMatrix2x3& A, const VmVector3& v);
//! Computes result-=A*v
inline void vmMultiplySub (VmVector3& result, const VmMatrix3x2& A, const VmVector2& v);
//! Computes result-=A*v
inline void vmMultiplySub (VmVector3& result, const VmMatrix3x3& A, const VmVector3& v);

//! Computes result+=factor*A*v
inline void vmMultiplyAddScale (VmVector2& result, const VmMatrix2x2& A, const VmVector2& v, double factor);
//! Computes result+=factor*A*v
inline void vmMultiplyAddScale (VmVector2& result, const VmMatrix2x3& A, const VmVector3& v, double factor);
//! Computes result+=factor*A*v
inline void vmMultiplyAddScale(VmVector3& result, const VmMatrix3x2& A, const VmVector2& v, double factor);
//! Computes result+=factor*A*v
inline void vmMultiplyAddScale (VmVector3& result, const VmMatrix3x3& A, const VmVector3& v, double factor);

//! Computes \c result = A*B
inline void vmMultiply (VmMatrix2x2& result, const VmMatrix2x2& A, const VmMatrix2x2& B);
//! Computes \c result = A*B
inline void vmMultiply (VmMatrix2x3& result, const VmMatrix2x2& A, const VmMatrix2x3& B);
//! Computes \c result = A*B
inline void vmMultiply (VmMatrix2x2& result, const VmMatrix2x3& A, const VmMatrix3x2& B);
//! Computes \c result = A*B
inline void vmMultiply (VmMatrix2x3& result, const VmMatrix2x3& A, const VmMatrix3x3& B);
//! Computes \c result = A*B
inline void vmMultiply (VmMatrix3x2& result, const VmMatrix3x2& A, const VmMatrix2x2& B);
//! Computes \c result = A*B
inline void vmMultiply (VmMatrix3x3& result, const VmMatrix3x2& A, const VmMatrix2x3& B);
//! Computes \c result = A*B
inline void vmMultiply (VmMatrix3x2& result, const VmMatrix3x3& A, const VmMatrix3x2& B);
//! Computes \c result = A*B
inline void vmMultiply (VmMatrix3x3& result, const VmMatrix3x3& A, const VmMatrix3x3& B);



//! Computes \c result+=A*B
inline void vmMultiplyAdd (VmMatrix2x2& result, const VmMatrix2x2& A, const VmMatrix2x2& B);
//! Computes \c result+=A*B
inline void vmMultiplyAdd (VmMatrix3x3& result, const VmMatrix3x3& A, const VmMatrix3x3& B);

//! Computes \c result=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKt (double&    result, const   VmVector2& J, const VmMatrix2x2& C, const   VmVector2& K);
//! Computes \c result=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKt (double&    result, const   VmVector3& J, const VmMatrix3x3& C, const   VmVector3& K);
//! Computes \c result=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKt (VmMatrix2x2& result, const   VmVector2& J, double C, const   VmVector2& K);
//! Computes \c result=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKt (VmMatrix2x2& result, const VmMatrix2x2& J, const VmMatrix2x2& C, const VmMatrix2x2& K);
//! Computes \c result=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKt (VmMatrix2x2& result, const VmMatrix2x3& J, const VmMatrix3x3& C, const VmMatrix2x3& K);
//! Computes \c result=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKt (VmMatrix2x3& result, const VmMatrix2x2& J, const VmMatrix2x3& C, const VmMatrix3x3& K); // special
//! Computes \c result=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKt (VmMatrix3x3& result, const   VmVector3& J, double C, const VmVector3& K);
//! Computes \c result=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKt (VmMatrix3x3& result, const VmMatrix3x2& J, const VmMatrix2x2& C, const VmMatrix3x2& K);
//! Computes \c result=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKt (VmMatrix3x3& result, const VmMatrix3x3& J, const VmMatrix3x3& C, const VmMatrix3x3& K);

//! Computes \c result+=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKtAdd (double&    result, const   VmVector2& J, const VmMatrix2x2& C, const   VmVector2& K);
//! Computes \c result+=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKtAdd (double&    result, const   VmVector3& J, const VmMatrix3x3& C, const   VmVector3& K);
//! Computes \c result+=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKtAdd (VmMatrix2x2& result, const   VmVector2& J, double C, const   VmVector2& K);
//! Computes \c result+=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKtAdd (VmMatrix2x2& result, const VmMatrix2x2& J, const VmMatrix2x2& C, const VmMatrix2x2& K);
//! Computes \c result+=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKtAdd (VmMatrix2x2& result, const VmMatrix2x3& J, const VmMatrix3x3& C, const VmMatrix2x3& K);
//! Computes \c result+=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKtAdd (VmMatrix2x3& result, const VmMatrix2x2& J, const VmMatrix2x3& C, const VmMatrix3x3& K); // special
//! Computes \c result+=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKtAdd (VmMatrix3x3& result, const   VmVector3& J, double C, const   VmVector3& K);
//! Computes \c result+=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKtAdd (VmMatrix3x3& result, const VmMatrix3x2& J, const VmMatrix2x2& C, const VmMatrix3x2& K);
//! Computes \c result+=J*C*K^T
/*! \c result may not be identical to \c J, \c C and \c K
 */
inline void vmJCKtAdd (VmMatrix3x3& result, const VmMatrix3x3& J, const VmMatrix3x3& C, const VmMatrix3x3& K);

//! Computes \c result=A^-1
/*! The function returns the determinant. If the determinant is 0, the
    resulting matrix contains 'NaN' The symmetric functions are a
    little bit faster but naturally work only on symmetric matrices \c
    refine specifies, whether one step of iterative refinement shall
    be used. This is especially advisable with 3x3 matrices, since
    there is a lot of cancellation in the closed inverse
    formulas. However it takes more than twice the time.
*/
VmReal vmInverse (VmMatrix2x2& result, const VmMatrix2x2& A, bool refine=true);
//! Computes \c A=A^-1
/*! The function returns the determinant. If the determinant is 0, the
    resulting matrix contains 'NaN' The symmetric functions are a
    little bit faster but naturally work only on symmetric matrices \c
    refine specifies, whether one step of iterative refinement shall
    be used. This is especially advisable with 3x3 matrices, since
    there is a lot of cancellation in the closed inverse
    formulas. However it takes more than twice the time.
*/
VmReal vmInverse (VmMatrix2x2& A, bool refine=true);
//! Computes \c result=A^-1
/*! The function returns the determinant. If the determinant is 0, the
    resulting matrix contains 'NaN' The symmetric functions are a
    little bit faster but naturally work only on symmetric matrices \c
    refine specifies, whether one step of iterative refinement shall
    be used. This is especially advisable with 3x3 matrices, since
    there is a lot of cancellation in the closed inverse
    formulas. However it takes more than twice the time.
*/
VmReal vmInverse (VmMatrix3x3& result, const VmMatrix3x3& A, bool refine=true);
//! Computes \c A=A^-1
/*! The function returns the determinant. If the determinant is 0, the
    resulting matrix contains 'NaN' The symmetric functions are a
    little bit faster but naturally work only on symmetric matrices \c
    refine specifies, whether one step of iterative refinement shall
    be used. This is especially advisable with 3x3 matrices, since
    there is a lot of cancellation in the closed inverse
    formulas. However it takes more than twice the time.
*/
VmReal vmInverse (VmMatrix3x3& A, bool refine=true);

/*! The system is solved using determinants. 
    It must not be singular. The determinant
    is returned and may be checked for being to 
    close to 0 by the application.
*/
inline VmReal vmSolve (const VmMatrix2x2& A, VmVector2& x, const VmVector2& b);
//! Solves a 3x3 \c Ax=b system of linear equations 
/*! \sa vmSolve */
inline VmReal vmSolve (const VmMatrix3x3& A, VmVector3& x, const VmVector3& b);

//! Computes \c result=A^-1 for symmetric matrices
/*! The function returns the determinant. If the determinant is 0, the
    resulting matrix contains 'NaN' The symmetric functions are a
    little bit faster but naturally work only on symmetric matrices \c
    refine specifies, whether one step of iterative refinement shall
    be used. This is especially advisable with 3x3 matrices, since
    there is a lot of cancellation in the closed inverse
    formulas. However it takes more than twice the time.
*/
VmReal vmInverseSymmetric (VmMatrix2x2& result, const VmMatrix2x2& A, bool refine=true);
//! Computes \c A=A^-1 for symmetric matrices
/*! The function returns the determinant. If the determinant is 0, the
    resulting matrix contains 'NaN' The symmetric functions are a
    little bit faster but naturally work only on symmetric matrices \c
    refine specifies, whether one step of iterative refinement shall
    be used. This is especially advisable with 3x3 matrices, since
    there is a lot of cancellation in the closed inverse
    formulas. However it takes more than twice the time.
*/
VmReal vmInverseSymmetric (VmMatrix2x2& A, bool refine=true);
//! Computes \c result=A^-1 for symmetric matrices
/*! The function returns the determinant. If the determinant is 0, the
    resulting matrix contains 'NaN' The symmetric functions are a
    little bit faster but naturally work only on symmetric matrices \c
    refine specifies, whether one step of iterative refinement shall
    be used. This is especially advisable with 3x3 matrices, since
    there is a lot of cancellation in the closed inverse
    formulas. However it takes more than twice the time.
*/
VmReal vmInverseSymmetric (VmMatrix3x3& result, const VmMatrix3x3& A, bool refine=true);
//! Computes \c A=A^-1 for symmetric matrices
/*! The function returns the determinant. If the determinant is 0, the
    resulting matrix contains 'NaN' The symmetric functions are a
    little bit faster but naturally work only on symmetric matrices \c
    refine specifies, whether one step of iterative refinement shall
    be used. This is especially advisable with 3x3 matrices, since
    there is a lot of cancellation in the closed inverse
    formulas. However it takes more than twice the time.
*/
VmReal vmInverseSymmetric (VmMatrix3x3& A, bool refine=true);

//! Replace \c A by \c (A+A^T)/2
/*! This is useful if in theory a matrix should be symmetric. Then
    applying \c symmetrize after a computation may greatly reduce
    numerical problems.
*/
inline void vmSymmetrize (VmMatrix2x2& A);
//! Replace \c A by \c (A+A^T)/2
/*! This is useful if in theory a matrix should be symmetric. Then
    applying \c symmetrize after a computation may greatly reduce
    numerical problems.
*/
inline void vmSymmetrize (VmMatrix3x3& A);

//! Prints \c A
void vmPrint (const VmMatrix2x2& A);
//! Prints \c A
void vmPrint (const VmMatrix2x3& A);
//! Prints \c A
void vmPrint (const VmMatrix3x2& A);
//! Prints \c A
void vmPrint (const VmMatrix3x3& A);

//! Rotates a vector in the plane 
/*! 
     \f{equation}
        result = {\left(\begin{array}{cc}
            \cos  & -\sin \\ \sin  & \cos 
        \end{array}\right)} * v
     \f}

    \c result and \c v may be identical.
  
 */
inline void vmRotate (VmVector2& result, const VmVector2& v, double angle);

//! Computes \f$ R(\alpha)^TAR(\alpha) \f$ for a (not necessaryily symmetric) matrix A
/*! This so to say rotates a matrix both on the column and on the row side.
    
     \f{equation}
        R(\alpha) = {\left(\begin{array}{cc}
            \cos  & -\sin \\ \sin  & \cos 
        \end{array}\right)}
     \f}
*/
inline void vmRotateSym (VmMatrix2x2& A, double angle);
//! \c rotateSym(A,cos(angle),sin(angle)) is the same as \c rotateSym(A, angle)
inline void vmRotateSym (VmMatrix2x2& A, double cAngle, double sAngle);

//! Computes a cholesky decomposition \c A=LLt for a symmetric positive semidefinite A.
/*! If the computation fails, the matrix is not S.P.D. and 'false' is returned.
    Diagonal entries between \c [-eps..0] are treated as 0 leading to a valid
    decomposition. Even diagonal entries infinity are allowed.
*/
bool vmCholesky (const VmMatrix2x2& A, VmMatrix2x2& L, double eps=0);
//! Computes a cholesky decomposition \c A=LLt for a symmetric positive semidefinite A.
/*! If the computation fails, the matrix is not S.P.D. and 'false' is returned.
    Diagonal entries between \c [-eps..0] are treated as 0 leading to a valid
    decomposition. Even diagonal entries infinity are allowed.
*/
bool vmCholesky (const VmMatrix3x3& A, VmMatrix3x3& L, double eps=0);

//! Computes the inverse of the Cholesky decomposition of A
/*! Even diagonal entries infinity are allowed leading to 0 in the inverse.
 */
bool vmCholeskyInverse (const VmMatrix2x2& A, VmMatrix2x2& LI, double eps=0);
//! Computes the inverse of the Cholesky decomposition of A
/*! Even diagonal entries infinity are allowed leading to 0 in the inverse.
 */
bool vmCholeskyInverse (const VmMatrix3x3& A, VmMatrix3x3& LI, double eps=0);

//! Computes an eigenvalue decomposition of a symmetric 2x2 matrix 
/*! \c A must be symmetric. The decomposition is stable even for
    degenerated eigenvalues \f$ \left(\begin{array}{c} \cos\phi\\
    \sin\phi \end{array}\right) \f$ will be an eigenvector with
    eigenvalue \c l0 and \f$ \left(\begin{array}{c} -\sin\phi\\
    \cos\phi \end{array}\right) \f$ will be an eigenvector with
    eigenvalue \c l1. \c l0>=l1.  If \c posDef is set, \c A must be
    positive definite. The algorithm asserts, that \c l0 and \c l1 are \c >=0
    with some suitable tolerance and clips them to 0 if they are a
    little below.
*/
void vmEigenDecompositionSymmetric (const VmMatrix2x2& A, VmReal& phi, VmReal& l0, VmReal& l1, bool posDef=false);
//! Computes an eigenvalue decomposition for a symmetric 3x3 matrix \c A
/*! The eigenvalues are returned in \c d with decreasing values. \c
    A=QDQ^T. \c Q is a orthogonal matrix, where the i-th column
    contains the normalized eigenvector for \c d[i].  This is done by
    Jacobi rotations specialized for dimension 3 (thus much
    faster). The algorithm terminates, if the sum of the off-diagonal
    elements is <= \c tolerance. \c tolerance may be 0, which will
    result in convergence to numerical underflow. 
*/
void vmEigenDecompositionSymmetric (const VmMatrix3x3& A, VmMatrix3x3& Q, VmReal d[3], VmReal tolerance=0);


//! Computes a singular value decomposition \c A=U*diag(s)*V 
/*! The singular values in \c s are descending and \c U and \c V are
    orthonormal. This is done forming \c AtA and calling
    \c eigenDecompositionSymmetric, which is numerically no wise but
    should work for 3x3 matrices. It has to be expected however, that
    the precision of \c s is smaller than achievable by better
    algorithm.
*/
void vmSingularValueDecomposition (const VmMatrix3x3& A, VmMatrix3x3& U, VmMatrix3x3& V, VmReal s[3], VmReal tolerance=0);

//! Computes the QR decomposition of A
/*! \c A=QR, where \c R is upper triangular and \c Q is orthonormal
 */
void vmQRDecomposition (const VmMatrix3x3& A, VmMatrix3x3& Q, VmMatrix3x3& R);
//@}


//!\name Functions operating on rigid body transformations
//@{

//! Interpretes \c T as a transform and applies it to a vector describing a point
inline void vmApplyTransformToPoint (VmVector3& result, const VmMatrix4x4& T, const VmVector3& p);
//! Interpretes \c T as a transform and applies \c T^{-1} to a vector describing a point
inline void vmApplyInverseTransformToPoint (VmVector3& result, const VmMatrix4x4& T, const VmVector3& p);
//! Interpretes \c T as a transform and applies it to a vector describing a direction
inline void vmApplyTransformToDirection (VmVector3& result, const VmMatrix4x4& T, const VmVector3& d);
//! Interpretes \c T as a transform and applies \c T^{-1} to a vector describing a direction
inline void vmApplyInverseTransformToDirection (VmVector3& result, const VmMatrix4x4& T, const VmVector3& d);
//! Computes the inverse of \c A
/*! Note, that this routine compute \f$ A^{-1} \f$ only, if \c A is a valid transform,
    i.e. the left/upper 3*3 matrix is orthonormal and the lower row is \c (0,0,0,1)
*/
inline void vmInverseTransform (VmMatrix4x4& result, const VmMatrix4x4& A);
//! Multiplies to transformations \f$ result=A*B \f$
inline void vmMultiply (VmMatrix4x4& result, const VmMatrix4x4& A, const VmMatrix4x4& B);

//@}


//! Returns a rigid body transformation, that maps \c p0a to \c p1a and \c p0b to \c p1b
/*! The transformation consist of first rotating by \c dO around the origin and than
    translating by \c (dX,dY). The transformation yield only an exact result, if the distance
    \c p0a-p0b is the same as \c p1a-p1b. Otherwise a least square solution is returned.
*/
void vmRigidBodyTransformation (const VmVector2& p0a, const VmVector2& p0b, const VmVector2& p1a, const VmVector2& p1b, 
				       VmReal& dX, VmReal& dY, VmReal& dO);

//! Overloaded function
void vmRigidBodyTransformation (VmReal p0aX, VmReal p0aY, VmReal p0bX, VmReal p0bY, VmReal p1aX, VmReal p1aY, VmReal p1bX, VmReal p1bY,
				       VmReal& dX, VmReal& dY, VmReal& dO);






//! Internal routine for \c eigenDecompositionSymmetric3x3
/*! It computes the parameter \c s, tau, t for one Jacobi rotation
    eliminating an off diagonal element of value \c acpq=ac[p][q], 
    with \c deltaD=d[q]-d[p]
*/
void vmOneJacobiParameter (VmReal acpq, VmReal deltaD, VmReal& s, VmReal& tau, VmReal& t);








#include "vectormath_i.h"


#endif  /* VECTORMATH_H */
