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

/*! \file vectormath.cc
    \brief Implementation of \c vectormath.h
    \author Udo Frese

    Contains the implementation of small matrix and vector routines.
 */
#include "vectormath.h"

#include <assert.h>
#include <stdexcept>

void vmPrint (const VmVector2& v)
{
    printf("{ %10e %10e }\n", v[0], v[1]);
}


void vmPrint (const VmVector3& v)
{
    printf("{ %10e %10e %10e }\n", v[0], v[1], v[2]);
}


void vmPrint (const VmMatrix2x2& A)
{
    printf("{{ %10e %10e }\n", A[0][0], A[0][1]);
    printf(" { %10e %10e }}\n", A[1][0], A[1][1]);
}


void vmPrint (const VmMatrix2x3& A)
{
    printf("{{ %10e %10e %10e }\n", A[0][0], A[0][1], A[0][2]);
    printf(" { %10e %10e %10e }}\n", A[1][0], A[1][1], A[1][2]);
}


void vmPrint (const VmMatrix3x2& A)
{
    printf("{{ %10e %10e }\n" , A[0][0], A[0][1]);
    printf(" { %10e %10e }\n" , A[1][0], A[1][1]);
    printf(" { %10e %10e }}\n", A[2][0], A[2][1]);
}


void vmPrint (const VmMatrix3x3& A)
{
    printf("{{ %10e %10e %10e }\n", A[0][0], A[0][1], A[0][2]);
    printf(" { %10e %10e %10e }\n", A[1][0], A[1][1], A[1][2]);
    printf(" { %10e %10e %10e }}\n", A[2][0], A[2][1], A[2][2]);
}


void vmRotateSym (VmMatrix2x2& A, double angle)
{
    vmRotateSym (A, cos(angle), sin(angle));
}


//! Internal subroutine for \c vmEigenDecompositionSymmetric
/*! Computes the parameters for one Jacobi rotation */
void vmOneJacobiParameter (VmReal acpq, VmReal deltaD, VmReal& s, VmReal& tau, VmReal& t, VmReal& h)
{
    VmReal g, theta, c;
    g = 100*fabs(acpq);
    if (g+fabs(h)==fabs(deltaD)) t = acpq/deltaD;
    else {
        theta = 0.5*deltaD/(acpq);
        t     = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
        if (theta<0) t=-t;
    }
    c=1.0/sqrt(1+t*t);
    s=t*c;
    tau=s/(1.0+c);
}




bool vmCholesky (const VmMatrix3x3& A, VmMatrix3x3& L, double eps)
{
    double sum;
    sum = A[0][0];
    if (sum<-eps) return false;
    if (sum>0) {
        // (0,0)
        L[0][0] = sqrt(sum);
    
        // (1,0)
        sum = A[1][0];
        L[1][0] = sum/L[0][0];

        // (2,0)
        sum = A[2][0];
        L[2][0] = sum/L[0][0];
    }
    else L[0][0] = L[1][0] = L[2][0] = 0;

    // (0,1)
    L[0][1] = 0;

    sum = A[1][1] - L[1][0]*L[1][0];
    if (sum<-eps) return false;
    if (sum>0) {
        // (1,1)
        L[1][1] = sqrt(sum);

        // (2,1)
        sum = A[2][1] - L[1][0]*L[2][0];
        L[2][1] = sum / L[1][1];
    }
    else L[1][1] = L[2][1] = 0;
    
    // (0,2)
    L[0][2] = 0;

    // (1,2)
    L[1][2] = 0;
    
    sum = A[2][2] - L[2][0]*L[2][0] - L[2][1]*L[2][1];
    if (sum<-eps) return false;
    if (sum>0) {
        // (2,2)
        L[2][2] = sqrt(sum);
    }
    else L[2][2] = 0;

    return true;
}



bool vmCholeskyInverse (const VmMatrix3x3& A, VmMatrix3x3& LI, double eps)
{
  if (!vmCholesky (A, LI, eps)) return false;
  // Inverse of triangular matrix
  if (LI[0][0]==0 || LI[1][1]==0 || LI[2][2]==0) return false;  
  LI[0][0] = 1/LI[0][0];
  LI[1][0] = -LI[0][0]*LI[1][0]/LI[1][1];
  LI[2][0] = (-LI[0][0]*LI[2][0]-LI[1][0]*LI[2][1])/LI[2][2];
  LI[0][1] = 0;
  LI[1][1] = 1/LI[1][1];
  LI[2][1] = -LI[1][1]*LI[2][1]/LI[2][2];
  LI[0][2] = 0;
  LI[1][2] = 0;
  LI[2][2] = 1/LI[2][2];  
  return true;  
}


bool vmCholesky (const VmMatrix2x2& A, VmMatrix2x2& L, double eps)
{
    double sum;
    sum = A[0][0];
    if (sum<-eps) return false;
    if (sum>0) {
        // (0,0)
        L[0][0] = sqrt(sum);
        // (1,0)
        sum = A[1][0];
        L[1][0] = sum/L[0][0];
    }
    else L[0][0] = L[1][0] = 0;
    
    
    // (0,1)
    L[0][1] = 0;
    
    
    // (1,1)
    sum = A[1][1] - L[1][0]*L[1][0];
    if (sum<-eps) return false;

    if (sum>0) L[1][1] = sqrt(sum);
    else L[1][1] = 0;

    return true;
}


bool vmCholeskyInverse (const VmMatrix2x2& A, VmMatrix2x2& LI, double eps)
{
  if (!vmCholesky (A, LI, eps)) return false;
  if (LI[0][0]==0 || LI[1][1]==0) return false;  
  // Inverse of triangular matrix
  LI[0][0] = 1/LI[0][0];
  LI[1][0] = -LI[0][0]*LI[1][0]/LI[1][1];
  LI[0][1] = 0;
  LI[1][1] = 1/LI[1][1];  
  return true;  
}


//! Internal MACRO for eigenDecompositionSymmetric3x3 performing one JACOBI rotation
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);
//! Internal MACRO for eigenDecompositionSymmetric3x3 
#define ADAPT(ip, iq) h = t*acpq; z[ip] -= h; z[iq] += h; d[ip] -= h; d[iq] += h; ac[ip][iq]=0.0;
//! Internal MACRO for eigenDecompositionSymmetric3x3 
#define POSTADAPT         b[0] += z[0]; d[0] = b[0]; z[0] = 0; b[1] += z[1]; d[1] = b[1]; z[1] = 0; b[2] += z[2]; d[2] = b[2]; z[2] = 0;
//! Internal MACRO for eigenDecompositionSymmetric3x3 swapping to numbers
#define SWAP(x1, x2) park=(x1); (x1)=(x2); (x2)=park;


void vmEigenDecompositionSymmetric (const VmMatrix3x3& a, VmMatrix3x3& q, VmReal d[3], VmReal tolerance)
{
    assert (a[0][1]==a[1][0] && a[1][2]==a[2][1] && a[0][2]==a[2][0]); // Must be symmetric
    // Using Jacobi rotations, which will only need 3 rotations per sweep
    // See numerical recipes chapter 11.1 for documentation of the method
    q[0][1] = q[0][2] = q[1][0] = q[1][2] = q[2][0] = q[2][1] = 0;
    q[0][0] = q[1][1] = q[2][2] = 1.0;
    VmReal ac[3][3];
    ac[0][0] = a[0][0];
    ac[0][1] = a[0][1];
    ac[0][2] = a[0][2];
    ac[1][1] = a[1][1];
    ac[1][2] = a[1][2];
    ac[2][2] = a[2][2]; // We use only above the diagonal

    VmReal b[3];
    b[0] = d[0] = ac[0][0];
    b[1] = d[1] = ac[1][1];
    b[2] = d[2] = ac[2][2];

    VmReal z[3];
    z[0] = z[1] = z[2] = 0;
    
    int ctr;
    VmReal t, tau, sm, s, h, g, acpq;
    
    for (ctr=0; ctr<20; ctr++) {
        sm = fabs(ac[0][1]) + fabs(ac[0][2]) + fabs(ac[1][2]);  // Off-Diagonal sum
        if (sm<=tolerance) break; // 0.0 is o.K. due to quadratic convergence to underflow
        // Eliminate ac[0][1]
        acpq = ac[0][1];
        vmOneJacobiParameter (acpq, d[1]-d[0], s, tau, t, h);
        ADAPT(0, 1)
            ROTATE(ac,0,2,1,2)
            ROTATE(q,0,0,0,1)
            ROTATE(q,1,0,1,1)
            ROTATE(q,2,0,2,1)
        // Eliminate ac[0][2]
        acpq = ac[0][2];
        vmOneJacobiParameter (acpq, d[2]-d[0], s, tau, t, h);
        ADAPT(0, 2)
            ROTATE(ac,0,1,1,2)
            ROTATE(q,0,0,0,2)
            ROTATE(q,1,0,1,2)
            ROTATE(q,2,0,2,2)
        // Eliminate ac[1][2]
        acpq = ac[1][2];
        vmOneJacobiParameter (acpq, d[2]-d[1], s, tau, t, h);
        ADAPT(1, 2)
            ROTATE(ac,0,1,0,2)
            ROTATE(q,0,1,0,2)
            ROTATE(q,1,1,1,2)
            ROTATE(q,2,1,2,2)
        POSTADAPT
    };
    assert (ctr<50);
// Error: Too many iterations without convergence
    // Now sort the eigenvalues
    if (d[0]<d[1] && d[2]<d[1]) { // Swap max=1 to 0
        VmReal park;
        SWAP (d[0], d[1])
        SWAP (q[0][0], q[0][1])
        SWAP (q[1][0], q[1][1])
        SWAP (q[2][0], q[2][1])
    }
    if (d[0]<d[2] && d[1]<d[2]) { // Swap max=2 to 0
        VmReal park;
        SWAP (d[0], d[2])
        SWAP (q[0][0], q[0][2])
        SWAP (q[1][0], q[1][2])
        SWAP (q[2][0], q[2][2])
    }
    if (d[1]<d[2]) { // Swap larger of 1,2 to 1
        VmReal park;
        SWAP (d[1], d[2])
        SWAP (q[0][1], q[0][2])
        SWAP (q[1][1], q[1][2])
        SWAP (q[2][1], q[2][2])
    }
    
#ifndef NDEBUG
    // Check, whether the decomposition is correct by vmMultiplying the equation
    if (!(d[0]>=d[1] && d[1]>=d[2])) {     // SPECIAL DEBUG CODE to be deleted if I found the bug
/*        printf ("A=\n%10.3f %10.3f %10.3f\n%10.3f %10.3f %10.3f\n%10.3f %10.3f %10.3f\n",
                 a[0][0], a[0][1], a[0][2],
                 a[1][0], a[1][1], a[1][2],
                 a[2][0], a[2][1], a[2][2]);
                 printf ("d=%10.3f %10.3f %10.3f\n", d[0], d[1], d[2]);*/
    }
    assert (d[0]>=d[1] && d[1]>=d[2]);
    VmReal res[3][3];
    VmReal diff=0;
    for (int i=0; i<3; i++) for (int j=0; j<3; j++) {
        res[i][j] = d[0]*q[i][0]*q[j][0] + d[1]*q[i][1]*q[j][1] + d[2]*q[i][2]*q[j][2];
        diff += fabs(res[i][j]-a[i][j]);
    }
    assert (diff<=1E-6*(fabs(d[0])+fabs(d[1])+fabs(d[2])));
#endif
}

//! Internal subroutine for vmQRDecomposion
/*! Applies a householder transform paramerized by \c u to A
     \f{equation}
           A \cdot \left(I-2\frac{uu^T}{u^Tu}\right) = A + (\frac{-2}{u^Tu})\cdot(Au)u^T
     \f}
 */
void vmATimesHouseholder (VmMatrix3x3& A, const VmVector3& u)
{
  // 
  double gamma = (u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);  
  if (gamma==0) return;
  gamma = -2/gamma;  
  VmVector3 v;
  vmMultiply (v, A, u);  
  vmJCKtAdd (A, v, gamma, u);  
}

//! Internal MACRO for vmQRDecomposion making a diagonal entry of R positive by negating the column of Q
#define MAKEPOSITIVE(i) if (R[i][i]<0) {R[i][0] = -R[i][0]; R[i][1] = -R[i][1]; R[i][2] = -R[i][2]; \
Q[0][i] = -Q[0][i]; Q[1][i] = -Q[1][i]; Q[2][i] = -Q[2][i];}                                                                     \
  


void vmQRDecomposition (const VmMatrix3x3& A, VmMatrix3x3& Q, VmMatrix3x3& R)
{
  VmMatrix3x3 Rt;  
  vmTranspose (Rt, A);
  vmOne (Q);  
  VmVector3 u;
  double alpha;

  // Column 0 (i.e. row 0 of Rt)
  alpha = sqrt(Rt[0][0]*Rt[0][0] + Rt[0][1]*Rt[0][1] + Rt[0][2]*Rt[0][2]);
  if (Rt[0][0]>0) alpha = -alpha;
  u[0] = Rt[0][0] - alpha;
  u[1] = Rt[0][1];
  u[2] = Rt[0][2];
  vmATimesHouseholder (Q, u);
  vmATimesHouseholder (Rt, u);

  // Column 1 (i.e. row 1 of Rt)
  alpha = sqrt(Rt[1][1]*Rt[1][1] + Rt[1][2]*Rt[1][2]);
  if (Rt[1][1]>0) alpha = -alpha;
  u[0] = 0;
  u[1] = Rt[1][1] - alpha;
  u[2] = Rt[1][2];
  vmATimesHouseholder (Q, u);
  vmATimesHouseholder (Rt, u);

  vmTranspose (R, Rt);  

  // Now make R's diagonal positive
  MAKEPOSITIVE (0);  
  MAKEPOSITIVE (1);
  MAKEPOSITIVE (2);  
}

//! Internal MACRO for vmSingularValueDecomposition
/*! It swaps to eigenvalues in s and the corresponding
    rows V and columns in U, such that U*S*V stays
    the same.
*/
#define SWAPROWANDCOL(i,j)                     \
  {                                             \
    SWAP(s[i],s[j]);                            \
    SWAP(U[0][i], U[0][j]);                     \
    SWAP(U[1][i], U[1][j]);                     \
    SWAP(U[2][i], U[2][j]);                     \
    SWAP(V[i][0], V[j][0]);                     \
    SWAP(V[i][1], V[j][1]);                     \
    SWAP(V[i][2], V[j][2]);                     \
  }                                             \
  


void vmSingularValueDecomposition (const VmMatrix3x3& A, VmMatrix3x3& U, VmMatrix3x3& V, VmReal s[3], VmReal tolerance)
{
  VmMatrix3x3 Sigma;
  vmCopy (A, Sigma);
  vmOne (U);
  vmOne (V);
  VmReal diff, norm;  
  int ctr=0;
  
  do {
    VmMatrix3x3 Q, R, hv;
    vmQRDecomposition (Sigma, Q, R);
    vmMultiply (hv, U, Q);
    vmCopy (hv, U);
    vmTranspose (Sigma, R);

    vmQRDecomposition (Sigma, Q, R);
    vmTranspose (Q);
    vmMultiply (hv, Q, V);
    vmCopy (hv, V);
    vmTranspose (Sigma, R);    
    
    diff = fabs(Sigma[0][1]) + fabs(Sigma[0][2]) + fabs(Sigma[1][0]) + fabs(Sigma[1][2]) + fabs(Sigma[2][0]) + fabs(Sigma[2][1]);    
    norm = fabs(Sigma[0][0]) + fabs(Sigma[1][1]) + fabs(Sigma[2][2]);    
    ctr++;    
  } while (diff>norm*1E-10);

  s[0] = Sigma[0][0];
  s[1] = Sigma[1][1];
  s[2] = Sigma[2][2];  

  double park;  
  if (s[0]<s[1]) SWAPROWANDCOL(0,1);
  if (s[0]<s[2]) SWAPROWANDCOL(0,2);
  if (s[1]<s[2]) SWAPROWANDCOL(1,2);  
  

  //if (ctr>100) fprintf (stderr, "SVD needed %d iterations for %f/%f (s=%f,%f,%f)\n", ctr, diff, norm, s[0], s[1], s[2]);  

#ifndef NDEBUG
    // Check, whether the decomposition is correct by multiplying the equation
    assert (s[0]>=s[1] && s[1]>=s[2]);
    VmReal res[3][3];
    diff=0;
    for (int i=0; i<3; i++) for (int j=0; j<3; j++) {
      res[i][j] = s[0]*U[i][0]*V[0][j] + s[1]*U[i][1]*V[1][j] + s[2]*U[i][2]*V[2][j];
      double vxx = V[0][i]*V[0][j] + V[1][i]*V[1][j] +V[2][i]*V[2][j];      
      double uxx = U[0][i]*U[0][j] + U[1][i]*U[1][j] +U[2][i]*U[2][j];
      if (i==j) {
        vxx-=1;
        uxx-=1;
      }      
      diff += fabs(res[i][j]-A[i][j]) + fabs(vxx) + fabs(uxx);
    }
    assert (diff<=1E-4*(fabs(s[0])+fabs(s[1])+fabs(s[2])));
#endif
}



void vmEigenDecompositionSymmetric (const VmMatrix2x2& a, VmReal& phi, VmReal& l0, VmReal& l1, bool posDef)
{
    if (!isFinite(a)) {
        phi = vmNaN();
        l0  = vmNaN();
        l1  = vmNaN();
        return;
    }
    VmReal ll0, ll1;
    phi = atan2(-2*a[1][0], a[1][1]-a[0][0])/2;
    if (!isfinite(phi)) phi=0;
    VmReal c = cos(phi);
    VmReal s = sin(phi);
    VmReal c2 = c*c, s2=s*s, cs = c*s;
    
    ll0 = c2*a[0][0]+2*cs*a[1][0]+s2*a[1][1];
    ll1 = s2*a[0][0]-2*cs*a[1][0]+c2*a[1][1];
    if (ll0<ll1) {
        if (phi>0) phi-=M_PI/2;
        else phi+=M_PI/2;
        l0 = ll1;
        l1 = ll0;
    }
    else {
        l0 = ll0;
        l1 = ll1;
    }
    double limit = (fabs(l0)+fabs(l1))*1E-3;
    if (posDef) {
        assert (l0>=-limit && l1>=-limit);
        if (l0<0) l0=0;
        if (l1<0) l1=0;
    }
    c = cos(phi);
    s = sin(phi);
    VmReal dd = fabs(a[0][0]*c + a[0][1]*s - c*l0);
    assert (dd<=limit);
    dd = fabs(a[1][0]*c + a[1][1]*s - s*l0);
    assert(dd<=limit); // 'c,s' is  eigenvector with eigenvalue l0

    dd = fabs(-a[0][0]*s + a[0][1]*c + s*l1);
    assert (dd<=limit);
    dd = fabs(-a[1][0]*s + a[1][1]*c - c*l1); // '-s,c' is  eigenvector with eigenvalue l1
    assert (dd<=limit);
}


void vmRigidBodyTransformation (const VmVector2& p0a, const VmVector2& p0b, const VmVector2& p1a, const VmVector2& p1b, 
                                                VmReal& dX, VmReal& dY, VmReal& dO)
{
    VmVector2 d0, d1;
    vmSub (d0, p0b, p0a);
    vmSub (d1, p1b, p1a);
    dO = atan2 (d1[1], d1[0]) - atan2(d0[1], d0[0]);
    double c = cos(dO), s = sin(dO);
    VmVector2 c0, c1, c0R;
    c0[0] = (p0a[0]+p0b[0])/2;
    c0[1] = (p0a[1]+p0b[1])/2;
    c1[0] = (p1a[0]+p1b[0])/2;
    c1[1] = (p1a[1]+p1b[1])/2;
    c0R[0]   = c*c0[0] - s*c0[1];
    c0R[1]   = s*c0[0] + c*c0[1];
    dX = c1[0] - c0R[0];
    dY = c1[1] - c0R[1];
}


void vmRigidBodyTransformation (VmReal p0aX, VmReal p0aY, VmReal p0bX, VmReal p0bY, VmReal p1aX, VmReal p1aY, VmReal p1bX, VmReal p1bY,
                                                VmReal& dX, VmReal& dY, VmReal& dO)
{
    VmVector2 p0a, p0b, p1a, p1b;
    vmSet(p0a, p0aX, p0aY);
    vmSet(p0b, p0bX, p0bY);
    vmSet(p1a, p1aX, p1aY);
    vmSet(p1b, p1bX, p1bY);
    vmRigidBodyTransformation (p0a, p0b, p1a, p1b, dX, dY, dO);
}






