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
/*! \file vmLSTransform.cc
    \brief Implementation of \c vmLSTransform.h
    \author Udo Frese

    Contains the implementation of least square routines for matching
    point-pair-clouds.
*/
#include "vmLSTransform.h"
#include <assert.h>

double vmLeastSquareTransform (const VmPointPair3List& p, VmMatrix4x4& T, double& evRatio)
{
  // First compute the weighted mean to move the COG of both point clouds to 0
  VmVector3 meanA, meanB;
  vmZero (meanA);
  vmZero (meanB);  
  double weightSum = 0;
  for (int i=0; i<(int) p.size(); i++) {
    const VmPointPair3& pp = p[i];    
    vmAddScale (meanA, pp._a, pp._weight);
    vmAddScale (meanB, pp._b, pp._weight);
    weightSum += pp._weight;    
  }
  vmScale (meanA, meanA, 1/weightSum);
  vmScale (meanB, meanB, 1/weightSum);  

  // Now compute the correlation matrix H and the fixed part of the chi^2 error
  VmMatrix3x3 H;
  vmZero (H);  
  double ssd=0;  
  for (int i=0; i<(int) p.size(); i++) {
    const VmPointPair3& pp = p[i];
    VmVector3 aCentered, bCentered;
    vmSub (aCentered, pp._a, meanA);
    vmSub (bCentered, pp._b, meanB);
    vmJCKtAdd (H, aCentered, pp._weight, bCentered);
    ssd += pp._weight * (vmLength2 (aCentered) + vmLength2 (bCentered));    
  }

  // And find the optimal rotation matrix maximizing Trance(RH)
  // by computing the SVD decomposition H=USV and setting R=VU
  VmMatrix3x3 U, V, R;
  double s[3];  
  vmSingularValueDecomposition (H, U, V, s);
  if (s[0]>0) evRatio = s[1]/s[0];
  else evRatio = 0;  
  vmMultiply (R, U, V);
  vmTranspose (R);
  
  // Handle the special case of a left-hand R
  if (vmDeterminant (R)<0) {
    U[0][2] = -U[0][2];
    U[1][2] = -U[1][2];
    U[2][2] = -U[2][2];
    vmMultiply (R, U, V);
    vmTranspose (R);
  }
  
  // And find the optimal translation
  VmVector3 RAMean;
  vmMultiply (RAMean, R, meanA);  
  T[0][0] = R[0][0]; T[0][1] = R[0][1]; T[0][2] = R[0][2];
  T[1][0] = R[1][0]; T[1][1] = R[1][1]; T[1][2] = R[1][2];
  T[2][0] = R[2][0]; T[2][1] = R[2][1]; T[2][2] = R[2][2];
  T[0][3] = meanB[0] - RAMean[0];
  T[1][3] = meanB[1] - RAMean[1];
  T[2][3] = meanB[2] - RAMean[2];
  T[3][0] = T[3][1] = T[3][2] = 0; T[3][3] = 1;  

  // Finally compute the chi^2 error
  return ssd - 2*(s[0]+s[1]+s[2]);  
}



double vmLeastSquareTransform ( const VmPointPair2List& p, VmVector2& d, double& phi)
{
    d[0] = d[1]  = phi = 0;
    int n = (int) p.size();
    if (n==0) return 0;
    
    // The optimal least square solution under pure translation is
    // to move the COGs onto each other. So we move the COG both of
    // a and b to 0 by removing its mean. Then rotation around 0 has
    // no influence on the COG and thus finding the optimal rotation
    // is decoupled from finding the optimal translatio.
    VmVector2 aMean, bMean, optRot;
    aMean[0] = aMean[1] = 0;
    bMean[0] = bMean[1] = 0;
    double totalWeight = 0;
    for (int i=0; i<n; i++) {
        const VmPointPair& pp = p[i];
        if (pp._weight>0) {
            aMean[0] += pp._weight * pp._a[0];
            aMean[1] += pp._weight * pp._a[1];
            bMean[0] += pp._weight * pp._b[0];
            bMean[1] += pp._weight * pp._b[1];
            totalWeight += pp._weight;
        }
    }
    aMean[0] /= totalWeight;
    aMean[1] /= totalWeight;
    bMean[0] /= totalWeight;
    bMean[1] /= totalWeight;
    if (totalWeight<=0) return 0;

    double ssd = 0;
    optRot[0] = optRot[1] = 0;
    if (n>1) {
        // To minimize 'sum (b[i]-Rot(phi)a[i])^2) we expand:
        //    = sum b[i]^2 + sum a[i]^2 + sum b[i]Rot(phi)a[i]
        // The first two sums are constant. We write the last sum in
        // terms of c = cos(phi), s=sin(phi)
        //  = sum (b[i][0]*c*a[i][0] + b[i][0]*(-s)*a[i][1] + b[i][1]*s*a[i][0] + b[i][1]*c*a[i][1]
        // grouping by sin and cos yields:
        //  = c*(sum b[i][0]*a[i][0]+b[i][1]*a[i][1]) + s*(sum b[i][1]*a[i][0] -b[i][0]*a[i][1])
        // Minimizing this expression means, making (c,s) antiparralel to the sums.
        for (int i=0; i<n; i++) {
            const VmPointPair& pp = p[i];
            if (pp._weight>0) {
                VmVector2 aS, bS;
                aS[0] = pp._a[0] - aMean[0];
                aS[1] = pp._a[1] - aMean[1];
                bS[0] = pp._b[0] - bMean[0];
                bS[1] = pp._b[1] - bMean[1];
                ssd += pp._weight * (aS[0]*aS[0] + aS[1]*aS[1] + bS[0]*bS[0] + bS[1]*bS[1]);
                optRot[0] += pp._weight * (-bS[0]*aS[0] - bS[1]*aS[1]);
                optRot[1] += pp._weight * (-bS[1]*aS[0] + bS[0]*aS[1]);
            }
        }
        if (optRot[0] != 0 || optRot[1] != 0) phi = atan2 (-optRot[1], -optRot[0]);
    }
    // Compute 'd': First shift aMean to 0. Then 0 to bMean.
    double c = cos(phi), s=sin(phi);
    d[0] = -c*aMean[0] +s*aMean[1] + bMean[0];
    d[1] = -s*aMean[0] -c*aMean[1] + bMean[1];
    // Compute SSD error
    ssd += 2*(c*optRot[0] + s*optRot[1]);
    assert (!isnan(ssd));
    return ssd;
}




//****************** VmPointPair2 ***************

VmPointPair2::VmPointPair2 ()
  :_weight(0)
{
  vmZero (_a);
  vmZero (_b);  
}

VmPointPair2::VmPointPair2 (const VmVector2& a, const VmVector2& b, double weight)
        :_weight(weight)
{
    _a[0] = a[0];
    _a[1] = a[1];
    _b[0] = b[0];
    _b[1] = b[1];
}


VmPointPair2::VmPointPair2 (double ax, double ay, double bx, double by, double weight)
        :_weight(weight)
{
    _a[0] = ax;
    _a[1] = ay;
    _b[0] = bx;
    _b[1] = by;    
}


//****************** VmPointPair3 ***************

VmPointPair3::VmPointPair3 ()
  :_weight(0)
{
  vmZero (_a);
  vmZero (_b);  
}

VmPointPair3::VmPointPair3 (const VmVector3& a, const VmVector3& b, double weight)
        :_weight(weight)
{
  vmCopy (a, _a);
  vmCopy (b, _b);  
}


