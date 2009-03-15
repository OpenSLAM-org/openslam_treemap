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
/*! \file vmLSTransform.h 
    \brief Least square matching of clouds of point pairs
    \author Udo Frese

    Contains routines to least square match
    clouds of point pairs
*/

#ifndef VMLSTRANSFORM_H
#define VMLSTRANSFORM_H

#include "vectormath.h"
#include <vector>

//! Class encapsulating a pair of planar points together with a weight
/*! The function vmLeastSquareTransform computes a transformation that
  applied to the point \c _a of a set of point pairs minimizes the
  weighted sum of distance to the corresponding \c _b points.
*/
class VmPointPair2 
{
 public:
  //! First point
  VmVector2 _a;
  //! Second point
  VmVector2 _b;
  //! The weight, with which the squared distance is multiplied
  double _weight;

  //! Empty pair  
  VmPointPair2 ();  
  //! Pair of \c a and \c b with a given \c weight
  VmPointPair2 (const VmVector2& a, const VmVector2& b, double weight=1);
  //! Pair of \c (ax,ay) and \c (bx,by) with a given \c weight
  VmPointPair2 (double ax, double ay, double bx, double by, double weight=1);
};

//! A pair of points in 2D
typedef VmPointPair2 VmPointPair;

//! List of pairs of point
typedef std::vector<VmPointPair2> VmPointPair2List;

//! Class encapsulating a pair of 3D points together with a weight
/*! The function vmLeastSquareTransform computes a transformation that
  applied to the point \c _a of a set of point pairs minimizes the
  weighted sum of distance to the corresponding \c _b points.
*/
class VmPointPair3 
{
 public:
  //! First point
  VmVector3 _a;
  //! Second point
  VmVector3 _b;
  //! The weight, with which the squared distance is multiplied
  double _weight;
  
  //! Empty pair
  VmPointPair3 ();  
  //! Pair of two points with a given weight
  VmPointPair3 (const VmVector3& a, const VmVector3& b, double weight=1);  
};

//! A pair of points in 3D
typedef std::vector<VmPointPair3> VmPointPair3List;


//! Computes the best fitting rigid body transform for a set of 2D point pairs
/*! Determines for two corresponding sets of planar point \c p[i]._a
    and \c p[i]._b the optimal rigid body transform of \c p._a[i] to fit
    \c p._b[i]. It is optimal in the least square sense, minimizing 
    \f{equation}
        \sum_i p[i]._weight ((Rot(phi)*p._a[i]+d)-p._b[i])
    \f}
    In the case of underdetermined DOF the translation and / or rotation is 0.
    The return value is the ssd error.
*/
double vmLeastSquareTransform (const VmPointPair2List& p, VmVector2& d, double& phi);

//! Computes the best fitting rigid body transform for a set of 3D point pairs
/*! The rotation is returned as a homogenous transform \c T \c evRatio
    is the ratio between the 1. and 2. singular value. If it is, say
    <0.01 the problem is rank deficient (points on a line) or very ill
    conditioned.

    See D.W. Eggert, A. Lorusso, R.B. Fisher: Estimating 3-D rigid
    body transformations: a compurison of four major algorithms,
    Machine Vision and Applications (1997) 9: 272-290
*/
double vmLeastSquareTransform (const VmPointPair3List& p, VmMatrix4x4& T, double& evRatio);

#endif
