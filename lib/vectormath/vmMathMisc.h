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
/*! \file vmMathMisc.h
    \author Udo Frese
    \brief Collection of mathematical routines.

    A collection of mathematical routines. This routines are stand-alone.
    If you need only one, you can copy it as well instead of including
    vmMathMisc.h

*/

#ifndef VMMATHMISC_H
#define VMMATHMISC_H

#include <limits>
#include <math.h>

/*! Usually the matrices consist of \c double.
     If \c USEFLOAT is defined, the matrices consist of \c float.
*/
#ifndef USEFLOAT
typedef double VmReal;
#else
typedef float VmReal;
#endif


//! returns v*v
inline double vmSqr (double v) 
{
  return v*v;
}


//! Returns IEEE 854 infinity
inline VmReal vmInf() 
{
  return std::numeric_limits<VmReal>::infinity();
}

//! Returns IEEE 854 NaN
inline VmReal vmNaN()
{
  return std::numeric_limits<float>::quiet_NaN(); // doesn't work with double
}

//! Swaps a and b
inline void vmSwap (VmReal& a, VmReal& b)
{
    VmReal park = a;
    a = b;
    b = park;
}



//! Converts radian to degree
inline VmReal vmDeg(VmReal angleRad) 
{
  return 180*angleRad/M_PI;
}

//! Converts degree to radian
inline VmReal vmRad(VmReal angleDeg)  
{
  return angleDeg/180*M_PI;
}


//! Solves a quadratic equation \f$ ax^2+bx+c=0 \f$
/*!If there are solutions, the smaller is returned in \c xLow, the larger in 
   \c xHigh and the return value is \c true. Otherwise \c false and \c NaN are returned.
*/
inline bool vmSolvePolynom2 (VmReal a, VmReal b, VmReal c, VmReal& xLow, VmReal& xHigh)
{
    // Numerically more precise formula (reimplemented after Numerical Recipes 5.6)
    double r, q;
    r = b*b-4*a*c;
    if (r<0) {
        xLow = xHigh = vmNaN();
        return false;
    }
    if (b>0) {
        q = (-b-sqrt(r))/2;
        xLow  = q/a;
        xHigh = c/q;
    }
    else {
        q = (-b+sqrt(r))/2;
        xHigh = q/a;
        xLow  = c/q;
    }
    return true;
}

//! Solves a quadratic equation \f$ a[2]x^2+a[1]x+a[0]=0 \f$
/*! 
   \warning In this version \c a[0] is the coefficient for \f$ x^0 \f$ and
   \c a[2] for \f$ x^2 \f$.

   \sa vmSolvePolynom2
*/
inline bool vmSolvePolynom2 (VmReal a[3], VmReal& xLow, VmReal& xHigh)
{
    return vmSolvePolynom2 (a[2], a[1], a[0], xLow, xHigh);
}


//! Returns the maximum of \f$ a[0] + a[1]*x + a[2]*x^2 \f$ in \c [xLow..xHigh] 
/*! The maximum itself is returned in \c y and the corresponding \c x value 
    as function result.
*/
inline VmReal vmP2Maximum (VmReal a[2], VmReal& y, VmReal xLow=-0.5, VmReal xHigh=+0.5)
{
    VmReal x;
    if (a[2]<0) {
        // f''<0
        VmReal xMax = -a[1]/(2*a[2]);
        if (xMax>xHigh) x=xHigh;
        else if (xMax<xLow) x=xLow;
        else x=xMax;
    }
    else if (a[2]>0) {
        VmReal xMin = -a[1]/(2*a[2]);
        if (xMin>(xLow+xHigh)/2) x=xLow;
        else x=xHigh;
    }
    else {
        if (a[1]>=0) x=xHigh;
        else x=xLow;
    }
    y = a[0]+x*(a[1]+x*a[2]);
    return x;
}


//! Normalize \c angle into the interval \c [reference-M_PI..reference+M_PI)
inline double vmNormalizedAngle (double angle, double reference=0);

//! Normalize \c x into the interval \c [x-ref/2..x+ref/2]
/*! \c periodicity is the periodicity assumed. */
inline double vmNormalized (double x, double ref, double periodicity);

//! Returns, whether \c angle is inside the interval \c [low..high] modulo \c 2*M_PI
/*! That means, whether there exists an integer 'i', so that \c
    low<=angle+i*2*M_PI<=high
*/
inline bool vmInPeriodicInterval (double angle, double low, double high);

//! Returns, whether the intervals \c [low1..high1] and \c [low2..high2] intersect modulo \c 2*M_PI,
/*! That means, whether there exist VmReal a,b and integer i, so that
    \c a=b+2*M_PI*i and \c low1<=a<=high1 and \c low2<=b<=high2
    . Setting \c low1==high1 can be used to test, whether an angle is
    contained in a certain angle interval.
*/
inline bool vmPeriodicIntervalsIntersect (double low1, double high1, double low2, double high2);

//! Computes the intersection between two periodic intervals. 
/*! Formally this means the set of angles \c alpha, which are both
    inside \c [low1..high1] and \c [low2..high2] modulo \c 2*M_PI as
    in the definition of \c inPeriodicInterval. If both intervals
    together have at most a length of 2*M_PI, the resulting
    intersection is again a periodic interval (\c component==1) or
    empty (\c components==0). Otherwise, it may happen that the
    intersection consists of two intervals (\c components==2). The the
    larger one is returned. The return value specifies, whether the
    intersection is not empty (equivalent to \c components>0).
*/
inline bool vmPeriodicIntervalsIntersection (double& iLow, double& iHigh, int& components,
					     double low1, double high1, double low2, double high2);

#include "vmMathMisc_i.h"

#endif
