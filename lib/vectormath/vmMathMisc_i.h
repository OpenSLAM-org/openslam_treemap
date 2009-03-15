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
//! \file vmMathMisc_i.h Implementation of \c vmMathMisc.h
/*! \file vmMathMisc_i.h
    \author Udo Frese

    Contains the implementation of inline function in \c vmMathMisc.h
*/
inline double vmNormalizedAngle (double angle, double reference)
{
    if (!isfinite(angle) || !isfinite(reference)) return vmNaN();
    angle -=reference;
    angle -=floor(angle/(2*M_PI)+0.5)*2*M_PI;
    return angle+reference;
}



inline double vmNormalized (double x, double reference, double periodicity)
{
    if (!isfinite(x) || !isfinite(reference)) return vmNaN();
    x -=reference;
    x -=floor(x/periodicity+0.5)*periodicity;
    return x+reference;
}



inline bool vmInPeriodicInterval (double angle, double low, double high)
{
    int i = (int) ceil ((low-angle)/(2*M_PI));
    return angle+i*(2*M_PI)<=high;
}



inline bool vmPeriodicIntervalsIntersect (double low1, double high1, double low2, double high2)
{
    double dLow, dHigh; // Compute the interval [low1..high1]-[low2..high2]
    dLow  = low1-high2;
    dHigh = high1-low2;
    // Look, whether a multiple of 2*M_PI exists inside [dLow..dHigh]
    double norm = floor(dLow/(2*M_PI)+0.5)*2*M_PI;
    dLow -= norm;
    dHigh -= norm;
    // dLow is in [-M_PI..+M_PI]
    if (dLow<=0 && dHigh>=0) return true;
    else if (dLow<=2*M_PI && dHigh>=2*M_PI) return true;
    else return false;
}



inline bool vmPeriodicIntervalsIntersection (double& iLow, double& iHigh, int& components,
				      double low1, double high1, double low2, double high2)
{
    iLow = 0;
    iHigh = -1;
    components = 0;

    // Special cases
    if (low1>high1 || low2>high2) return false;

    if (high1>=low1+2*M_PI) {
        iLow = low2;
        iHigh = high2;
        components = 1;
        return true;
    }

    if (high2>=low2+2*M_PI) {
        iLow = low1;
        iHigh = high1;
        components = 1;
        return true;
    }
    

    // Normalize low2 to just below low1
    double norm = floor ((low1-low2)/(2*M_PI))*(2*M_PI);
    low2  += norm;
    high2 += norm;
    if (high2>=low1) {
        // Intersection
        iLow = low1;
        if (high1<high2) iHigh = high1;
        else iHigh = high2;
        components++;
    }

    // Normalize low1 to just below low2
    norm = floor ((low2-low1)/(2*M_PI))*(2*M_PI);
    low1  += norm;
    high1 += norm;
    if (high1>=low2) {
        // Intersection
        double i2Low, i2High;
        i2Low = low2;
        if (high2<high1) i2High = high2;
        else i2High = high1;
        if (components==0 || i2High-i2Low>iHigh-iLow) {
            iLow = i2Low;
            iHigh = i2High;
        }
        components++;
    }
    
    return components>0;
}
