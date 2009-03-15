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
/**
* @file vmTest.cc
*
* Unit test functions for (some) \c vectormath.cc function
*
* @author <a href="mailto:udo.frese@dfki.de">Udo Frese</a>
*/

#define SEGMENTATION_FAULT_ASSERT
#include <vectormath/vectormath.h>
#include <vectormath/vmRandom.h>
#include <vectormath/vmGeometry.h>
#include <algorithm>
#include <assert.h>
using namespace std;

 
void testVmDistancePointToLine ()
{
  double eps =  1E-7;  
  VmRandom rnd;  
  int n = 100000;
  for (int i=0; i<n; i++) {
    // sample two random points
    VmVector2 a, b, d, y, p;
    a[0] = rnd.gauss (0,1);
    a[1] = rnd.gauss (0,1);
    d[0] = rnd.gauss (0,1);
    d[1] = rnd.gauss (0,1);
    vmAdd (b, a, d);
    double dL = sqrt(d[0]*d[0] + d[1]*d[1]);    

    // Test 'b' as closest point
    {      
      double alpha = rnd.random (-M_PI/2, M_PI/2);
      double length = rnd.random (0, 2);
      p[0] = b[0] + length/dL*cos(alpha)*d[0] + length/dL*sin(alpha)*d[1];
      p[1] = b[1] + length/dL*cos(alpha)*d[1] - length/dL*sin(alpha)*d[0];
      VmVector2 pClosest;    
      double lambda;    
      double dist = vmDistancePointToLine (lambda, pClosest, p, a, b);
      assert (pClosest[0]==b[0]);
      assert (pClosest[1]==b[1]);
      assert (lambda==1);
      assert (fabs(dist-length)<eps);
    }
    

    // Test 'a' as closest point
    {      
      double alpha = rnd.random (-M_PI/2, M_PI/2);
      double length = rnd.random (0, 2)/dL;
      p[0] = a[0] - length*cos(alpha)/dL*d[0] - length/dL*sin(alpha)*d[1];
      p[1] = a[1] - length*cos(alpha)/dL*d[1] + length/dL*sin(alpha)*d[0];
      VmVector2 pClosest;    
      double lambda;    
      double dist = vmDistancePointToLine (lambda, pClosest, p, a, b);
      assert (pClosest[0]==a[0]);
      assert (pClosest[1]==a[1]);
      assert (lambda==0);
      assert (fabs(dist-length)<eps);
    }
    

    // Test interior point as closest
    {      
      double lambdaTrue = rnd.random (0,1);
      double length = rnd.random (-2, 2);
      p[0] = a[0] + lambdaTrue*d[0] - length/dL*d[1];
      p[1] = a[1] + lambdaTrue*d[1] + length/dL*d[0];
      VmVector2 pClosest;    
      double lambda;    
      double dist = vmDistancePointToLine (lambda, pClosest, p, a, b);
      assert (fabs(pClosest[0]-(a[0]+lambda*d[0]))<eps);
      assert (fabs(pClosest[1]-(a[1]+lambda*d[1]))<eps);
      assert (fabs(lambda-lambdaTrue)<eps);
      assert (fabs(dist-fabs(length))<eps);
    }    
  }  
  
}


void testVmDistancePointToLine2 ()
{
  double eps =  1E-7;  
  VmRandom rnd;  
  int n = 100000;
  for (int i=0; i<n; i++) {
    // sample two random points
    VmVector2 a, b, p;
    a[0] = rnd.gauss ();
    a[1] = rnd.gauss ();
    b[0] = rnd.gauss ();
    b[1] = rnd.gauss ();
    p[0] = rnd.gauss ();
    p[1] = rnd.gauss ();
    
    double lambda;
    VmVector2 pClosest;    
    double dist = vmDistancePointToLine (lambda, pClosest, p, a, b);
    VmVector2 d;
    vmSub (d, pClosest, p);
    assert (fabs(vmLength(d)-dist)<eps);
    double scpA, scpB, scpP;
    scpA = d[0]*a[0] + d[1]*a[1];    
    scpB = d[0]*b[0] + d[1]*b[1];    
    scpP = d[0]*p[0] + d[1]*p[1];

    double dd;    
    if (scpA>scpB) swap(scpA, scpB);
    if (scpP<scpA) dd = scpA-scpP;
    else if (scpB<scpP) dd = scpP-scpB;
    else dd = 0;

    assert (fabs(dd-dist*vmLength(d))<eps);    
  }
}


double intervalDistance (double aLow, double aHigh, double bLow, double bHigh)
{
  if (aHigh<bLow) return bLow-aHigh;
  else if (bHigh<aLow) return aLow-bHigh;
  else return 0;  
}


void testVmDistanceLineToLine()
{  
  double eps =  1E-7;  
  VmRandom rnd;  
  int n = 100000;
  for (int i=0; i<n; i++) {
    // sample two random points
    VmVector2 a1, a2, b1, b2;
    a1[0] = rnd.gauss ();
    a1[1] = rnd.gauss ();
    a2[0] = rnd.gauss ();
    a2[1] = rnd.gauss ();
    b1[0] = rnd.gauss ();
    b1[1] = rnd.gauss ();
    b2[0] = rnd.gauss ();
    b2[1] = rnd.gauss ();

    double lambdaA, lambdaB;
    VmVector2 pCA, pCB;
    double dist = vmDistanceLineToLine (lambdaA, pCA, lambdaB, pCB,
					a1, a2, b1, b2);
    
    
    VmVector2 dir;
    vmSub (dir, pCA, pCB);
    double scpA1, scpA2, scpB1, scpB2;
    scpA1 = vmDot (a1, dir);
    scpA2 = vmDot (a2, dir);
    if (scpA1>scpA2) swap (scpA1, scpA2);
    
    scpB1 = vmDot (b1, dir);
    scpB2 = vmDot (b2, dir);
    if (scpB1>scpB2) swap (scpB1, scpB2);

    VmVector2 d;
    vmSub (d, pCA, pCB);    
    assert (fabs(vmLength(d)-dist)<eps);    
    
    assert (fabs((1-lambdaA)*a1[0]+lambdaA*a2[0]-pCA[0])<eps);
    assert (fabs((1-lambdaA)*a1[1]+lambdaA*a2[1]-pCA[1])<eps);

    assert (fabs((1-lambdaB)*b1[0]+lambdaB*b2[0]-pCB[0])<eps);
    assert (fabs((1-lambdaB)*b1[1]+lambdaB*b2[1]-pCB[1])<eps);

    double dd = intervalDistance (scpA1, scpA2, scpB1, scpB2);
    assert (fabs(dd - dist*vmLength(dir))<eps);    
  }

}


int main (int argc, char** argv)
{
  testVmDistancePointToLine ();
  testVmDistancePointToLine2();  
  testVmDistanceLineToLine ();  
  return 0;  
}
