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


/*!\file tmSlamDriver2DP.cc 
   \brief Implementation of \c TmSlamDriver2DP
   \author Udo Frese

   Contains the implementation of class \c
   TmSlamDriver2DP representing treemap driver for 2D pose relation based
   SLAM.
*/

#include "tmSlamDriver2DP.h"
#include <xymatrix/xymMatrixC.h>
#include <xymatrix/xymOperations.h>
#include <assert.h>

TmSlamDriver2DP::TmSlamDriver2DP ()
 :TmTreemap(), pose2Feature(), p(0)
{}


void TmSlamDriver2DP::clear ()
{
  TmTreemap::clear();
  pose2Feature.clear();
  p = 0;  
}


void TmSlamDriver2DP::create (int nrOfPosesReserved, int nrOfMovesPerStep, int maxNrOfUnsuccessfulMoves)
{
  TmTreemap::create (nrOfMovesPerStep, maxNrOfUnsuccessfulMoves);
  feature.resize (nrOfPosesReserved*3);  
  pose2Feature.clear();
  pose2Feature.reserve (nrOfPosesReserved);  
  p = 0;  
}


void TmSlamDriver2DP::allocatePose (int id)
{
  assert (id>=0);
  if (id>=pose2Feature.size()) pose2Feature.resize (id+1, -1);
  else if (pose2Feature[id]>=0) return;  
  int pid = pose2Feature[id] = newFeatureBlock (3);
}


void TmSlamDriver2DP::setInitialEstimate (const Link& link)
{
  allocatePose (link.poseA);
  if (link.poseB>=0) allocatePose (link.poseB);  
  int pfidA = pose2Feature [link.poseA];
  if (!isfinite(feature[pfidA].est)) {
    p++;
    VmVector3 pBEst;
    poseEstimate (link.poseB, pBEst);
    assert (finite(pBEst[0])); // otherwise B was undefined too
    
    double c = cos(pBEst[2]), s = sin(pBEst[2]);    
    TmTreemap::setInitialEstimate (pfidA  , pBEst[0] + c*link.d[0] - s*link.d[1]);
    TmTreemap::setInitialEstimate (pfidA+1, pBEst[1] + s*link.d[0] + c*link.d[1]);
    TmTreemap::setInitialEstimate (pfidA+2, pBEst[2] + link.d[2]);
  }
  else if (link.poseB>=0) {
    int pfidB = pose2Feature[link.poseB];
    if (!isfinite(feature[pfidB].est)) {
      p++;
      VmVector3 pAEst;
      poseEstimate (link.poseA, pAEst);
      assert (finite (pAEst[0])); // otherwise B was undefined too
    
      double bTheta = pAEst[2] - link.d[2], c = cos(bTheta), s = sin(bTheta);    
      TmTreemap::setInitialEstimate (pfidB  , pAEst[0] - c*link.d[0] + s*link.d[1]);
      TmTreemap::setInitialEstimate (pfidB+1, pAEst[1] - s*link.d[0] - c*link.d[1]);
      TmTreemap::setInitialEstimate (pfidB+2, bTheta);
    }
  }  
}


void TmSlamDriver2DP::addLink (const Link& link)
{
  setInitialEstimate (link);  

  // Create the list of features involved. The order of the features
  // corresponds to the columns of the SRI matrix passed.
  TmExtendedFeatureList fl;
  int poseAFeature = pose2Feature [link.poseA]; 
  fl.push_back (TmExtendedFeatureId (poseAFeature  )); // a_x
  fl.push_back (TmExtendedFeatureId (poseAFeature+1)); // a_y 
  fl.push_back (TmExtendedFeatureId (poseAFeature+2)); // a_theta
  if (link.poseB>=0) {    
    int poseBFeature = pose2Feature[link.poseB];  
    fl.push_back (TmExtendedFeatureId (poseBFeature  )); // b_x
    fl.push_back (TmExtendedFeatureId (poseBFeature+1)); // b_y   
    fl.push_back (TmExtendedFeatureId (poseBFeature+2)); // b_theta
  }  

  // Compute the linearization point
  // We use the existing estimate for b (in particular b[2])
  // and choose a not from the initial estimate but according
  // to the measurement itself. That's better when closing loops
  // for the first iteration of least square.
  VmVector3 poseALinPoint, poseAEst, poseBLinPoint;
  poseEstimate (link.poseB, poseBLinPoint);
  poseEstimate (link.poseA, poseAEst);  
  double c = cos(poseBLinPoint[2]), s = sin(poseBLinPoint[2]);  
  poseALinPoint[0] = poseBLinPoint[0] + c*link.d[0] - s*link.d[1];
  poseALinPoint[1] = poseBLinPoint[1] + s*link.d[0] + c*link.d[1];
  poseALinPoint[2] = vmNormalizedAngle (poseBLinPoint[2] + link.d[2], poseAEst[2]);

  // Now linearize the measurement equation and scale according to the covariance
  XymMatrixC LInv, A;
  linearizeLink (A, poseALinPoint, poseBLinPoint, link.d, link.poseB>=0);
  inverseCholeskyFactor (LInv, XymMatrixC(link.dCov));  

  // Add the link as a single Gaussian leaf
  TmGaussian linearizedLink (LInv*A, fl, false);
  addLeaf (linearizedLink);
}


void TmSlamDriver2DP::linearizeLink (XymMatrixC& A, const VmVector3& poseALinPoint, 
				     const VmVector3& poseBLinPoint, const VmVector3& d, bool includeJacobianForB)
{
  assert (isFinite (poseALinPoint) && isFinite(poseBLinPoint));  
  if (includeJacobianForB) A.create (3, 7);
  else A.create (3, 4);  
  double c = cos(poseBLinPoint[2]), s = sin(poseBLinPoint[2]);  
  VmMatrix3x3 Ja = {{ c,  s, 0},
                    {-s,  c, 0},
		    { 0,  0, 1}}; // df(a,b)/da | a=poseALinPoint, b=poseBLinPoint
  VmMatrix3x3 Jb = {{-c, -s, -s*(poseALinPoint[0]-poseBLinPoint[0]) + c*(poseALinPoint[1]-poseBLinPoint[1])},
		    { s, -c, -c*(poseALinPoint[0]-poseBLinPoint[0]) - s*(poseALinPoint[1]-poseBLinPoint[1])},
                    { 0,  0, -1}}; // df(a,b)/db | a=poseALinPoint, b=poseBLinPoint
  VmVector3 y = { c*(poseALinPoint[0]-poseBLinPoint[0]) + s*(poseALinPoint[1]-poseBLinPoint[1]) -d[0],
		 -s*(poseALinPoint[0]-poseBLinPoint[0]) + c*(poseALinPoint[1]-poseBLinPoint[1]) -d[1],
		  vmNormalizedAngle((poseALinPoint[2]-poseBLinPoint[2])- d[2])};   // f(poseALinPoint, poseBLinPoint) - d  
  
  vmMultiplySub (y, Ja, poseALinPoint);
  if (includeJacobianForB) vmMultiplySub (y, Jb, poseBLinPoint);
  A.store (Ja, 0, 0);
  if (includeJacobianForB){
    A.store (Jb, 0, 3);
    A.storeCol (y, 0, 6);
  }
  else A.storeCol (y, 0, 3);
}



void TmSlamDriver2DP::poseEstimate (int idx, double& x, double& y, double &theta) const
{
  if (idx>=0) {    
    int poseFeat = pose2Feature [idx];  
    x     = feature[poseFeat  ].est;
    y     = feature[poseFeat+1].est;
    theta = feature[poseFeat+2].est;
  }
  else x = y = theta = 0;  
}

void TmSlamDriver2DP::nameOfFeature (char* txt, int featureId, int& n) const
{
  for (int i=0; i<pose2Feature.size(); i++) if (pose2Feature[i]==featureId) {    
    sprintf(txt, "p%d", i);
    n=3;  
    return;    
  }
  strcpy (txt, "??");
  n = 1;  
}


int TmSlamDriver2DP::memory () const
{
  return TmTreemap::memory()+ sizeof(TmSlamDriver2DP) - sizeof (TmTreemap) + pose2Feature.capacity()*sizeof(int);
}

