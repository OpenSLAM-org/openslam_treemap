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

/*! \file tmSlamDriver3D.h
    \brief Implementation of \c TmSlamDriver3D
    \author Udo Frese

    Contains the implementation of \c TmSlamDriver3D a driver for 3D
    landmark based SLAM with 3D landmarks and 3D/6DOF poses, but without
    any connection between successive poses.
*/
#include "tmSlamDriver3D.h"

#include "vectormath/vmLSTransform.h"


TmSlamDriver3D::TmSlamDriver3D ()
  :TmTreemap (), isFirstPoseX (true)
{
  clear ();  
}

void TmSlamDriver3D::create (const VmMatrix4x4& initialPose, int nrOfMovesPerStep, int maxNrOfUnsuccessfulMoves, int nrOfNonlinearLeaves)
{
  TmTreemap::create (nrOfMovesPerStep, maxNrOfUnsuccessfulMoves);  
  vmCopy (initialPose, this->initialPose);
  isFirstPoseX = true;
  dummyPoseId = newFeatureBlock (6);  
  this->nrOfNonlinearLeaves = nrOfNonlinearLeaves;  
}


void TmSlamDriver3D::registerRandomVariable (RandomVariable* rv)
{
  assert (rv!=NULL && rv->userId>=0);
  assert (variablesByUserId.find (rv->userId)==variablesByUserId.end());  
  variablesByUserId[rv->userId] = rv;
  rv->tree = this;
  assert(rv->dimension>0);
  rv->featureId = newFeatureBlock (rv->dimension);
  variablesByFeatureId[rv->featureId] = rv;  
}

  
void TmSlamDriver3D::deleteFeature (TmFeatureId id)
{
  // Remove from \c variablesByFeatureId
  IntRVMap::iterator it = variablesByFeatureId.find (id);
  if (it!=variablesByFeatureId.end()) {  
    // When the first 1D feature of multi-D feature is deleted, the multi-D feature is deleted too
    RandomVariable* rv = (*it).second;
    assert (rv != NULL);  
    variablesByFeatureId.erase (it);  
    
    // Remove from \c variablesByUserId
    IntRVMap::iterator it = variablesByUserId.find (rv->userId);
    assert (it!=variablesByUserId.end());  
    variablesByUserId.erase (it);

    delete rv;    
  }
}


TmSlamDriver3D::~TmSlamDriver3D ()
{
  clear();
}

void TmSlamDriver3D::clear()
{
  TmTreemap::clear ();
  for (IntRVMap::iterator it=variablesByFeatureId.begin();it!=variablesByFeatureId.end();it++)
    delete ((*it).second);  
  variablesByFeatureId.clear();
  variablesByUserId.clear();
  recentlyObservedLandmarks.clear();  
  vmOne (initialPose);
  dummyPoseId = -1;  
  statistic.clear();  
}


void TmSlamDriver3D::step ()
{
  // We dont need to do anything since we dont use odometry
}


TmTreemap::SlamStatistic TmSlamDriver3D::slamStatistics () const
{
  return statistic;
}


int TmSlamDriver3D::memory () const
{
  int mem = TmTreemap::memory()+ sizeof(TmSlamDriver3D) - sizeof (TmTreemap);
  for (IntRVMap::const_iterator it = variablesByUserId.begin(); it!=variablesByUserId.end(); it++)
    mem += (*it).second->memory ();
  mem += nonlinearLeaf.size () * sizeof(int);  
  return mem;  
}


void TmSlamDriver3D::observe (const LandmarkObservationList& obs)
{
  statistic.p++;
  statistic.pMarginalized++;
  statistic.m += obs.size();  
  for (int i=0; i<obs.size(); i++)
    if (landmarkFromUserId (obs[i].landmarkId)==NULL) {      
      registerRandomVariable (new Landmark (obs[i].landmarkId));
      statistic.n++;
    }  

  NonlinearLeaf* leaf;
  if (isFirstPoseX) leaf = new NonlinearLeaf (this, obs, initialPose);  
  else leaf = new NonlinearLeaf (this, obs);  
  leaf->resetFlag (TmNode::CAN_BE_INTEGRATED);  

  VmMatrix4x4 pose;
  if (isFirstPoseX) vmCopy (initialPose, pose);
  else leaf->estimatePose (pose, false, recentlyObservedLandmarks);
  
  leaf->setInitialEstimate (pose);  
  // when introducing a new leaf, use the measurement but only those to
  // landmarks already observed before.
  leaf->linearize (pose, true);
  addNonlinearLeaf (leaf);
  nonlinearLeaf.push_back (leaf->index);

  isFirstPoseX = false;
  setRecentlyObservedLandmarks (obs);  

  while ((int) nonlinearLeaf.size()>nrOfNonlinearLeaves) {
    TmNode* n = getNode (nonlinearLeaf.front());
    n->setFlag (TmNode::CAN_BE_INTEGRATED);
    n->resetFlagUpToRoot (TmNode::IS_FEATURE_PASSED_VALID);    
    if (n->parent!=NULL) n->setToBeOptimizedUpToRoot ();
    nonlinearLeaf.pop_front();
  }
  updateFeaturePassed();
}


void TmSlamDriver3D::computeLinearEstimate ()
{
  TmTreemap::computeLinearEstimate ();  
  for (int i=0; i<(int) nonlinearLeaf.size(); i++) {
    NonlinearLeaf* n = dynamic_cast<NonlinearLeaf*> (getNode (nonlinearLeaf[i]));
    if (!n->isFlag(TmNode::DONT_UPDATE_ESTIMATE)) n->hasBeenIncorporated = true;    
  }  
}


void TmSlamDriver3D::computeNonlinearEstimate ()
{
  for (int i=0; i<(int) nonlinearLeaf.size(); i++) {
    NonlinearLeaf* n = dynamic_cast<NonlinearLeaf*> (getNode (nonlinearLeaf[i]));
    if (n->hasBeenIncorporated) {
      VmMatrix4x4 pose;
      n->estimatePose (pose);      
      n->linearize (pose, false); // relinearize using the estimate    
    }    
  }
  computeLinearEstimate ();
}


void TmSlamDriver3D::setRecentlyObservedLandmarks (const LandmarkObservationList& obs)
{
  recentlyObservedLandmarks.clear();
  recentlyObservedLandmarks.reserve (obs.size());
  for (int i=0; i<obs.size(); i++) recentlyObservedLandmarks.push_back (obs[i].landmarkId);  
}


//**************** TmSlamDriver3D::RandomVariable

int TmSlamDriver3D::RandomVariable::memory () const
{
  return sizeof (TmSlamDriver3D);
}

TmSlamDriver3D::RandomVariable::~RandomVariable ()
{}





//**************** TmSlamDriver3D::Landmark ************

TmSlamDriver3D::Landmark::~Landmark ()
{}


void TmSlamDriver3D::Landmark::getEstimate (VmVector3& p) const
{
  p[0] = tree->feature[featureId  ].est;
  p[1] = tree->feature[featureId+1].est;
  p[2] = tree->feature[featureId+2].est;  
}


void TmSlamDriver3D::Landmark::setInitialEstimate (const VmVector3& p)
{
  tree->setInitialEstimate (featureId  , p[0]);
  tree->setInitialEstimate (featureId+1, p[1]);
  tree->setInitialEstimate (featureId+2, p[2]);  
}

int TmSlamDriver3D::Landmark::memory () const
{
  return TmSlamDriver3D::RandomVariable::memory() - sizeof (RandomVariable) + sizeof (Landmark);
}


//*************** TmSlamDriver3D::NonlinearLeaf *************

TmSlamDriver3D::NonlinearLeaf::NonlinearLeaf (TmSlamDriver3D* tree)
  :TmNode (), hasBeenIncorporated(false), hasKnownPose (false)
{
  this->tree = tree;  
  vmOne (knownPose);  
}


TmSlamDriver3D::NonlinearLeaf::NonlinearLeaf (TmSlamDriver3D* tree, const LandmarkObservationList& obs)
  :TmNode (), obs(obs), hasBeenIncorporated(false), hasKnownPose (false)
{
  this->tree = tree;
  vmOne (knownPose);  
}


TmSlamDriver3D::NonlinearLeaf::NonlinearLeaf (TmSlamDriver3D* tree, const LandmarkObservationList& obs, const VmMatrix4x4& knownPose)
  :TmNode (), obs(obs), hasBeenIncorporated(false), hasKnownPose (true)
{
  this->tree = tree;
  vmCopy (knownPose, this->knownPose);  
}


TmNode* TmSlamDriver3D::NonlinearLeaf::duplicate() const
{
  return new NonlinearLeaf (*this);  
}


int TmSlamDriver3D::NonlinearLeaf::memory() const
{
  return TmNode::memory() - sizeof(TmNode) + sizeof (NonlinearLeaf) +
    obs.capacity() * sizeof (LandmarkObservation);
}



void TmSlamDriver3D::NonlinearLeaf::setInitialEstimate (const VmMatrix4x4& pose)
{
  for (int i=0; i<obs.size(); i++) {
    const LandmarkObservation& o = obs[i];
    assert (o.landmarkId>=0);
    Landmark* rv = ((TmSlamDriver3D*)tree)->landmarkFromUserId (o.landmarkId);    
    if (rv!=NULL) {
      VmVector3 pW;
      vmApplyTransformToPoint (pW, pose, o.pos);
      rv->setInitialEstimate (pW); // Only done if there is no estimate yet      
    }
  }  
}


void TmSlamDriver3D::NonlinearLeaf::estimatePose (VmMatrix4x4& pose, bool allLandmarks, const XycVector<int>& landmarks)
{
  assert (!allLandmarks || landmarks.empty());  
  VmPointPair3List pl;
  for (int i=0; i<obs.size(); i++) {
    const LandmarkObservation& o = obs[i];
    int j;    
    for (j=0; j<landmarks.size(); j++) if (o.landmarkId==landmarks[j]) break;
    if (allLandmarks || j<landmarks.size()) { // is overlapping
      VmVector3 pW;
      Landmark* lm = ((TmSlamDriver3D*)tree)->landmarkFromUserId (o.landmarkId);   
      lm->getEstimate (pW);
      if (finite(pW[0])) pl.push_back (VmPointPair3 (o.pos, pW, 1/vmTrace(o.posCov)));      
    }
  }
  double evRatio;  
  vmLeastSquareTransform (pl, pose, evRatio);
  if (evRatio<0.005) {
    char txt[200];
    sprintf (txt, "Overlapping landmarks form a rank deficient constellation (%f)", evRatio);    
    throw runtime_error (txt);
  }  
}


void TmSlamDriver3D::NonlinearLeaf::linearize (VmMatrix4x4& pose, bool useMeasurement)
{
  TmSlamDriver3D* tree = dynamic_cast<TmSlamDriver3D*> (this->tree);  
  if (index!=-1) beforeChange ();

  TmExtendedFeatureList fl;
  if (!hasKnownPose) 
    for (int i=0; i<6; i++) 
      fl.push_back (TmExtendedFeatureId (tree->dummyPoseId+i, 1)); // first 6 DOF are pose
  int nrOfRows = 0; // no row for pose
  for (int i=0; i<obs.size(); i++) {
    const LandmarkObservation& o = obs[i];
    Landmark* lm = tree->landmarkFromUserId (o.landmarkId);
    fl.push_back (TmExtendedFeatureId (lm->featureId  , 1));
    fl.push_back (TmExtendedFeatureId (lm->featureId+1, 1));
    fl.push_back (TmExtendedFeatureId (lm->featureId+2, 1));
    nrOfRows+=3;    
  }
//  for (int i=0; i<(int) fl.size(); i++) assert (!tree->feature[fl[i].id].isEmpty());

  // Construct gaussian
  TmGaussian G; // first make a Gaussian that represents the pose, than marginalize the pose out
  G.create (fl, nrOfRows); 
  int jLm, jPose; // column of the pose and current landmark  
  if (hasKnownPose) {
    jLm = 0;
    jPose = -1;
  }
  else {
    jLm = 6;
    jPose = 0;
  }  
  for (int i=0; i<obs.size(); i++) {
    const LandmarkObservation& o = obs[i];
    VmVector3 pW;
    if (useMeasurement) {
      vmApplyTransformToPoint (pW, pose, o.pos);      
    }
    else {
      Landmark* lm = tree->landmarkFromUserId (o.landmarkId);
      lm->getEstimate (pW);
    }
    o.addToGaussian (G, jPose, jLm, pW, pose);    
    jLm+=3;    
  }

  if (!hasKnownPose) { // Normal case: The pose is unknown and will be marginalized out
    G.triangularize ();
    G.computeMarginal (gaussian, 6);    
  }
  else gaussian = G; // first pose: We know the pose exactly and implicitly condition on it

  if (index!=-1) afterChange ();
}


void TmSlamDriver3D::NonlinearLeaf::linearize ()
{
  if (hasKnownPose) linearize (knownPose, false);
  else {    
    VmMatrix4x4 pose;
    estimatePose (pose);
    linearize (pose, false);
  }  
}

//************** TmSlamDriver3D::LandmarkObservation ************

void TmSlamDriver3D::LandmarkObservation::jacobian (double JPL[3][3], double JPP[3][6],
                                                    const VmVector3& pL0, VmMatrix4x4& T0)
{
  // f = (T0*T(pP))^{-1}*pL
  // f linearized = T(pP)^{-1}*T0^{-1}*pL

  // dF/dPL = T0^{-1}
  for (int i=0;i<3;i++) for (int j=0;j<3;j++) JPL[i][j] = T0[j][i];  

  // dF/dPP = M |T0^{-1}*pL
  VmVector3 v;
  vmApplyInverseTransformToPoint (v, T0, pL0);  
  double M[3][6] = {
    {-1, 0, 0, -v[1],  v[2],     0},
    {0, -1, 0,  v[0],     0, -v[2]},
    {0, 0, -1,     0, -v[0],  v[1]}
  };
  for (int i=0;i<3;i++) for (int j=0; j<6; j++) JPP[i][j] = M[i][j];
}

void TmSlamDriver3D::LandmarkObservation::addToGaussian (TmGaussian& G, int jPose, int jLm, const VmVector3& pL0, VmMatrix4x4& T0) const
{
  double JPL[3][3], JPP[3][6];
  jacobian (JPL, JPP, pL0, T0);
  VmMatrix3x3 L;
  vmCholeskyInverse (posCov, L);

  // Our equation is: L*(p0+JPL*(pL-pL0)+JPP*(pP)-pos)=0 with pL and pP being estimated
  // or  L*JPL*pL + L*JPP*pP + L*(p0-JPL*pL0-pos)

  VmVector3 b, hv, lB;
  vmApplyInverseTransformToPoint (b, T0, pL0);
  vmSub (b, b, pos);  
  vmMultiply (hv, JPL, pL0);
  vmSub (b, b, hv);
  vmMultiply (lB, L, b);
 
  int i = G.R.rows();  
  G.R.appendRow (3, true);
  if (jPose>=0) for (int ii=0;ii<3;ii++) for (int jj=0; jj<6; jj++) {
    double sum = L[ii][0]*JPP[0][jj] + L[ii][1]*JPP[1][jj] + L[ii][2]*JPP[2][jj];
    G.R (i+ii, jPose+jj) = sum;
  }  
  if (jLm>=0) for (int ii=0;ii<3;ii++) for (int jj=0; jj<3; jj++) {
    double sum = L[ii][0]*JPL[0][jj] + L[ii][1]*JPL[1][jj] + L[ii][2]*JPL[2][jj];
    G.R (i+ii, jLm+jj) = sum;
  }  
  G.R.storeCol (lB, i, G.feature.size());
}
                                                          
