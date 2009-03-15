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
/*!\file tmSlamDriver2DL.cc 
   \brief Implementation of class \c TmSlamDriver2DL
   \author Udo Frese

   Contains the implementation of class \c
   TmSlamDriver2DL representing treemap driver for 2D landmark based
   SLAM.
*/

#include "tmSlamDriver2DL.h"

TmSlamDriver2DL::TmSlamDriver2DL ()
  :TmTreemap(), pose(), landmark(), level(),
   statistic(), nrOfLandmarks(0), 
   landmarkBaseFeature(-1), poseFeature(-1), poseIndex(-1), 
   sparsificationDistance (vmInf()), 
   nrOfNonlinearLeaves (-1), nonlinearLeaf(), estimateIncludes(-1),
   currentLevel (0)
{
  clearRelativePose ();  
}


void TmSlamDriver2DL::clear()
{
  TmTreemap::clear();  
  clearRelativePose ();
  currentLevel = 0;
  pose.clear();
  pose.compact();
  landmark.clear();  
  landmark.compact();  
  level.clear();
  level.compact();  
  statistic.clear();  
}


TmTreemap::SlamStatistic TmSlamDriver2DL::slamStatistics () const
{
  return statistic;
}

void TmSlamDriver2DL::clearRelativePose ()
{  
  relativePose[0] = relativePose[1] = relativePose[2] = 0;
  relativePoseCov[0][0] = relativePoseCov[0][1] = relativePoseCov[0][2] = 0;
  relativePoseCov[1][0] = relativePoseCov[1][1] = relativePoseCov[1][2] = 0;
  relativePoseCov[2][0] = relativePoseCov[2][1] = relativePoseCov[2][2] = 0;
  relativeDistance = 0;  
}


void TmSlamDriver2DL::create (int nrOfLandmarks, int reservePoses, const VmVector3& initialPose, const VmMatrix3x3& initialPoseCov,
                              int nrOfMovesPerStep, int maxNrOfUnsuccessfulMoves, int nrOfNonlinearLeaves, double sparsificationDistance)
{
  clear();
  TmTreemap::create (nrOfMovesPerStep, maxNrOfUnsuccessfulMoves);  
  assert (nrOfNonlinearLeaves>=1);  
  this->nrOfNonlinearLeaves    = nrOfNonlinearLeaves;  
  this->sparsificationDistance = sparsificationDistance;  
  feature.reserve (2*nrOfLandmarks+3*reservePoses);  
  this->nrOfLandmarks = nrOfLandmarks;
  landmarkBaseFeature = newFeatureBlock (2*nrOfLandmarks);
  for (int i=0; i<nrOfLandmarks; i++) setFlagsForLandmark (landmarkFeature(i), 0);
  landmark.resize (nrOfLandmarks);  
  clearRelativePose ();
  poseBaseFeature = feature.size();  // Is this clean?
  poseFeature = newPose (0);
  setInitialEstimateForPose (poseFeature, initialPose);
  nonlinearLeaf.clear();  
  estimateIncludes = -1;  
  statistic.clear();
  statistic.p = 1;

  TmFeatureList fl;
  fl.push_back (poseFeature);
  fl.push_back (poseFeature+1);
  fl.push_back (poseFeature+2);

  NonlinearLeaf* leaf = new NonlinearLeaf (this);
  leaf->resetFlag (TmNode::CAN_BE_INTEGRATED);  
  leaf->poseFeature = poseFeature;  
  leaf->absolutePose = AbsolutePose (poseFeature, initialPose, initialPoseCov);
  leaf->linearize (true);
  addNonlinearLeaf (leaf);
  nonlinearLeaf.push_back (leaf->index);
}


void TmSlamDriver2DL::step (const VmVector3& odometry, const VmMatrix3x3& odometryCov)
{
  VmVector3 newRelativePose;  
  double theta = relativePose[2], c=cos(theta), s=sin(theta);
  newRelativePose[0] = relativePose[0] + c*odometry[0] - s*odometry[1];
  newRelativePose[1] = relativePose[1] + s*odometry[0] + c*odometry[1];
  newRelativePose[2] = relativePose[2] + odometry[2];
  
  VmMatrix3x3 J1 = 
    {{1, 0, -s*odometry[0] - c*odometry[1]},
     {0, 1,  c*odometry[0] - s*odometry[1]},
     {0, 0,  1
     }};
  VmMatrix3x3  J2 = 
    {{c, -s, 0},
     {s, c, 0},
     {0, 0, 1}};
     
  // relativePoseCov = J1*relativePoseCov*J1^T + J2*odometry*J2^T
  VmMatrix3x3 newCov;
  vmJCKt (newCov, J1, relativePoseCov, J1);
  vmJCKtAdd (newCov, J2, odometryCov, J2);
  vmCopy (newCov, relativePoseCov);
  vmCopy (newRelativePose, relativePose);  

  relativeDistance += sqrt(odometry[0]*odometry[0] + odometry[1]*odometry[1]);  
}


void TmSlamDriver2DL::setLevel (int level)
{
  assert(0<=level);  
  this->currentLevel = level;
  if (level>=this->level.size()) this->level.resize(level+1);  
}


int TmSlamDriver2DL::sharedLandmarks (const ObservationList& a, const ObservationList& b) 
{
  int ctr=0;
  for (int i=0; i<(int) a.size(); i++)
    for (int j=0; j<(int) b.size(); j++) 
      if (a[i].id==b[j].id) {
        ctr++;
        break;
      }
  return ctr;  
}


void TmSlamDriver2DL::setFlagsForPose (int fId, int flags)
{
  int whichFlags = TmFeature::CAN_BE_MARGINALIZED_OUT | TmFeature::CAN_BE_SPARSIFIED | TmFeature::USER_FLAGS;  
  feature[fId  ].setFlag (whichFlags, POSEX     | flags);
  feature[fId+1].setFlag (whichFlags, POSEY     | flags);
  feature[fId+2].setFlag (whichFlags, POSETHETA | flags);
}


void TmSlamDriver2DL::setFlagsForLandmark (int fId, int flags)
{
  int whichFlags = TmFeature::CAN_BE_MARGINALIZED_OUT | TmFeature::CAN_BE_SPARSIFIED | TmFeature::USER_FLAGS;  
  feature[fId  ].setFlag (whichFlags, LANDMARKX | flags);
  feature[fId+1].setFlag (whichFlags, LANDMARKY | flags);
}


void TmSlamDriver2DL::initLandmark (int id)
{
  if (id>=landmark.size()) landmark.resize (id+1);
  Landmark& lm = landmark[id];
  assert (lm.featureId==-1);  
  lm.featureId = landmarkFeature (id);
  if (lm.level==-1) {    
    lm.level                   = currentLevel;  
    lm.nextLandmarkInSameLevel = level[currentLevel].firstLandmarkInLevel;    
    level[currentLevel].firstLandmarkInLevel = id;
  }  
  setFlagsForLandmark (lm.featureId, 0);  
}


void TmSlamDriver2DL::setInitialEstimateForPose (int fId, const VmVector3& pose)
{
  setInitialEstimate (fId  , pose[0]);
  setInitialEstimate (fId+1, pose[1]);
  setInitialEstimate (fId+2, pose[2]);  
}


int TmSlamDriver2DL::newPose (double distance)
{
  int id = newFeatureBlock (3);

  int poseId = (id-poseBaseFeature)/3;
  assert (0<=poseId && poseId<=(int) pose.size());
  if (poseId==(int) pose.size()) pose.resize (poseId+1);

  Pose& p = pose[poseId];
  p.featureId = id;
  p.nextPoseId = -1;
  p.prevPoseId = poseIndex;  
  p.level      = currentLevel;  
  if (currentLevel>=level.size()) level.resize (currentLevel+1);
  p.nextPoseInSameLevel = level[currentLevel].firstPoseInLevel;  
  level[currentLevel].firstPoseInLevel = poseId;  
  if (poseIndex>=0) {    
    Pose& oldPose = pose[poseIndex];  
    p.distance = oldPose.distance + distance;
    p.poseNr   = oldPose.poseNr+1;
    oldPose.nextPoseId = poseId;
  }  
  else {
    p.distance = 0;
    p.poseNr   = 0;
  }  
  poseFeature = id;  
  poseIndex   = poseId;  
  int flag = TmFeature::CAN_BE_MARGINALIZED_OUT;
  setFlagsForPose (id, flag);  
  return id;  
}


void TmSlamDriver2DL::deleteFeature (TmFeatureId id)
{
  if (id>=poseBaseFeature) {
    int poseId = (id-poseBaseFeature)/3;
    Pose& p = pose[poseId];
    if (p.poseNr>=0) {
      statistic.pMarginalized++;      
      p.poseNr    = -1;
      p.featureId = -1;
      p.distance  = 0;
      int prevLevelPoseId;
      if (p.prevPoseId>=0 && pose[p.prevPoseId].nextPoseInSameLevel == poseId) 
        prevLevelPoseId = p.prevPoseId;
      else {
        int i=level[p.level].firstPoseInLevel;
        prevLevelPoseId = -1;        
        while (i!=poseId) {
          prevLevelPoseId = i;          
          i = pose[i].nextPoseInSameLevel;
        }        
      }
      if (prevLevelPoseId>=0) pose[prevLevelPoseId].nextPoseInSameLevel = p.nextPoseInSameLevel;
      else level[p.level].firstPoseInLevel = p.nextPoseInSameLevel;      
      if (p.prevPoseId>=0) pose[p.prevPoseId].nextPoseId = p.nextPoseId;
      if (p.nextPoseId>=0) pose[p.nextPoseId].prevPoseId = p.prevPoseId;      
    }    
  }  
  TmTreemap::deleteFeature (id);    
}


void TmSlamDriver2DL::hasBeenSparsifiedOut (TmFeatureId id)
{
  TmTreemap::hasBeenSparsifiedOut (id);
  if (feature[id].userFlag()==POSEX) statistic.pSparsified++;    
}


void TmSlamDriver2DL::observe (const ObservationList& observation)
{
#if ASSERT_LEVEL>=3
  assertIt();  
#endif
//  if (observation.empty()) return;
  
  // allocate current robot pose as 3 features and set initial estimate
  VmVector3 poseEst;  
  robotEstimate (poseEst);  

  int oldPoseFeature = poseFeature;
  poseFeature = newPose (relativeDistance);

  setInitialEstimateForPose (poseFeature, poseEst);


  for (int i=0; i<(int) observation.size(); i++) {
    int id = observation[i].id;
    if (id>=0) {      
      assert (0<=id && id<nrOfLandmarks);    
      int j;
      for (j=0; j<i && observation[j].id!=id; j++);      
      if (j==i && (!landmark.idx(id) || landmark[id].featureId==-1)) {
        initLandmark(id);  
        statistic.n++;
      }
      statistic.m++;
    }    
  }  
  statistic.p++;  



  NonlinearLeaf* leaf = new NonlinearLeaf (this);
  leaf->resetFlag (TmNode::CAN_BE_INTEGRATED);  
  leaf->poseFeature = poseFeature;  
  leaf->odometry = Odometry (oldPoseFeature, poseFeature, relativePose, relativePoseCov);
  leaf->observation = observation;  
  leaf->linearize (true); // use measurements when first linearizing a leaf because estimate can be way off  
  addNonlinearLeaf (leaf);
  nonlinearLeaf.push_back (leaf->index);


  clearRelativePose();  


  while ((int) nonlinearLeaf.size()>nrOfNonlinearLeaves) {
    TmNode* n = getNode (nonlinearLeaf.front());
    n->setFlag (TmNode::CAN_BE_INTEGRATED);
    n->resetFlagUpToRoot (TmNode::IS_FEATURE_PASSED_VALID);    
    if (n->parent!=NULL) n->setToBeOptimizedUpToRoot ();
    nonlinearLeaf.pop_front();
    estimateIncludes--;    
  }
  updateFeaturePassed();
}


void TmSlamDriver2DL::computeLinearEstimate ()
{
  TmTreemap::computeLinearEstimate ();  
  estimateIncludes = nonlinearLeaf.size()-1; // All leaves are incorporated now
}


void TmSlamDriver2DL::computeNonlinearEstimate ()
{
  for (int i=0; i<=estimateIncludes; i++) {
    NonlinearLeaf* n = dynamic_cast<NonlinearLeaf*> (getNode (nonlinearLeaf[i]));
    n->linearize (false); // relinearize using the estimate    
  }
  computeLinearEstimate ();
}


void TmSlamDriver2DL::robotEstimate (double& x, double& y, double &theta)
{
  double xR = feature[poseFeature].est;
  double yR = feature[poseFeature+1].est;
  double tR = feature[poseFeature+2].est;
  double c = cos(tR), s = sin(tR);  
  x = xR + c*relativePose[0] - s*relativePose[1];
  y = yR + s*relativePose[0] + c*relativePose[1];
  theta = tR + relativePose[2];
}


void TmSlamDriver2DL::robotEstimate (VmVector3& pose)
{
  robotEstimate (pose[0], pose[1], pose[2]);  
}

 
int TmSlamDriver2DL::nrOfLandmarksDefined () const
{
  int ctr=0;  
  for (int i=0; i<(int) feature.size(); i++)
    if (feature[i].isDefined() && feature[i].userFlag()==LANDMARKX) ctr++;
  return ctr;  
}

  
void TmSlamDriver2DL::nameOfFeature (char* txt, int featureId, int& n) const
{
  int status = feature[featureId].userFlag();
  if (status==LANDMARKX) {
    n = 2;
    featureToString (txt, featureId/2-landmarkBaseFeature);
  }
  else if (status==POSEX) {
    n = 3;
    sprintf(txt, "p%d", pose[(featureId-poseBaseFeature)/3].poseNr);
  }
  else {
    n = 1;
    txt[0]='?';
    txt[1]=0;
  }  
}


///***** TmSlamDriver2DL::NonlinearLeaf

TmSlamDriver2DL::NonlinearLeaf::NonlinearLeaf (TmSlamDriver2DL* tree)
  :TmNode(), poseFeature(-1), absolutePose(), odometry(), observation()
{
  this->tree = tree;  
}


TmNode* TmSlamDriver2DL::NonlinearLeaf::duplicate() const
{
  return new NonlinearLeaf (*this);
}


int TmSlamDriver2DL::NonlinearLeaf::memory() const
{
  return TmNode::memory () - sizeof(TmNode) + sizeof (TmSlamDriver2DL) + 
    observation.memory();  
}


void TmSlamDriver2DL::NonlinearLeaf::linearize (bool useMeasurement)
{  
  TmSlamDriver2DL* tree = dynamic_cast<TmSlamDriver2DL*> (this->tree);
#if ASSERT_LEVEL>=2
  tree->assertIt();
#endif
  assert (!useMeasurement || index==-1);  
  if (index!=-1) beforeChange ();
  
  
  // Create the list of features involved. We exploit that
  // TmGaussian's can have duplicate feature entries
  TmExtendedFeatureList fl;
  fl.push_back (TmExtendedFeatureId (poseFeature  , 1)); // x
  fl.push_back (TmExtendedFeatureId (poseFeature+1, 1));  // y 
  fl.push_back (TmExtendedFeatureId (poseFeature+2, 1));  // theta
  int nrOfRows = 0; // no row for new pose
  for (int i=0; i<(int) observation.size(); i++) {
    int id = observation[i].id;    
    if (id>=0) {      
      assert (0<=id && id<tree->nrOfLandmarks);    
      int f = tree->landmarkFeature(id);
      fl.push_back (TmExtendedFeatureId (f  , 1));   // x
      fl.push_back (TmExtendedFeatureId (f+1, 1)); // y
      nrOfRows += 2;      
    }    
  }
  if (!odometry.isEmpty()) {    
    assert (odometry.newPoseFeature==poseFeature);    
    fl.push_back (TmExtendedFeatureId (odometry.oldPoseFeature  , 1));  // x
    fl.push_back (TmExtendedFeatureId (odometry.oldPoseFeature+1, 1)); // y
    fl.push_back (TmExtendedFeatureId (odometry.oldPoseFeature+2, 1)); // theta 
    nrOfRows += 3;    
  }  
  for (int i=0; i<(int) fl.size(); i++) assert (!tree->feature[fl[i].id].isEmpty());

  // Construct gaussian
  gaussian.create (fl, nrOfRows); 
  if (!absolutePose.isEmpty()) addAbsolute (0, absolutePose.pose, absolutePose.poseCov);
  int j = 3;  // Column of the current landmark
  for (int i=0; i<(int) observation.size(); i++) {
    const Observation& obs = observation[i];
    if (obs.id>=0) {      
      addLandmarkMeasurement (0, j, obs.pos, obs.posCov, useMeasurement);    
      j += 2;
    }    
  }  
  if (!odometry.isEmpty()) addOdometry (j, 0, odometry.pose, odometry.poseCov);

  if (index!=-1) afterChange ();  
}


void TmSlamDriver2DL::NonlinearLeaf::addOdometry (int cOld, int cNew, const VmVector3& odo, const VmMatrix3x3& odoCov)
{
  TmSlamDriver2DL* tree = dynamic_cast<TmSlamDriver2DL*> (this->tree);  
  int fOld = gaussian.feature[cOld].id;  
  int fNew = gaussian.feature[cNew].id;  
  double oldTheta = tree->feature[fOld+2].est; // use latest estimate for linearization
  double c = cos(oldTheta), s = sin(oldTheta);
  double dX = tree->feature[fNew  ].est - tree->feature[fOld  ].est;
  double dY = tree->feature[fNew+1].est - tree->feature[fOld+1].est;  
  VmVector3 relPos = 
    {  c*dX + s*dY,
      -s*dX + c*dY,
      tree->feature[fNew+2].est - tree->feature[fOld+2].est};  

  // Odometry
  VmMatrix3x3 J3 = // Jacobian of odometry measurement with respect to old pose
    {{-c, -s,  relPos[1]},
     { s, -c, -relPos[0]},
     { 0,  0, -1 }};  
  VmMatrix3x3 J4 = // Jacobian of odometry measurement with respect to new pose
    {{ c, s, 0},
     {-s, c, 0},
     {0 , 0,1}};
  VmVector3 b =   // minus the measurement
    {-odo[0]-relPos[1]*oldTheta,
     -odo[1]+relPos[0]*oldTheta,
     -odo[2]};

  // Nonlinear odometry: 
  //    odo[0] =  c*(newPose[0]-oldPose[0]) + s*(newPose[1]-oldPose[1])
  //    odo[1] = -s*(newPose[0]-oldPose[0]) + c*(newPose[1]-oldPose[1])
  //    odo[2] = newPose[2]-oldPose[2]
  // Odometry means observing J3*oldPose + J4*newPose +b = 0
  
  VmMatrix3x3 L, LJ3, LJ4;
  VmVector3 lB;
  vmCholeskyInverse (odoCov, L); // C = J^TJ = L^-TL^-1 = {LL^T)^{-1
  vmMultiply (LJ3, L, J3); 
  vmMultiply (LJ4, L, J4);
  vmMultiply (lB, L, b);
  int i = gaussian.R.rows();  
  gaussian.R.appendRow (3, true);
  gaussian.R.store (LJ3, i, cOld);
  gaussian.R.store (LJ4, i, cNew);
  gaussian.R.storeCol (lB, i, gaussian.feature.size());
}


void TmSlamDriver2DL::NonlinearLeaf::addAbsolute (int cP, const VmVector3& pos, const VmMatrix3x3& posCov)
{
  // absolute information means observing pose = pos
  VmMatrix3x3 L;
  VmVector3 lB;  
  vmCholeskyInverse (posCov, L);
  vmMultiply (lB, L, pos);
  vmScale (lB, lB, -1);

  int i = gaussian.R.rows();  
  gaussian.R.appendRow (3, true);
  gaussian.R.store (L, i, cP);
  gaussian.R.storeCol (lB, i, gaussian.feature.size());
}


void TmSlamDriver2DL::NonlinearLeaf::addLandmarkMeasurement (int cP, int cL, const VmVector2& pos, const VmMatrix2x2& posCov, bool useMeasurement)
{
  TmSlamDriver2DL* tree = dynamic_cast<TmSlamDriver2DL*> (this->tree);  
  int fP = gaussian.feature[cP].id;
  int fL = gaussian.feature[cL].id;  
  double theta = tree->feature[fP+2].est;
  double c = cos (theta), s = sin (theta);
  VmVector2 landmark2Robot; // The landmark position in robot coordinates that is used for linearization
  if (useMeasurement) vmCopy (pos, landmark2Robot);  
  else {
    // Compute from estimate
    landmark2Robot[0] = tree->feature[fL+0].est - tree->feature[fP+0].est;
    landmark2Robot[1] = tree->feature[fL+1].est - tree->feature[fP+1].est;
    vmRotate (landmark2Robot, landmark2Robot, -theta);    
  }

  VmMatrix2x3 J5 = // Jacobian of landmark measurement with respect to new pose
    {{-c, -s, landmark2Robot[1]},
     {s,  -c, -landmark2Robot[0]}};
  VmMatrix2x2 J6 = // Jacobian of landmark measurement with respect to landmark position
    {{c, s},
     {-s, c}};
  VmVector2 b = // minus the measurement
    {-pos[0]-landmark2Robot[1]*theta, 
     -pos[1]+landmark2Robot[0]*theta};
  
  // Landmark observation means, observing J5*robotPose + J6*landmarkPosition +b = 0
  VmMatrix2x2 L;
  VmMatrix2x3 LJ5;
  VmMatrix2x2 LJ6;
  VmVector2 lB;    
  vmCholeskyInverse (posCov, L);
  vmMultiply (LJ5, L, J5);
  vmMultiply (LJ6, L, J6);
  vmMultiply (lB, L, b);    
  int i = gaussian.R.rows();
  gaussian.R.appendRow (2, true);
  gaussian.R.store (LJ5, i, cP);
  gaussian.R.store (LJ6, i, cL);
  gaussian.R.storeCol (lB, i, gaussian.feature.size()); 
}


bool TmSlamDriver2DL::canBeSparsifiedOut (TmFeatureId id)
{
  updateFeaturePassed ();  
  if (id<0 || id>=(int) feature.size()) return false;
  if (feature[id].userFlag()!=POSEX) return false;
  if (!feature[id].isFlag (TmFeature::CAN_BE_MARGINALIZED_OUT)) return false;
  int poseId = (id-poseBaseFeature)/3;      
  if (pose[poseId].distance > pose[poseIndex].distance-sparsificationDistance) return false;  
  XycVector<TmNode*> node;
  TmNode* mag = feature[id].marginalizationNode;
  assert (mag!=NULL);  
  if (!mag->isFlag (TmNode::IS_OPTIMIZED)) return false;  
  findLeavesInvolving (id, node);
  assert (node.size()>=2);  
  for (int i=0; i<(int) node.size(); i++) for (int j=i+1; j<(int) node.size(); j++) 
    if (sharedLandmarks (node[i]->gaussian.feature, node[j]->gaussian.feature)<2) return false;      
  for (int i=0; i<(int) node.size(); i++) if (!node[i]->isFlag(TmNode::CAN_BE_INTEGRATED)) return false;      
  return true;  
}


void TmSlamDriver2DL::checkForSparsification (TmNode* n)
{
  if (!n->isFlag (TmNode::IS_OPTIMIZED)) return;
  TmExtendedFeatureList fl;
  n->computeFeaturesInvolved (fl);
  for (int i=0; i<(int) fl.size(); i++) {
    TmFeature& f = feature[fl[i].id];    
    if (f.marginalizationNode==n && f.userFlag()==POSEX && canBeSparsifiedOut (fl[i].id))
      sparsifyOut (fl[i].id, 3); // Always sparsify out POSEX, POSEY, POSETHETA
  }
}



int TmSlamDriver2DL::sharedLandmarks (const TmExtendedFeatureList& a, const TmExtendedFeatureList& b)
{
  int ctr = 0;  
  for (int i=0; i<(int) a.size(); i++) {
    int id = a[i].id;    
    int j;
    for (j=0; j<i; j++) if (a[j].id==id) break;
    if (j>=i && feature[id].userFlag()==LANDMARKX) {      
      for (j=0; j<(int) b.size(); j++) 
        if (b[j].id==id) break;
      if (j<(int) b.size()) ctr++;
    }    
  }
  return ctr;  
}


void TmSlamDriver2DL::assertIt () const
{
#ifndef NDEBUG
  if (poseIndex==-1) return;  
  TmTreemap::assertIt();
  assert (0<=poseIndex && poseIndex<(int) feature.size());
  assert (poseBaseFeature + 3*poseIndex==poseFeature);  
  const TmFeature &p = feature[poseFeature];  
  assert (!p.isEmpty());
  assert (p.isFlag (TmFeature::CAN_BE_MARGINALIZED_OUT));
  assert (!p.isFlag (TmFeature::CAN_BE_SPARSIFIED));
  assert (feature[poseFeature+0].userFlag()==POSEX);
  assert (feature[poseFeature+1].userFlag()==POSEY);
  assert (feature[poseFeature+2].userFlag()==POSETHETA);  
#endif
}


int TmSlamDriver2DL::memory () const
{
  int mem = TmTreemap::memory()+ sizeof(TmSlamDriver2DL) - sizeof (TmTreemap);
  mem += pose.capacity() * sizeof(Pose);
  // We ignore the nonlinearleaves (they are not many)
  return mem;  
}


void TmSlamDriver2DL::onlyUpdateEstimatesForLandmarks (int from, int to, bool setDontUpdateFlag)
{
  onlyUpdateEstimatesFor (landmarkFeature(from), landmarkFeature(to), setDontUpdateFlag);
}


void TmSlamDriver2DL::onlyUpdateLevel (int level, bool setDontUpdateFlag)
{
  if (root==NULL) return;  
  updateFeaturePassed();  
  if (setDontUpdateFlag) root->setFlagEverywhere (TmNode::DONT_UPDATE_ESTIMATE);  
  if (!this->level.idx(level)) return;  
  int i=this->level[level].firstPoseInLevel;
  while (i>=0) {
    Pose& p = pose[i];    
    assert (p.level==level);    
    if (p.featureId>=0) {
      TmFeature& f = feature[p.featureId];
      if (f.marginalizationNode!=NULL) f.marginalizationNode->resetFlagUpToRoot (TmNode::DONT_UPDATE_ESTIMATE);
    }
    i = p.nextPoseInSameLevel;    
  }
  i=this->level[level].firstLandmarkInLevel;  
  while (i>=0) {    
    Landmark& lm = landmark[i];    
    assert (lm.level==level);    
    if (lm.featureId>=0) {
      TmFeature& f = feature[lm.featureId];
      if (f.marginalizationNode!=NULL) f.marginalizationNode->resetFlagUpToRoot (TmNode::DONT_UPDATE_ESTIMATE);
    }
    i = lm.nextLandmarkInSameLevel;    
  }  
}


double TmSlamDriver2DL::chi2 (const ObservationList& observation)
{
  double chi2E = 0;
  VmVector3 pose;
  robotEstimate (pose);  
  for (int i=0; i<(int) observation.size(); i++) {
    const Observation& obs = observation[i];    
    VmVector2 lm;
    if (!landmarkEstimate (obs.id, lm)) return vmInf();    
    if (!finite(lm[0]) || !finite(lm[1])) return vmInf();    
    VmVector2 posErr;
    double c = cos(pose[2]), s = sin(pose[2]);    
    posErr[0] =  c*(lm[0]-pose[0]) + s*(lm[1]-pose[1]) - obs.pos[0];
    posErr[1] = -s*(lm[0]-pose[0]) + c*(lm[1]-pose[1]) - obs.pos[1];
    VmMatrix2x2 cInv;
    vmInverse (cInv, obs.posCov);    
    vmJCKtAdd (chi2E, posErr, cInv, posErr);    
  }  
  return chi2E;  
}

