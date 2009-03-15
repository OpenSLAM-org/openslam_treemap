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

/*!\file tmNode.cc 
   \brief Implementation of \c TmNode
   \author Udo Frese

   Contains the implementation of class \c TmNode.
*/

#include "tmNode.h"
#include "tmTreemap.h"
#include <algorithm>


TmNode::TmNode ()
  : index (-1), tree (NULL), parent(0), updateCost(0), worstCaseUpdateCost (0), 
    featurePassed(), linearizationPointFeature(-1),
    gaussian(), firstFeaturePassed(-1), status (0) 
{
  child[0] = child[1] = NULL;
}


TmNode::~TmNode() {};
  

TmNode* TmNode::duplicate () const
{
  return new TmNode (*this);
}


void TmNode::resetFlagEverywhere (int flag)
{
  resetFlag (flag);
  if (!isLeaf()) {
    child[0]->resetFlagEverywhere (flag);
    child[1]->resetFlagEverywhere (flag);
  }
}


void TmNode::setFlagEverywhere (int flag, bool stopIfFound)
{
  if (stopIfFound && isFlag(flag)) return;  
  setFlag (flag);
  if (!isLeaf()) {
    child[0]->setFlagEverywhere (flag);
    child[1]->setFlagEverywhere (flag);
  }
}


int TmNode::mergeFeatureLists ()
{
  int n = 0;  
  const TmExtendedFeatureList& a = child[0]->featurePassed;
  const TmExtendedFeatureList& b = child[1]->featurePassed;
  featurePassed.resizeWithUndefinedData (a.size()+b.size());
  
  const TmExtendedFeatureId* fA    = a.begin();
  const TmExtendedFeatureId* fAEnd = a.end();  
  const TmExtendedFeatureId* fB    = b.begin();
  const TmExtendedFeatureId* fBEnd = b.end();
  TmExtendedFeatureId* result = featurePassed.begin();
  int lpF = linearizationPointFeature;  
  TmFeature* tF = tree->feature.begin();  

  if (fA==fAEnd) goto fAEmpty;
  if (fB==fBEnd) goto fBEmpty;  
// Merge a and b
  while (true) {
    int aId = fA->id;
    int bId = fB->id;    
    n++;      
    if (bId<aId) { // only in 'b'
      int bCount = fB->count;
      if (bId==lpF || bCount<tF[bId].totalCountOfExistingFeature()) {
        result->id    = bId;
        result->count = bCount;        
        result++;
      }
      else tF[bId].marginalizationNode = this;        
      fB++;
      if (fB==fBEnd) goto fBEmpty;
    }
    else if (bId==aId) { // in both
      int countSum = fA->count + fB->count;
      if (bId==lpF || countSum<tF[bId].totalCountOfExistingFeature()) {
        result->id    = bId;
        result->count = countSum;
        result++;
      }
      else tF[bId].marginalizationNode = this;      
      fB++;
      fA++;
      if (fB==fBEnd) goto fBEmpty;
      if (fA==fAEnd) goto fAEmpty;
    }
    else { // only in 'a'
      int aCount = fA->count;      
      if (aId==lpF || aCount<tF[aId].totalCountOfExistingFeature()) {
        result->id    = aId;
        result->count = aCount;        
        result++;
      }
      else tF[aId].marginalizationNode = this;        
      fA++;
      if (fA==fAEnd) goto fAEmpty;
    }
  };  
  // We will never reach this line

 fAEmpty: 
  // Copy the remains of b
  while (fB!=fBEnd) {
    n++;    
    int bId = fB->id;    
    int bCount = fB->count;    
    if (bId==lpF || bCount<tF[bId].totalCountOfExistingFeature()) {
      result->id    = bId;
      result->count = bCount;      
      result++;
    }
    else tF[bId].marginalizationNode = this;        
    fB++;
  }
  featurePassed.eraseAfter (result);
  return n;

 fBEmpty:
  // Copy the remains of a
  while (fA!=fAEnd) {
    n++;    
    int aId = fA->id;    
    int aCount = fA->count;    
    if (aId==lpF || aCount<tF[aId].totalCountOfExistingFeature()) {
      result->id = aId;
      result->count = aCount;      
      result++;
    }
    else tF[aId].marginalizationNode = this;        
    fA++;
  }
  featurePassed.eraseAfter (result);
  return n;
  
/* Old Code
    const TmExtendedFeatureList& a = child[0]->featurePassed;
    const TmExtendedFeatureList& b = child[1]->featurePassed;

    // merge the feature list and remove eliminated ones
    int szA = a.size();
    int szB = b.size();
    n=0;    
    featurePassed.reserve (szA+szB);
    featurePassed.clear();    
    int i=0, j=0;
    while (i<szA || j<szB) {
        if (j<szB && (i>=szA || b[j].id<a[i].id)) {
            featurePassed.push_back (b[j]); // Only in 'b'
            j++;
            n++;            
        }        
        else if (j<szB && i<szA && a[i].id==b[j].id) {  // is contained in both
          featurePassed.push_back (a[i]);
          featurePassed.back().count += b[j].count;          
          i++;          
          j++;  
          n++;          
        }
        else {
          featurePassed.push_back (a[i]); // Only in 'a'
          i++;          
          n++;
        }
        TmExtendedFeatureId& f = featurePassed.back(); 
        assert (f.count<=tree->feature[f.id].totalCount());      
        if (f.id != linearizationPointFeature && f.count == tree->feature[f.id].totalCount()) {
          tree->feature[f.id].marginalizationNode = this;          
          featurePassed.pop_back();
        }        
    }
*/
}


void TmNode::updateFeaturePassed ()
{
  if (isFlag(IS_FEATURE_PASSED_VALID)) return;
  int n;  
  featurePassed.clear(); 
  if (isLeaf()) {
    // is leaf so update from \c gaussian
    linearizationPointFeature = gaussian.linearizationPointFeature;
    n = gaussian.feature.size();    
    updateCost = worstCaseUpdateCost = updateGaussianCost (n); 
    tree->stat.accumulatedOptimizationCost += updateFeaturePassedCost (n);  
    featurePassed = gaussian.feature;
    sort (featurePassed.begin(), featurePassed.end());
    // Add up duplicates and remove features marginalized out
    int j=-1;  // j is the last entry of the already merged part of list
    for (int i=0; i<(int) featurePassed.size(); i++) {
      if (j>=0 && featurePassed[i].id==featurePassed[j].id) 
        featurePassed[j].count += featurePassed[i].count; // merge into featurePassed[j]
      else {
        j++; // add new entry to merged part
        featurePassed[j] = featurePassed[i];
      }
      TmExtendedFeatureId& f = featurePassed[j];      
      assert (f.count<=tree->feature[f.id].totalCount());      
      if (f.id!=linearizationPointFeature && f.count==tree->feature[f.id].totalCount()) {
        tree->feature[f.id].marginalizationNode = this;        
        j--;      
      }      
    }
    featurePassed.resize (j+1);
  }
  else {   // update recursively
    child[0]->updateFeaturePassed();
    child[1]->updateFeaturePassed();

    // set linearizationPointFeature
    if (child[0]->linearizationPointFeature>=0 && child[1]->linearizationPointFeature>=0) // both are rotation invariant
      linearizationPointFeature = child[0]->linearizationPointFeature; // so take the first one
    else linearizationPointFeature = -1; // the distribution is not rotation invariant    

    // set worstCaseUpdateCost
    double wc0 = child[0]->worstCaseUpdateCost;
    double wc1 = child[1]->worstCaseUpdateCost;
    if (wc0<wc1) worstCaseUpdateCost = wc1;
    else worstCaseUpdateCost = wc0;

    n = mergeFeatureLists ();    

    updateCost = updateGaussianCost (n);    
    worstCaseUpdateCost += updateCost;    
    tree->stat.accumulatedOptimizationCost += updateFeaturePassedCost (n);  

    // set feature flags
    if (child[0]->isFlag (CAN_BE_INTEGRATED) && child[1]->isFlag (CAN_BE_INTEGRATED))
      setFlag (CAN_BE_INTEGRATED);
    else resetFlag (CAN_BE_INTEGRATED);
    if (child[0]->isFlag (DONT_UPDATE_ESTIMATE) && child[1]->isFlag (DONT_UPDATE_ESTIMATE))
      setFlag (DONT_UPDATE_ESTIMATE);    
    else resetFlag (DONT_UPDATE_ESTIMATE);    
  }
  setFlag (IS_FEATURE_PASSED_VALID);
}


void TmNode::updateGaussian ()
{
  if (isFlag(IS_GAUSSIAN_VALID)) return;  
  updateFeaturePassed ();
  if (isLeaf()) {
    TmExtendedFeatureList fl;
    fl.reserve (gaussian.feature.size());
    addMarginalizedFeatures (fl, gaussian.feature);
    firstFeaturePassed = fl.size();    
    for (int i=0; i<(int) featurePassed.size(); i++)
      fl.push_back (TmExtendedFeatureId (featurePassed[i].id, 0));
    TmGaussian myGaussian (fl, gaussian.rows());
    myGaussian.multiply (gaussian, 0);
    myGaussian.setLinearizationPoint (linearizationPointFeature, 0); // TODO 0 is wrong
    myGaussian.triangularize (tree->workspace);
    myGaussian.compress ();
    gaussian.transferFrom (myGaussian);
  }  
  else {
        // update recursively
    child[0]->updateGaussian ();
    child[1]->updateGaussian ();

    TmExtendedFeatureList fl;
    fl.reserve (child[0]->featurePassed.size()+child[1]->featurePassed.size());
    addMarginalizedFeatures (fl, child[0]->featurePassed);
    addMarginalizedFeatures (fl, child[1]->featurePassed);
    firstFeaturePassed = fl.size();    
    for (int i=0; i<(int) featurePassed.size(); i++)
      fl.push_back (featurePassed[i]);
    TmGaussian myGaussian( fl, child[0]->gaussian.rows() + child[1]->gaussian.rows());
    myGaussian.multiply (child[0]->gaussian, child[0]->firstFeaturePassed);
    // TODO we should rotate here
    myGaussian.multiply (child[1]->gaussian, child[1]->firstFeaturePassed);    
    // TODO we should rotate here
    myGaussian.setLinearizationPoint (linearizationPointFeature, 0); // TODO 0 is wrong
    myGaussian.triangularize (tree->workspace);    
    myGaussian.compress ();
    gaussian.transferFrom (myGaussian);
  }
#if ASSERT_LEVEL>=1
  gaussian.assertIt ();  
#endif
  tree->stat.accumulatedUpdateCost += updateCost;
  tree->stat.nrOfGaussianUpdates++;  
  setFlag (IS_GAUSSIAN_VALID);  
}


void TmNode::computeFeaturesInvolved (TmExtendedFeatureList& list) const
{
  if (isLeaf()) {
    // take from gaussian
    list.clear();
    list.reserve (gaussian.feature.size());
    for (int i=0; i<(int) gaussian.feature.size(); i++)
      list.push_back (gaussian.feature[i]);
    sumUp (list);    
  }
  else merge (list, child[0]->featurePassed, child[1]->featurePassed);  
}


int TmNode::getHeight (int max) const
{
  if (isLeaf() || max==1) return 1;
  else {
    int h1 = child[0]->getHeight();
    int h2 = child[1]->getHeight();
    if (h1<h2) return h2;
    else return h1;
  }
}


bool TmNode::isAncestor(TmNode* n) const
{
  const TmNode* n2 = this;  
  while (n2!=NULL) {
    if (n2==n) return true;
    n2 = n2->parent;
  }
  return false;
}


void TmNode::addMarginalizedFeatures (TmExtendedFeatureList& fl, TmExtendedFeatureList& fl2)
{
  for (int i=0; i<(int) fl2.size(); i++) {
    int feat = fl2[i].id;
    int j;    
    for (j=0; j<(int) fl.size(); j++) if (fl[j].id==feat) break;
    if (j==(int) fl.size() && !isElement (featurePassed, feat)) fl.push_back (TmExtendedFeatureId (fl2[i].id, 0));
  }  
}


void TmNode::estimate ()
{
  assert (isFlag(IS_GAUSSIAN_VALID));
  if (isFlag (DONT_UPDATE_ESTIMATE)) return;  
  XymVector& v = tree->workspace;  
  v.resize (gaussian.feature.size());
  // Fill lower part of v with estimates for features already passed
  for (int i=firstFeaturePassed; i<v.size(); i++)  
    v[i] = tree->feature[gaussian.feature[i].id].est;
  // We should rotate here
  // Compute estimate for upper part as conditioned mean
  gaussian.mean (v, firstFeaturePassed);
  for (int i=0; i<firstFeaturePassed; i++) {    
    tree->feature[gaussian.feature[i].id].est = v[i];
    assert (finite(v[i]));
  }  
  if (!isLeaf()) {
    child[0]->estimate ();
    child[1]->estimate ();
  }  
}


void TmNode::estimateUsingRCompressed ()
{
  if (isFlag (DONT_UPDATE_ESTIMATE)) return;  
  tree->workspaceFloat.resize (gaussian.feature.size()+5);  // We need 5 as additional space for the SSE implementation
  // Fill v with estimates for features already passed in reverse order
  float* v  = tree->workspaceFloat.begin();
  float* vE = v + gaussian.feature.size() + 1;  
  *v = 1; // homogenous 1
  v++;    
  TmExtendedFeatureId* srcP  = gaussian.feature.end() -1;
  TmExtendedFeatureId* srcPE = gaussian.feature.begin() + firstFeaturePassed-1;
  while (srcP!=srcPE) {
    *v = tree->feature[srcP->id].est;
    srcP--;
    v++;
  }
  // Compute mean of remaining features conditioned on the one stored in v
  gaussian.meanCompressed (tree->workspaceFloat.begin(), firstFeaturePassed);
  // Store the result in the feature estimates
  while (v!=vE) {
    tree->feature[srcP->id].est = *v;
    srcP--;
    v++;
  }  
  // Go recursively down
  if (!isLeaf()) {
    child[0]->estimateUsingRCompressed ();
    child[1]->estimateUsingRCompressed ();
  }  
}



void TmNode::checkLinearizationPointFeature ()
{
}


void TmNode::changeGaussian (const TmGaussian& gaussian)
{
  assert (isLeaf());  
  beforeChange ();  
  this->gaussian = gaussian;
  afterChange ();  
}


void TmNode::beforeChange ()
{
  assert (isLeaf());  
  // invalidate all old marginalization nodes, subtract from totalCount
  for (int i=0; i<(int) this->gaussian.feature.size(); i++) {
    TmFeatureId id = this->gaussian.feature[i].id;
    TmFeature& feat = tree->feature[id];
    feat.addTotalCount (-this->gaussian.feature[i].count);
    if (feat.marginalizationNode!=NULL) {
      feat.marginalizationNode->resetFlagUpToRoot (IS_FEATURE_PASSED_VALID | IS_GAUSSIAN_VALID);
      feat.marginalizationNode=NULL;
    }    
  }
}


void TmNode::afterChange ()
{
  assert (isLeaf());  
  // invalidate all new marginalization nodes, extend feature vector if necessary
  // and add to totalCount
  for (int i=0; i<(int) this->gaussian.feature.size(); i++) {
    TmFeatureId id = this->gaussian.feature[i].id;
    if (id>=(int) tree->feature.size()) tree->feature.resize (id+1);    
    TmFeature& feat = tree->feature[id];
    feat.addTotalCount (this->gaussian.feature[i].count);    
    if (feat.marginalizationNode!=NULL) {
      feat.marginalizationNode->resetFlagUpToRoot (IS_FEATURE_PASSED_VALID | IS_GAUSSIAN_VALID);
      feat.marginalizationNode = NULL;
    }    
  }  
  resetFlagUpToRoot (IS_FEATURE_PASSED_VALID | IS_GAUSSIAN_VALID);
  // Ensure, that at least for the linearization point features
  // an initial estimate is available.
  if (gaussian.linearizationPointFeature>=0)
    tree->setInitialEstimate (gaussian.linearizationPointFeature, gaussian.linearizationPoint);

  // If the featuresPassed list changed, set all ancestors as to be optimized.
  TmExtendedFeatureList oldFl = featurePassed;
  updateFeaturePassed ();
  if (oldFl!=featurePassed) setToBeOptimizedUpToRoot();
}


int TmNode::rowsBelow () const
{
  if (isLeaf()) return gaussian.rows();
  else return child[0]->rowsBelow () + child[1]->rowsBelow ();
}



int TmNode::degreeOfFeature (int feature) const
{
  assert (isFlag(IS_FEATURE_PASSED_VALID));  
  int ctr = 0;  
  if (child[0]!=NULL && isElement (child[0]->featurePassed, feature)) ctr++;
  if (child[1]!=NULL && isElement (child[1]->featurePassed, feature)) ctr++;
  return ctr;
}


void TmNode::assertIt () const
{
#ifndef NDEBUG
  assert (0<=index && index<(int) tree->node.size());  
  assert (tree->node[index]==this);
  assertInTree ();  
  if (isLeaf()) {
  }
  else {
    assert (child[0]->parent==this);
    assert (child[1]->parent==this);    
  }
  if (isFlag(IS_FEATURE_PASSED_VALID)) {
    if (!isLeaf()) {
      assert (worstCaseUpdateCost+costEps()>=child[0]->worstCaseUpdateCost+updateCost);
      assert (worstCaseUpdateCost+costEps()>=child[1]->worstCaseUpdateCost+updateCost);
//      assert (!isIntersectionEmpty (child[0]->featurePassed, child[1]->featurePassed)); // TODO remove      
    }    
    else assert (isFlag(IS_OPTIMIZED));    
    for (int i=0; i<(int) featurePassed.size(); i++) {
      const TmExtendedFeatureId& f = featurePassed[i];    
      assert (f.id>=0 && f.id<(int) tree->feature.size());
      assert (tree->feature[f.id].totalCount()>=f.count);
#if ASSERT_LEVEL >=4
      bool below, notBelow;
      isRepresentedWithRespectToNode (f.id, below, notBelow);
      assert (below && notBelow);      
#endif
    }
  }
  if (tree->isGaussianValidValid && isFlag(IS_GAUSSIAN_VALID)) {
    assert (gaussian.isTriangular);
    assert (gaussian.R.isValid() || !gaussian.RCompressed.empty());    
    if (gaussian.R.isValid()) {      
      assert (gaussian.R.cols()==(int) gaussian.feature.size()+1);
      assert (gaussian.R.rows()<=(int) gaussian.feature.size()+1);    
    }
    for (int i=0; i<(int) gaussian.feature.size(); i++) {
      int id = gaussian.feature[i].id;
      assert (id>=0 && id<(int) tree->feature.size());
      if (isFlag(IS_FEATURE_PASSED_VALID)) {
        if (isElement (featurePassed, id)) {
          assert (i>=firstFeaturePassed);
        }
        else {
          assert (i<firstFeaturePassed);
        }
      }      
    }
    gaussian.assertIt();    
  }  
#endif
}

  
void TmNode::recursiveAssertIt () const
{
  assertIt ();
  if (!isLeaf()) {
    child[0]->recursiveAssertIt();
    child[1]->recursiveAssertIt();
  }
}


void TmNode::isRepresentedWithRespectToNode (TmFeatureId id, bool& below, bool& notBelow) const
{
  below = notBelow = false;  
  tree->root->recursiveIsRepresentedWithRespectToNode (this , false, id, below, notBelow);  
}


void TmNode::recursiveIsRepresentedWithRespectToNode (const TmNode* n, bool isN2BelowN, TmFeatureId id, bool& below, bool& notBelow) const
{
  if (this==n) isN2BelowN = true;
  if (isLeaf()) {
    for (int i=0; i<(int) gaussian.feature.size(); i++) 
      if (gaussian.feature[i].id==id) {
        if (isN2BelowN) below = true;
        else notBelow = true;
      }    
  }  
  else {
    child[0]->recursiveIsRepresentedWithRespectToNode (n, isN2BelowN, id, below, notBelow);
    child[1]->recursiveIsRepresentedWithRespectToNode (n, isN2BelowN, id, below, notBelow);
  }  
}


int TmNode::nrOfNodes() const
{
  if (isLeaf()) return 1;
  else return 1+child[0]->nrOfNodes()+child[1]->nrOfNodes();
}


int TmNode::sprintTreeIndex (char* s)
{
  int ctr=0;  
  if (child[0]!=NULL) {
    s[ctr++] = '(';    
    ctr += child[0]->sprintTreeIndex (s+ctr);   
  }  
  sprintf (s+ctr, "%d ", index);  
  ctr+= strlen(s+ctr);  
  if (child[1]!=NULL) {
    ctr += child[1]->sprintTreeIndex (s+ctr);   
    s[ctr++] = ')';
  }
  s[ctr]='\0';  
  return ctr;  
}

void TmNode::recursivePrintTreeIndex ()
{
  if (child[0]!=NULL) {
    printf("(");    
    child[0]->recursivePrintTreeIndex ();   
  }  
  printf ("%d ", index);  
  if (child[1]!=NULL) {
    child[1]->recursivePrintTreeIndex ();   
    printf(")");    
  }
}

void TmNode::printTreeIndex ()
{
  recursivePrintTreeIndex ();
  printf("\n");
}


void TmNode::nodesBelow (XycVector<TmNode*>& nodes)
{
  nodes.push_back (this);
  if (!isLeaf()) {
    child[0]->nodesBelow (nodes);
    child[1]->nodesBelow (nodes);
  }
}


int TmNode::whichSideOf (const TmNode* n) const
{
  const TmNode* n2 = this;
  if (n==this) return 3;  
  while (n2!=NULL) {
    if (n2==n->child[0]) return 0;
    else if (n2==n->child[1]) return 1;
    n2 = n2->parent;
  }
  return false;  
}


void TmNode::makeChildNr (int childNr)
{
  if (parent==NULL) return;
  if (parent->child[childNr]==this) return;
  parent->child[1-childNr] = parent->child[childNr];
  parent->child[childNr]   = this;  
}


void TmNode::moveTo (TmNode* above, bool invalidateFeaturePassedOnly)
{  
  assert (parent!=NULL && this!=above);  
  TmNode* n = parent;  // That is the node that is actually moving

  // Invalidate everything from s's parent to above's parent
  if (invalidateFeaturePassedOnly) {
    n->resetFlagUpToRoot (TmNode::IS_FEATURE_PASSED_VALID);
    if (above->parent!=NULL) above->parent->resetFlagUpToRoot (TmNode::IS_FEATURE_PASSED_VALID);
  }
  else { 
    // invalidate from \c this and \c above up to the root
    n->resetFlagUpToRoot (TmNode::IS_FEATURE_PASSED_VALID | TmNode::IS_GAUSSIAN_VALID);
    if (above->parent!=NULL) above->parent->resetFlagUpToRoot (TmNode::IS_FEATURE_PASSED_VALID | TmNode::IS_GAUSSIAN_VALID);
    if (true) {       // false for oneKLMove optimization
      TmNode* lca = leastCommonAncestor (this, above);
      // mark nodes as to be optimized from \c this and \c above up to lca
      for (TmNode* n=parent; n!=lca; n=n->parent) n->setToBeOptimized();    
      if (above!=lca) for (TmNode* n=above->parent; n!=lca; n=n->parent) n->setToBeOptimized();
      lca->setToBeOptimized();  // TODO: comment this in for old code
    }    
  }  
  
  int wcS = whichChild();  
  int wcN = n->whichChild();

  // Change links
  if (n==above) above = n->child[1-wcS];  
  // remove n with s
  if (n->parent!=NULL) n->parent->child[wcN]     = n->child[1-wcS];
  else tree->root = n->child[1-wcS];  
  n->child[1-wcS]->parent   = n->parent;  
  // insert n with s above a
  n->parent                 = above->parent;
  int wcA = above->whichChild ();
  if (wcA==-1) wcA=0; // treat root as a left child
  if (above->parent!=NULL) above->parent->child[wcA] = n;
  else tree->root = n;  
  above->parent             = n;
  n->child[wcA]             = above;
  n->child[1-wcA]           = this;

#if ASSERT_LEVEL>=3
  tree->assertIt();
#endif
}


void TmNode::setToBeOptimized ()
{
  if (isFlag(IS_OPTIMIZED)) {
    resetFlag (IS_OPTIMIZED);
    tree->optimizer.optimizationQueue.push_back (this->index);    
  }
}


void TmNode::setToBeOptimizedUpToRoot ()
{
  TmNode* n=this;
  if (n->isLeaf()) n = n->parent;  
  while (n!=NULL) {
    n->setToBeOptimized ();
    n=n->parent;
  }  
}


TmNode* TmNode::leastCommonAncestor (TmNode* a, TmNode* b)
{
  if (a==NULL || b==NULL) return NULL;  
  int ctrA=0, ctrB=0;
  for (TmNode* n=a; n!=NULL; n=n->parent) ctrA++; 
  for (TmNode* n=b; n!=NULL; n=n->parent) ctrB++; 
  while (ctrA>ctrB) {
    a=a->parent;
    ctrA--;
  }
  while (ctrB>ctrA) {
    b=b->parent;
    ctrB--;
  }
  while (a!=b) {
    a=a->parent;
    b=b->parent;
  }  
  return a;  
}

void TmNode::recursiveAssertFlag (int flag)
{
  assert (isFlag(flag));
  if (!isLeaf()) {
    child[0]->recursiveAssertFlag (flag);
    child[1]->recursiveAssertFlag (flag);
  }  
}


void TmNode::recursiveAddLeavesInvolving (TmFeatureId id, XycVector<TmNode*>& node, int& ctr) 
{
  if (isLeaf()) {
//    if (isElement (gaussian.feature, id)) node.push_back (this);
    bool first = true;    
    for (int i=0; i<(int) gaussian.feature.size(); i++)
      if (gaussian.feature[i].id==id) {
        if (first) node.push_back (this);
        first = false;
        ctr += gaussian.feature[i].count;
      }    
  }
  else {
    if (isElement (child[0]->featurePassed, id)) child[0]->recursiveAddLeavesInvolving (id, node, ctr);
    if (isElement (child[1]->featurePassed, id)) child[1]->recursiveAddLeavesInvolving (id, node, ctr);
  }
}


void  TmNode::assertInTree () const
{
  const TmNode* n = this;  
  while (n!=NULL) {
    assert (n->tree == tree);    
    if (n->parent!=NULL) {
      assert (n->parent->child[0]==n || n->parent->child[1]==n);
    }
    else assert (tree->root==n);
    n = n->parent;
  }
}


void TmNode::recursivelyIdentifyFeature (int from, int to)
{
  if (isLeaf()) {
    for (int i=0; i<(int) gaussian.feature.size(); i++) {      
      if (gaussian.feature[i].id==from) {
        gaussian.feature[i].id = to;
        resetFlagUpToRoot (TmNode::IS_FEATURE_PASSED_VALID | TmNode::IS_GAUSSIAN_VALID);
        setToBeOptimizedUpToRoot ();
      }
    }
  }
  else {
    if (isElement (child[0]->featurePassed, from)) child[0]->recursivelyIdentifyFeature (from, to);
    if (isElement (child[1]->featurePassed, from)) child[1]->recursivelyIdentifyFeature (from, to);
  }  
}


int TmNode::memory() const
{
  int mem = sizeof (TmNode);
  mem += featurePassed.capacity() * sizeof(TmExtendedFeatureId);
  mem += gaussian.memory() - sizeof(TmGaussian); // TmGaussian itself is included in sizeof(*this);  
  return mem;  
}


int TmNode::recursiveMemory() const
{
  if (isLeaf()) return memory();
  else return memory() + child[0]->recursiveMemory() + child[1]->recursiveMemory();
}

