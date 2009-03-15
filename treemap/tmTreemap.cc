/*!
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


/*!\file tmTreemap.cc 
   \brief Implementation of \c TmTreemap
   \author Udo Frese
 
   Contains the implementation of class \c TmTreemap representing a
   a treemap with all algorithms involved.
*/

#include <utility>
#include <set>
#include <stdlib.h>
#include "tmTreemap.h"

#ifdef linux
#include <sys/time.h>
#include <sys/resource.h>
#endif

float tmNan = nan("NAN");

#define JOINFACTOR 1.0 //! TODO: originally we had 1.5*

TmTreemap::TmTreemap()
  :root (NULL), node(), unusedNodes(), isEstimateValid (true), isGaussianValidValid(true), feature(), 
   optimizer(), stat(), workspace(), workspaceFloat()
{
  for (int i=0; i<MAX_FEATURE_BLOCK_SIZE; i++) firstUnusedFeature[i]=-1;  
}

TmTreemap::TmTreemap (const TmTreemap& tm)
  :root (NULL), node(), unusedNodes (), isEstimateValid (true), isGaussianValidValid(true), feature(), 
   optimizer(), stat(), workspace(), workspaceFloat()
{
  *this = tm;
}


TmTreemap::TmTreemap (int nrOfMovesPerStep, int maxNrOfUnsuccessfulMoves)
  :root (NULL), node(), unusedNodes(), isEstimateValid (true), isGaussianValidValid(true), feature(), 
   optimizer(), stat(), workspace(), workspaceFloat()
{
  for (int i=0; i<MAX_FEATURE_BLOCK_SIZE; i++) firstUnusedFeature[i]=-1;  
  create (nrOfMovesPerStep, maxNrOfUnsuccessfulMoves);
}


void TmTreemap::create (int nrOfMovesPerStep, int maxNrOfUnsuccessfulMoves)
{
  clear ();  
  optimizer.create (this, nrOfMovesPerStep, maxNrOfUnsuccessfulMoves);
}


TmTreemap& TmTreemap::operator = (const TmTreemap& tm)
{
  if (&tm==this) return *this;
  isEstimateValid = tm.isEstimateValid;
  isGaussianValidValid = tm.isGaussianValidValid;  
  feature = tm.feature;
  optimizer = tm.optimizer;  
  optimizer.tree = this;  
  node = tm.node;  // Will be overwritten with correct pointer by \c recursiveCopyTreeFrom  
  unusedNodes = tm.unusedNodes;  
  stat = tm.stat;  
  workspace = tm.workspace;  
  workspaceFloat = tm.workspaceFloat;  
  for (int i=0; i<MAX_FEATURE_BLOCK_SIZE; i++) firstUnusedFeature[i] = tm.firstUnusedFeature[i];
  // We reset all marginalization node pointers to NULL for which we
  // cannot compute the involved features. This is necessary since
  // recursiveCopyTreeFrom can only change pointers to nodes where the
  // featurePassed lists are valid. It is allowed, since these nodes
  // are recomputed anyway thereby setting \c .marginalizationNode.
  for (int i=0; i<(int) feature.size(); i++) {
    TmNode* mn = feature[i].marginalizationNode;
    if (mn!=NULL && !mn->isLeaf() && 
        (!mn->child[0]->isFlag(TmNode::IS_FEATURE_PASSED_VALID) || !mn->child[1]->isFlag(TmNode::IS_FEATURE_PASSED_VALID)))
      feature[i].marginalizationNode = NULL;
  }  
  root = recursiveCopyTreeFrom (tm.root);
#if ASSERT_LEVEL>=2
  TmTreemap::assertIt (); 
   // Don't use the virtual one here, because member variables of a 
   // derived class have not been copied yet
#endif
  return *this;  
}


TmNode* TmTreemap::recursiveCopyTreeFrom (const TmNode* n2)
{
  if (n2==NULL) return NULL;  
  TmNode* nCopy = n2->duplicate();  
  nCopy->tree = this;
  TmExtendedFeatureList fl;
  n2->computeFeaturesInvolved (fl);
  for (int i=0; i<(int) fl.size(); i++) {
    TmFeature& feat = feature[fl[i].id];
    if (feat.marginalizationNode==n2) feat.marginalizationNode = nCopy;
  }
  node[n2->index] = nCopy;
  if (!n2->isLeaf()) {
    nCopy->child[0] = recursiveCopyTreeFrom (n2->child[0]);
    nCopy->child[0]->parent = nCopy;    
    nCopy->child[1] = recursiveCopyTreeFrom (n2->child[1]);    
    nCopy->child[1]->parent = nCopy;    
  }
  return nCopy;  
}


TmTreemap::~TmTreemap ()
{
  recursivelyDelete (root);
}

  
void TmTreemap::fullRecompute ()
{
  if (root!=NULL) root->resetFlagEverywhere (TmNode::IS_FEATURE_PASSED_VALID | TmNode::IS_GAUSSIAN_VALID);
  updateGaussians ();
  computeLinearEstimate ();  
}


void TmTreemap::recursiveAllNodes (TmNode* n, XycVector<TmNode*>& innerNode, XycVector<TmNode*>& leaf)
{
  if (n->isLeaf()) leaf.push_back (n);
  else {
    innerNode.push_back (n);
    recursiveAllNodes (n->child[0], innerNode, leaf);
    recursiveAllNodes (n->child[1], innerNode, leaf);
  }
}


void TmTreemap::optimalKLStep (TmNode* lca, double joinOnlyBelow, Move& move)
{
#if ASSERT_LEVEL>=3
    char treestruct1[10000], treestruct2[10000];
    lca->sprintTreeIndex (treestruct1);
#endif
  bool oldIsGaussianValidValid = isGaussianValidValid;  
  isGaussianValidValid = false;  
  lca->updateFeaturePassed ();
  
  move.clear();  
  if (!lca->isLeaf()) {
    recursiveOptimalKL (lca, lca->child[0], 0, joinOnlyBelow, move);
    recursiveOptimalKL (lca, lca->child[1], 1, joinOnlyBelow, move);  
  }
  updateFeaturePassed();  
  isGaussianValidValid =  oldIsGaussianValidValid;
#if ASSERT_LEVEL>=3
  assertIt();  
  Move dummy;
  safeOptimalKLStep (lca, joinOnlyBelow, dummy, move);  
  assert (move.cost==dummy.cost || fabs(move.cost-dummy.cost)<1E-6);
  lca->sprintTreeIndex (treestruct2);
  assert (strcmp(treestruct1, treestruct2)==0);    
#endif
}


void TmTreemap::recursiveOptimalKL (TmNode* lca, TmNode* subtreeBelow, int sideOfLca, double joinOnlyBelow, Move& bestMove)
{
  if (isIntersectionEmpty (subtreeBelow->featurePassed, lca->child[1-sideOfLca]->featurePassed)) return;

  if (subtreeBelow->isFlag(TmNode::CAN_BE_MOVED)) {    
    // try to move \c subtreeBelow to the other side
#if ASSERT_LEVEL>=3
    char treestruct1[10000], treestruct2[10000];
    lca->sprintTreeIndex (treestruct1);
#endif
    
    Move move (subtreeBelow, NULL); // Consider to move \c subtreeBelow
    if (subtreeBelow->parent!=lca) {
      // first move to the other side of lca
      subtreeBelow->moveTo (lca->child[1-sideOfLca], true);
      lca->updateFeaturePassed ();

      // Find the best place for it there
      recursiveOptimalDescend (lca->boundForChild (1-sideOfLca, bestMove.cost), lca->boundForChild (1-sideOfLca, joinOnlyBelow), move, true);
      move.cost = max(move.cost, lca->child[sideOfLca]->worstCaseUpdateCost) + lca->updateCost;

      // And move it back
      subtreeBelow->moveTo (move.oldAbove, true);
      subtreeBelow->makeChildNr (move.whichChild);
    }
    else {
      // We move one whole side of \c lca to the other. So just call recursiveOptimalDescend
      // But we do not allow subtree to stay where it is, because this would be a null-move.
      recursiveOptimalDescend (bestMove.cost, joinOnlyBelow, move, false);
    }    
    if (move.cost<bestMove.cost) bestMove = move;        
    
    lca->updateFeaturePassed ();      
  
#if ASSERT_LEVEL>=3
    lca->sprintTreeIndex (treestruct2);
    assert (strcmp(treestruct1, treestruct2)==0);    
#endif
  }
  
  assert (!bestMove.join || bestMove.cost<joinOnlyBelow+TmNode::costEps());  
  if (!subtreeBelow->isLeaf()) {    
    recursiveOptimalKL (lca, subtreeBelow->child[0], sideOfLca, joinOnlyBelow, bestMove);
    recursiveOptimalKL (lca, subtreeBelow->child[1], sideOfLca, joinOnlyBelow, bestMove);    
  }
  assert (!bestMove.join || bestMove.cost<joinOnlyBelow+TmNode::costEps());  
}


void TmTreemap::recursiveOptimalDescend (double bound, double joinOnlyBelow, Move& bestMove, bool mayStayHere)
{
  assert(bestMove.subtree!=NULL);
  TmNode* s = bestMove.subtree;
  bestMove.cost  = vmInf();
  bestMove.above = NULL;
  if (bound<0 || s->parent == NULL) return;  

  
  s->parent->updateFeaturePassed();
  TmNode* sibling;
  int whichChild;
  s->getSibling (sibling, whichChild);  

  if (isIntersectionEmpty (s->featurePassed, sibling->featurePassed)) return; // move only to somewhere across the border  

  if (mayStayHere) {    
    // let s stay where it is namely above \c sibling() 
    bestMove.above      = sibling;  
    bestMove.cost       = s->parent->worstCaseUpdateCost; // preliminary
    if (bestMove.cost<bound) bound = bestMove.cost;
  }

  if (s->parent!=NULL && s->isLeaf() && sibling->isLeaf()) {
    // Try whether \c subtree and \c above can be joined    
    double cost = JOINFACTOR*costOfJoining (s->parent); 
    if (cost<bestMove.cost && cost<joinOnlyBelow) {
      bestMove.above      = sibling;
      bestMove.join       = true;      
      bestMove.cost       = cost;
    }    
  }

  double lowerBound = s->worstCaseUpdateCost+TmNode::updateGaussianCost(s->featurePassed.size()); // lower bound for target function
  if (lowerBound>=bound) return;  

#if ASSERT_LEVEL>=3
    char treestruct1[10000], treestruct2[10000];
    s->parent->sprintTreeIndex (treestruct1);
#endif

  if (!sibling->isLeaf()) {
#ifndef NDEBUG
    double oldSiblingCost = sibling->worstCaseUpdateCost;
#endif
    Move bestAbove0 (bestMove), bestAbove1 (bestMove);

    // Check locations somewhere below sibling->child[0]
    s->moveTo (sibling->child[0], true);
    sibling->updateFeaturePassed();
    recursiveOptimalDescend (sibling->boundForChild (0, bound), sibling->boundForChild (0, joinOnlyBelow), bestAbove0, true);
    bestAbove0.cost = max(bestAbove0.cost, sibling->child[1]->worstCaseUpdateCost) + sibling->updateCost;
    if (bestAbove0.cost<bestMove.cost) {        
      bestMove  = bestAbove0;
      bound = bestMove.cost;
    }

    // Check locations somewhere below sibling->child[1]
    s->moveTo (sibling->child[1], true);
    sibling->updateFeaturePassed();
    recursiveOptimalDescend (sibling->boundForChild (1, bound), sibling->boundForChild (1, joinOnlyBelow), bestAbove1, true);
    bestAbove1.cost = max(bestAbove1.cost, sibling->child[0]->worstCaseUpdateCost) + sibling->updateCost;    
    if (bestAbove1.cost<bestMove.cost) {
      bestMove  = bestAbove1;
      bound = bestMove.cost;
    }

    // Move back
    s->moveTo (sibling, true);
    s->makeChildNr (whichChild);    
    sibling->parent->updateFeaturePassed();
    assert (oldSiblingCost == sibling->worstCaseUpdateCost);
  }

#if ASSERT_LEVEL>=3
  assert ((finite(bestMove.cost)) != (bestMove.above==NULL));
  assert (!bestMove.join || bestMove.cost<joinOnlyBelow);  
  s->parent->sprintTreeIndex (treestruct2);
  assert (strcmp(treestruct1, treestruct2)==0);    
#endif
}


void TmTreemap::safeOptimalKLStep (TmNode* lca, double joinOnlyBelow, Move& bestMove, const Move& specialMove)
{  
  
  XycVector<TmNode*> nodes;
  lca->nodesBelow (nodes);
  bestMove.clear();  
  
  for (int i=0; i<(int) nodes.size(); i++) 
    for (int j=0; j<(int) nodes.size(); j++) {
      updateFeaturePassed();      
      TmNode* s = nodes[i];
      TmNode* a = nodes[j];
      if (s==specialMove.subtree && a==specialMove.above) {
        // set a breakpoint here to stop when \c specialMove is considered
        assert (i>=0 && j>=0);        
      }      
      int whichSideS = s->whichSideOf (lca);
      int whichSideA = a->whichSideOf (lca);      
      if (s==lca || a==lca || s==a) continue;
      if (whichSideS>1 || whichSideA>1 || whichSideS==whichSideA) continue;      
      if (!s->isFlag(TmNode::CAN_BE_MOVED)) continue;
      if (isIntersectionEmpty (s->featurePassed, lca->child[whichSideA]->featurePassed)) continue;
      if (isIntersectionEmpty (a->featurePassed, s->featurePassed)) continue;

      TmNode* sibling;
      int whichChild;
      s->getSibling (sibling, whichChild);      
      TmNode* llca = lca;
      if (sibling!=a) { // no 'nop' moves, but 'nop' move plus join       
        if (s->parent==lca) llca = sibling;      

        s->moveTo (a, true);
        llca->updateFeaturePassed ();
        // effect of just moving there
        if (llca->worstCaseUpdateCost<bestMove.cost) { 
          bestMove.cost         = llca->worstCaseUpdateCost;
          bestMove.subtree      = s;
          bestMove.above        = a;
          bestMove.join         = false;        
        }
      }      
      if (s->parent!=NULL && s->isLeaf() && a->isLeaf()) { 
        // effect of joining \c s and \c a into one leaf
        llca->updateFeaturePassed ();        
        double cost = JOINFACTOR*costOfJoining (s->parent);
        assert (s->isAncestor (llca));      
        for (TmNode* n=s->parent; n!=llca; n=n->parent) {          
          cost = n->parent->updateCost + max(cost, n->sibling()->worstCaseUpdateCost);
        }        
        if (cost<bestMove.cost && cost<joinOnlyBelow) {
          bestMove.cost         = cost;
          bestMove.subtree      = s;
          bestMove.above        = a;
          bestMove.join         = true;          
        }        
      }      
      s->moveTo (sibling, true);
      s->makeChildNr (whichChild);      
    }
  updateFeaturePassed();
  bestMove.setOldAbove();
}


TmFeatureId TmTreemap::newFeatureBlock (int n)
{  
  int nn=n, id;
  while (nn<MAX_FEATURE_BLOCK_SIZE) {
    id = firstUnusedFeature[nn];    
    if (id>=0) {
      firstUnusedFeature[nn] = feature[id].nextUnusedFeature(); // discard block
      if (n<nn) { // but put remaining block into nn-n list
        feature[id+n].setNextUnusedFeature(firstUnusedFeature[nn-n]);
        firstUnusedFeature[nn-n] = id+n;
      }
      for (int j=id; j<id+n; j++) feature[j].setFlag (TmFeature::IS_EMPTY, 0);
      return id;      
    }
    nn += n;
  }
  id = feature.size();  
  feature.resize (id+n);
  for (int j=id; j<id+n; j++) feature[j].setFlag (TmFeature::IS_EMPTY, 0);
  return id;  
}

  
void TmTreemap::deleteFeature (TmFeatureId id)
{
  TmFeature& fId = feature[id];  
  fId.setFlag (TmFeature::IS_EMPTY, TmFeature::IS_EMPTY);
  fId.marginalizationNode = NULL;  
  // Try to join it to any existing block
  int i;  
  for (i=1; i<MAX_FEATURE_BLOCK_SIZE-1; i++) {
    int blockBase = id-i;
    if (blockBase>=0 && firstUnusedFeature[i]==blockBase) { // join this block and make it one longer
      // so remove it from firstUnusedFeature[i] and put it into firstUnusedFeature[i+1]
      firstUnusedFeature[i] = feature[blockBase].nextUnusedFeature();
      feature[blockBase].setNextUnusedFeature (firstUnusedFeature[i+1]);      
      firstUnusedFeature[i+1] = blockBase;
#if ASSERT_LEVEL>=3
  assertUnusedFeatureLists();  
#endif
      return;      
    }
  }
  // Else append it as a size 1 block
  fId.setNextUnusedFeature (firstUnusedFeature[1]);
  firstUnusedFeature[1] = id;
#if ASSERT_LEVEL>=3
  assertUnusedFeatureLists();  
#endif
}


void TmTreemap::printFeatureFragmentation () const
{
  printf ("feature.size()==%d feature.capacity()==%d\n", feature.size(), feature.capacity());
  for (int i=1; i<MAX_FEATURE_BLOCK_SIZE; i++) {
    int ctr=0, j=firstUnusedFeature[i];
    printf("size %d:", i);    
    while (j>=0) {
      if (ctr<10) printf ("%d ",j);      
      ctr++;
      j = feature[j].nextUnusedFeature();
    }    
    printf ("%5d free blocks\n", ctr);    
  }  
  printf("\n\n");  
}


void TmTreemap::hasBeenSparsifiedOut (TmFeatureId id)
{
  assert (feature[id].isFlag(TmFeature::CAN_BE_SPARSIFIED));  
}


TmNode* TmTreemap::addLeaf (const TmGaussian& gaussian, int flags)
{
  TmNode* leaf = new TmNode;
  leaf->index = -1;  
  leaf->gaussian = gaussian;
  leaf->status = (flags & TmNode::CAN_BE_INTEGRATED ) | TmNode::CAN_BE_MOVED;  
  addNonlinearLeaf (leaf);
  return leaf;  
}


void TmTreemap::newNodeIndex (TmNode* newNode)
{
  if (unusedNodes.empty()) {
    newNode->index = node.size();
    node.push_back (newNode);    
  }    
  else {
    newNode->index = unusedNodes.back();
    unusedNodes.pop_back();
    node[newNode->index] = newNode;      
  }
  stat.nrOfNodes++;  
}

  
void TmTreemap::addNonlinearLeaf (TmNode* newLeaf)
{
#if ASSERT_LEVEL >=3
  assertIt ();
#endif
#if ASSERT_LEVEL >=2
  assert (newLeaf->gaussian.R.isFinite());
#endif  
  newLeaf->tree = this;
  newLeaf->status &= ~(TmNode::IS_FEATURE_PASSED_VALID | TmNode::IS_GAUSSIAN_VALID);
  newLeaf->status |= TmNode::CAN_BE_MOVED | TmNode::IS_OPTIMIZED;
  if (newLeaf->index==-1) newNodeIndex (newLeaf);  
  isEstimateValid = false;
  if (root!=NULL) {
    // Make a new node root with newLeaf and root as children
    TmNode* n = new TmNode;
    newNodeIndex (n);    
    n->status = TmNode::CAN_BE_MOVED | TmNode::IS_OPTIMIZED;
       // We need \c IS_OPTIMIZED because otherwise \c setToBeOptimized
       // won't insert n into the queue
    n->tree = this;
    n->parent = NULL;
    n->child[0] = root;
    n->child[1] = newLeaf;
    newLeaf->parent = n;    
    root->parent = n;
    root = n;
    n->setToBeOptimized();
  }
  else {
    // just a single leaf
    root = newLeaf;    
    newLeaf->parent = NULL;    
  }  
  newLeaf->afterChange();
  updateFeaturePassed ();

  Move move;  
  optimalKLStep (root, root->worstCaseUpdateCost, move);
  if (!move.isEmpty() && move.cost<root->worstCaseUpdateCost) move.doIt ();
  updateFeaturePassed();
#if ASSERT_LEVEL >=2
  assertIt ();
#endif
}


void TmTreemap::clear()
{
  recursivelyDelete (root);  
  root = NULL;
  unusedNodes.clear();
  isEstimateValid = true;  
  isGaussianValidValid = true;
  feature.clear();
  node.clear();  
  for (int i=0; i<MAX_FEATURE_BLOCK_SIZE; i++) firstUnusedFeature[i]=-1;
  stat = TreemapStatistics();  
  workspace.clear();
  workspaceFloat.clear();  
}


void TmTreemap::optimizeFullRuns ()
{
  double oldAccumulatedOptimizationCost, factor=2;  
  if (root==NULL) return;  
  oldAccumulatedOptimizationCost = stat.accumulatedOptimizationCost;  
  while (updateGaussiansCost () + stat.accumulatedOptimizationCost - oldAccumulatedOptimizationCost
         <factor*root->worstCaseUpdateCost) {  //! TODO originally we had a do..while loop
    if (optimizer.optimizationQueue.empty()) return;    
    optimizer.oneKLRun ();
  } 
}


void TmTreemap::setInitialEstimate (TmFeatureId id, float est)
{
  if (!finite(feature[id].est)) feature[id].est = est;
}


void TmTreemap::updateFeaturePassed ()
{
  if (root!=NULL) root->updateFeaturePassed ();
}


void TmTreemap::onlyUpdateEstimatesFor (int from, int to, bool setDontUpdateFlag)
{
  updateFeaturePassed();  
  if (root==NULL) return;  
  if (setDontUpdateFlag) root->setFlagEverywhere (TmNode::DONT_UPDATE_ESTIMATE);
  for (int i=from; i<to; i++) {
    TmNode* n = feature[i].marginalizationNode;
    if (n!=NULL && n->isFlag (TmNode::DONT_UPDATE_ESTIMATE)) 
      n->resetFlagUpToRoot (TmNode::DONT_UPDATE_ESTIMATE);
  }
}


void TmTreemap::updateAllEstimates ()
{
  if (root!=NULL) root->resetFlagEverywhere (TmNode::DONT_UPDATE_ESTIMATE);
}


void TmTreemap::computeLinearEstimate ()
{
  if (isEstimateValid || root==NULL) return;  
  updateGaussians ();
  if (root->gaussian.RCompressed.empty()) root->estimate ();
  else root->estimateUsingRCompressed ();  
  isEstimateValid = true; 
#if ASSERT_LEVEL>=3
  assertEstimate ();
#endif
}


TmTreemap::SlamStatistic TmTreemap::slamStatistics () const
{
  return SlamStatistic();
}


void TmTreemap::updateGaussians ()
{
  if (root!=NULL) root->updateGaussian ();
}


void TmTreemap::computeNonlinearEstimate ()
{
  computeLinearEstimate ();
}


void TmTreemap::assertEstimate ()
{
  if (root==NULL) return;  
  computeLinearEstimate ();  
  TmExtendedFeatureList all;
  all.reserve (feature.size());  
  for (int i=0; i<(int) feature.size(); i++)
    if (feature[i].isDefined()) all.push_back(TmExtendedFeatureId (i,0));
  TmGaussian joined (all, root->rowsBelow());
  recursivelyMultiply (joined, root);
  joined.triangularize ();
  XymVector v(all.size());  
  joined.mean (v, all.size());
  double maxDiff = 0;  
  for (int i=0; i<(int) all.size(); i++) {
    double delta = feature[all[i].id].est-v[i];
    if (fabs(delta)>maxDiff) maxDiff = fabs(delta);    
    if (fabs(delta)>1E-3) {
      printf ("%d %f %f\n", i, feature[all[i].id].est, v[i]);
      assert (false);      
    }    
  }  
  if (maxDiff>1E-5) printf("Max error %e\n", maxDiff);  
}


void TmTreemap::computeStatistics ( TreemapStatistics& stat, bool expensive) const
{
  stat = this->stat;
  stat.nrOfNodesToBeOptimized = optimizer.optimizationQueue.size();
  if (expensive) stat.memory = memory ();
  else stat.memory = 0;  
}


void TmTreemap::computeEstimateByQR ()
{
  TmExtendedFeatureList all;
  all.reserve (feature.size());  
  for (int i=0; i<(int) feature.size(); i++)
    if (feature[i].isDefined()) all.push_back(TmExtendedFeatureId (i,0));
  int m;
  if (root!=NULL) m = root->rowsBelow();
  else m=0;  
  TmGaussian joined (all, m);
  recursivelyMultiply (joined, root);
  joined.triangularize ();
  XymVector v(all.size());  
  joined.mean (v, all.size());
  for (int i=0; i<(int) all.size(); i++) 
    feature[all[i].id].est = v[i];  
  isEstimateValid = true;  
}


double TmTreemap::costOfJoining (TmNode* subtree)
{
  subtree->updateFeaturePassed();  
  if (subtree->isLeaf()) return subtree->worstCaseUpdateCost;
  TmExtendedFeatureList fl;
  int nPM, nM, nP;  
  effectOfJoining (subtree, fl, nPM, nM, nP);
  if (fl.size()==0) return vmInf(); // not allowed to be joined  
  return TmNode::updateGaussianCost (nM+nP);
}


void TmTreemap::joinSubtree (TmNode* subtree)
{
  subtree->updateFeaturePassed (); 

  TmExtendedFeatureList fl;
  int nPM, nM, nP;  
  effectOfJoining (subtree, fl, nPM, nM, nP);
  assert (fl.size()>0);  // not allowed to be joined  
  
  // Now do it
  // collect all Gaussians
  for (int i=0; i<(int) fl.size(); i++) fl[i].count = 0;  
  TmGaussian joined (fl, subtree->rowsBelow());
  recursivelyMultiply (joined, subtree);  
  joined.triangularize();

  // Free features and adapt counter
  recursivelySubtractCount (subtree);  
  for (int i=0; i<nPM; i++) {
    // If there are no other copies it is marginalization otherwise sparsification.
    if (feature[fl[i].id].totalCount()==0) deleteFeature (fl[i].id); 
    else hasBeenSparsifiedOut (fl[i].id);    
  }  

  // Compute marginalized Gaussian
  joined.computeMarginal (subtree->gaussian, nPM);
  subtree->gaussian.compress ();  
#if ASSERT_LEVEL>=2
  subtree->gaussian.assertIt ();  
#endif
  for (int i=0; i<(int) subtree->gaussian.feature.size(); i++)
    feature[subtree->gaussian.feature[i].id].addTotalCount (subtree->gaussian.feature[i].count);  

  // Delete all children  
  recursivelyDelete (subtree->child[0]);
  recursivelyDelete (subtree->child[1]);
  subtree->child[0] = subtree->child[1] = NULL;  
  subtree->resetFlagUpToRoot (TmNode::IS_FEATURE_PASSED_VALID | TmNode::IS_GAUSSIAN_VALID);
  subtree->setFlag (TmNode::IS_OPTIMIZED | TmNode::CAN_BE_INTEGRATED);  
  for (TmNode* n=subtree->parent; n!=NULL; n=n->parent) n->setToBeOptimized ();  
  updateFeaturePassed ();
  assert ((int) subtree->gaussian.feature.size()==nP+nM);  

#if ASSERT_LEVEL>=2
  assertIt ();
#endif
}



bool TmTreemap::canBeSparsifiedOut (TmFeatureId id)
{
  return feature[id].isFlag (TmFeature::CAN_BE_SPARSIFIED);  
}


void TmTreemap::checkForSparsification (TmNode* n)
{
}


void TmTreemap::sparsifyOut (TmFeatureId id, int n)
{
  optimizer.report += " sparsified " + nameOfFeature(id) + " | ";  

  updateFeaturePassed ();  
  assert (0<=id && id+n<=(int) feature.size());  
  for (int i=0; i<n; i++) {
    assert (feature[id+i].isFlag (TmFeature::CAN_BE_MARGINALIZED_OUT));
    feature[id+i].setFlag (TmFeature::CAN_BE_SPARSIFIED, TmFeature::CAN_BE_SPARSIFIED);    
  }
  XycVector<TmNode*> node;
  findLeavesInvolving (id, node);
  // First check, whether \c CAN_BE_INTEGRATED is set
  for (int i=0; i<(int) node.size(); i++) 
    assert (node[i]->isFlag (TmNode::CAN_BE_INTEGRATED));  
  for (int i=0; i<(int) node.size(); i++)
    joinSubtree (node[i]);
  updateFeaturePassed ();
}


void TmTreemap::identifyFeatures (const XycVector<pair<int, int> >& assignment)
{
  updateFeaturePassed ();
  for (int i=0; i<(int) assignment.size(); i++) {
    int from = assignment[i].first;
    int to   = assignment[i].second;
    TmNode* m = feature[from].marginalizationNode;
    if (feature[to].marginalizationNode!=NULL)
      feature[to].marginalizationNode->resetFlagUpToRoot (TmNode::IS_FEATURE_PASSED_VALID | TmNode::IS_GAUSSIAN_VALID);    
    m->recursivelyIdentifyFeature (from, to);
    feature[to].addTotalCount (feature[from].totalCount());
    deleteFeature (from);    
  }  
}


void TmTreemap::recursivelySubtractCount (TmNode* n)
{
  if (n->isLeaf()) {
    const TmGaussian& g = n->gaussian;    
    for (int i=0; i<(int) g.feature.size(); i++)
      feature[g.feature[i].id].addTotalCount (-g.feature[i].count);
  }
  else {
    recursivelySubtractCount (n->child[0]);
    recursivelySubtractCount (n->child[1]);
  }
}


void TmTreemap::recursivelyDelete (TmNode* n)
{
  if (n!=NULL) {    
    if (!n->isLeaf()) {
      recursivelyDelete (n->child[0]);
      recursivelyDelete (n->child[1]);
    }
    node[n->index] = NULL;
    unusedNodes.push_back (n->index);    
    stat.nrOfNodes--;    
    delete n;
  }  
}


void TmTreemap::recursivelyMultiply (TmGaussian& join, TmNode* subtree)
{
  if (subtree->isLeaf()) {
    join.multiply (subtree->gaussian, 0); 
    // TODO rotate
  }
  else {
    recursivelyMultiply (join, subtree->child[0]);
    recursivelyMultiply (join, subtree->child[1]);
  }
}


void TmTreemap::effectOfJoining (TmNode* subtree, TmExtendedFeatureList& fl, int& nPM, int& nM, int& nP) const
{
  fl.clear();
  if (!subtree->isFlag(TmNode::CAN_BE_INTEGRATED)) {
    nPM = nM = nP = 0;
    return;
  }  

  // Build list of all features involved in this subtree
  TmExtendedFeatureList flx;  
  recursivelyAdd (flx, subtree);
  sumUp (flx);

  // First count
  nPM = 0; // how many can be marginalized out permanently (including those that are sparsified out)
  nM  = 0; // how many will be marginalized out at subtree
  nP  = 0; // how many will be passed to the parent
  for (int i=0; i<(int) flx.size(); i++) {
    const TmExtendedFeatureId& feat = flx[i];
    const TmFeature& feat2 = feature[feat.id];    
    if (feat.id==subtree->linearizationPointFeature) 
      nP++; // We must keep and pass linearization point features    
    else if (feat.count < feat2.totalCount()) {
      // Also involved outside subtree so we will have to pass it unless
      if (feat2.isFlag (TmFeature::CAN_BE_SPARSIFIED)) nPM++; // except may be when it can be sparsified out
      else nP++;      
    }
    else {
      // Only involved inside, so we might marginalize out
      if (feat2.isFlag (TmFeature::CAN_BE_MARGINALIZED_OUT)) nPM++; // permanently      
      else nM++; // or in subtree      
    }
  }

  // Then sort into fl making the same decisions again
  fl.resize (nPM+nM+nP);
  int iPM=0, iM=nPM, iP=nPM+nM;  // ctr
  for (int i=0; i<(int) flx.size(); i++) {
    const TmExtendedFeatureId& feat = flx[i];
    const TmFeature& feat2 = feature[feat.id];    
    if (feat.id==subtree->linearizationPointFeature) {
      fl[iP] = feat;
      iP++;
    }
    else if (feat.count < feat2.totalCount()) {
      if (feat2.isFlag (TmFeature::CAN_BE_SPARSIFIED)) {
        fl[iPM] = feat;
        iPM++;
      }
      else {
        fl[iP] = feat;
        iP++;
      }      
    }
    else {
      if (feat2.isFlag (TmFeature::CAN_BE_MARGINALIZED_OUT)) {
        fl[iPM] = feat;
        iPM++;
      }      
      else {
        fl[iM] = feat;
        iM++;
      }      
    }
  }
  assert (iPM==nPM && iM==nPM+nM && iP==nPM+nM+nP && iP==(int) flx.size());  
}


void TmTreemap::recursivelyAdd (TmExtendedFeatureList& fl, TmNode* subtree) const
{
  if (subtree->isLeaf()) {
    // add only features from leaves
    for (int i=0; i<(int) subtree->gaussian.feature.size(); i++)
      fl.push_back (subtree->gaussian.feature[i]);
  }
  else {
    recursivelyAdd (fl, subtree->child[0]);
    recursivelyAdd (fl, subtree->child[1]);
  }  
}


void TmTreemap::computeFeaturesInvolvedBelow (TmNode* subtree, TmExtendedFeatureList& fl)
{
  fl.clear();
  recursivelyAdd (fl, subtree);
  sumUp (fl);  
}


void TmTreemap::recursivelyCount (TmNode* n, XycVector<int>& count) const
{
  if (n==NULL) return;  
  if (n->isLeaf()) {
    for (int i=0; i<(int) n->gaussian.feature.size(); i++) {
      int ct = n->gaussian.feature[i].count;      
      count[n->gaussian.feature[i].id] += ct;
    }    
  }
  else {
    recursivelyCount (n->child[0], count);
    recursivelyCount (n->child[1], count);
  }  
}


TmNode* TmTreemap::lca (TmNode* a, TmNode* b) const
{
    if (a==NULL) return b;
    if (b==NULL) return a;

    // Count the level in which a and b are located
    int levelA=0;
    for (TmNode* n=a; n!=NULL; n=n->parent) levelA++;
    int levelB=0;
    for (TmNode* n=b; n!=NULL; n=n->parent) levelB++;

    // Move the lower one up to the same level
    while (levelA>levelB) {
      a = a->parent;
      levelA--;
    }
    while (levelB>levelA) {
      b = b->parent;
      levelB--;
    }

    // Move upward until they meet
    while (a!=b) {
      a=a->parent;
      b=b->parent;   
    }
    return a;
}


void TmTreemap::rotateGaussian (TmGaussian& gaussian, double angle) const
{
  assert (false); // not implemented
}


void TmTreemap::assertUnusedFeatureLists () const
{
  set<int> unused;  
  for (int i=0; i<feature.size(); i++) {
    if (feature[i].isEmpty()) unused.insert (i);
  }  
  for (int bs=1; bs<MAX_FEATURE_BLOCK_SIZE; bs++) {
    int i=firstUnusedFeature[bs];
    while (i>=0) {
      for (int j=i; j<i+bs; j++) {
        assert (feature[j].isEmpty());
        set<int>::iterator it = unused.find(j);        
        assert (it!=unused.end());
        unused.erase (it);        
      }
      i = feature[i].nextUnusedFeature();      
    }  
  }
  assert (unused.empty());  
}


void TmTreemap::assertIt () const
{
  assert (root==NULL || !root->isFlag(TmNode::IS_FEATURE_PASSED_VALID) ||
          root->featurePassed.empty());  
  XycVector<int> count (feature.size(), 0);
  recursivelyCount (root, count);
  for (int i=0; i<(int) feature.size(); i++) {
    assert (count[i]==feature[i].totalCount());
    const TmFeature& feat = feature[i];    
    TmNode* m = feat.marginalizationNode;
    if (m!=NULL) {
      assert (feat.totalCount()>0);      
      assert (m->tree==this);
      m->assertInTree();      
      assert (getNode(m->index) == m);
    }
  }
  int ctr=0;
  for (int i=0; i<(int) node.size(); i++) {
    TmNode* n = node[i];
    // Check that it is in the tree
    if (n!=NULL) {
      assert (n->tree==this);      
      n->assertInTree ();   
      ctr++;      
    }
    else {
      int j;
      for (j=0;j<unusedNodes.size();j++) if (unusedNodes[j]==i) break;
      assert (j<unusedNodes.size());
    }    
      
    bool isOptimized;    
    if (n!=NULL) {
      assert (n->index==i);
      isOptimized = n->isFlag (TmNode::IS_OPTIMIZED);
    }
    else isOptimized = true;
    int ctr =0;    
    for (int j=0; j<(int) optimizer.optimizationQueue.size(); j++)
      if (optimizer.optimizationQueue[j]==i) ctr++;
    assert (ctr>=(isOptimized?0:1));
  }  
  assert (stat.nrOfNodes==ctr);  
  assertUnusedFeatureLists ();  
#if ASSERT_LEVEL>=3
  assertConnectivity (3); // 3 is only o.K. for planar maps  
#endif
  if (root!=NULL) root->recursiveAssertIt ();
}


void TmTreemap::findLeavesInvolving (TmFeatureId id, XycVector<TmNode*>& node)
{
  updateFeaturePassed();
  node.clear();  
  int ctr=0;  
  if (0<=id && id<(int) feature.size()) {
    TmNode* lca = feature[id].marginalizationNode;
    if (lca!=NULL) lca->recursiveAddLeavesInvolving (id, node, ctr);
    assert (ctr==feature[id].totalCount());    
  }
}




double TmTreemap::time() 
{
#ifdef linux
    static bool isFirst = true;
    static timeval tv_base;
#ifdef USEPROCESSTIME
    timeval tv;
    rusage ruse;
    getrusage (RUSAGE_SELF, &ruse);
//    gettimeofday (&tv, NULL);
    tv = ruse.ru_utime;
#else
    timeval tv;
    gettimeofday (&tv, NULL);
#endif
    if (isFirst) {
        tv_base = tv;
        isFirst = false;
    }
    return (tv.tv_sec-tv_base.tv_sec)
        +1E-6*(tv.tv_usec-tv_base.tv_usec);
#else
#error no linux
    return 0;
#endif
}


void TmTreemap::calibrateGaussianPerformance (int nMax, double coef[4], int deactivate, char* filename)
{
/*  XycVector<double> data;  
  TmGaussian fitGaussian ("a b c d");
  int startAt = 4;  
  for (int i=startAt; i<nMax; i++) {
    TmFeatureList fl0, fl1, fl;
    for (int j=0;j<i;j++) fl.push_back(j);    
    for (int j=0;j<2*i/3;j++) fl0.push_back(j);
    for (int j=i/3; j<i; j++) fl1.push_back(j);
    TmGaussian A = TmGaussian::relative1D (fl0, true);
    TmGaussian B = TmGaussian::relative1D (fl1);    
    double t0 = time();    
    TmGaussian test (fl, A.R.rows()+B.R.rows());
    test.multiply (A, 0);
    test.multiply (B, 0);
    test.triangularize();    
    double t1 = time();
    data.push_back (t1-t0);    
    fitGaussian.R.appendRow (1,true);
    double f = 1/(t1-t0);
    int row = fitGaussian.R.rows()-1;    
    if ((1&deactivate)==0) fitGaussian.R (row,0) = f*1;
    if ((2&deactivate)==0) fitGaussian.R (row,1) = f*i;
    if ((4&deactivate)==0) fitGaussian.R (row,2) = f*i*i;
    if ((8&deactivate)==0) fitGaussian.R (row,3) = f*i*i*i;
    fitGaussian.R (row,4) = -f*(t1-t0);
  }
  if ((1&deactivate)!=0) {
    fitGaussian.R.appendRow (1,true);
    fitGaussian.R(fitGaussian.R.rows()-1,0)=1;
  }
  if ((2&deactivate)!=0) {
    fitGaussian.R.appendRow (1,true);
    fitGaussian.R(fitGaussian.R.rows()-1,1)=1;
  }
  if ((4&deactivate)!=0) {
    fitGaussian.R.appendRow (1,true);
    fitGaussian.R(fitGaussian.R.rows()-1,2)=1;
  }
  if ((8&deactivate)!=0) {
    fitGaussian.R.appendRow (1,true);
    fitGaussian.R(fitGaussian.R.rows()-1,3)=1;
  }  
  fitGaussian.triangularize();
  XymVector v (4);  
  fitGaussian.mean (v);
  for (int i=0; i<4; i++) coef[i] = v[i];  
  FILE* file;
  if (filename!=NULL) {
    file = fopen (filename, "w");  
    for (int i=0; i<(int) data.size(); i++) {
      int n = i+startAt;      
      fprintf (file, "%d %e %e\n", n, data[i], coef[0] + coef[1]*n + coef[2]*n*n + coef[3]*n*n*n);    
    }    
    fclose(file);  
  }  
*/
}



bool TmTreemap::isCompiledWithOptimization()
{
#if defined(NDEBUG)
  return xymIsCompiledWithOptimization();
#else
  return false;  
#endif
}

int TmTreemap::nrOfFeatures (int mustFlag, int mayNotFlag) const
{
  int ctr=0;
  for (int i=0; i<(int) feature.size(); i++) if (feature[i].isDefined()) ctr++;
  return ctr;
}


int TmTreemap::memory () const
{
  int mem = sizeof(TmTreemap);
  mem += node.capacity() * sizeof(TmNode);
  mem += unusedNodes.capacity() * sizeof(int);
  mem += feature.capacity() * sizeof(TmFeature);
  mem += optimizer.memory() - sizeof(Optimizer);
  mem += workspace.memoryUsage();
  mem += workspaceFloat.capacity() * sizeof(float);  
  if (root!=NULL) mem += root->recursiveMemory ();  
  return mem;  
}


void TmTreemap::nameOfFeature (char* txt, int featureId, int& n) const
{
  featureToString (txt, featureId);
  n = 1;  
}

string TmTreemap::nameOfFeature (int featureId)
{
  char txt[10];
  int n;  
  nameOfFeature (txt, featureId, n);
  return string(txt);  
}


double TmTreemap::updateGaussiansCost () const
{
  return recursiveUpdateGaussiansCost (root);
}


double TmTreemap::recursiveUpdateGaussiansCost (const TmNode* n) const
{
  if (n==NULL) return 0;
  else if (n->isFlag(TmNode::IS_GAUSSIAN_VALID)) return 0;
  else if (n->isLeaf()) return n->updateCost;
  else return n->updateCost + n->child[0]->updateCost + n->child[1]->updateCost;
}


void TmTreemap::TreemapStatistics::optimizationStatistics (int n, bool success)
{
  if (n>=(int) htp.size()) htp.resize (n+1);  
  if (success) htp[n].success++;  
  else htp[n].noSuccess++;  
}


double TmTreemap::TreemapStatistics::optimizationCondProb (int n)
{
  if (n>=(int) htp.size()) return 0;
  int nPos=0, nNeg=0;
  for (int i=n+1; i<(int) htp.size(); i++) {
    nPos += htp[i].success;  
    nNeg += htp[i].noSuccess;
  }
  if (nPos+nNeg>0) return double(nPos)/(nPos+nNeg);  
  else return 0;  
}


void TmTreemap::assertConnectivity (int minDOF) const
{
  // Initialize components
  XycVector <int> component;
  component.resize (node.size());
  int ctr=0;  
  for (int i=0; i<(int) node.size(); i++) {
    if (node[i]!=NULL && node[i]->isLeaf()) {
      component[i]=ctr;
      ctr++;
    }    
    else component[i]=-1;
  }
  // Compute connected components
  for (int i=0; i<(int) component.size(); i++) if (component[i]>=0) {
    for (int j=0; j<i; j++) if (component[j]>=0 && component[j]!=component[i]) {
      // Check, whether node[i] and node[j] share at least minDOF features
      int sharedDOF = nrOfIntersecting (node[i]->gaussian.feature, node[j]->gaussian.feature);
      if (sharedDOF>=minDOF) {
        // Join the components
        int oldC = component[j];
        int newC = component[i];
        if (oldC<newC) {
          int swap = newC;
          newC = oldC;
          oldC = swap;
        }        
        for (int k=0; k<=i; k++) if (component[k]==oldC) component[k]=newC;
      }      
    }    
  }
  // Assert that there is only one
  for (int i=0; i<(int) component.size(); i++) assert (component[i]<=0); // Either component 0 or -1  
}


void TmTreemap::printGaussian (TmGaussian& g)
{
  int i=0;
  TmExtendedFeatureList& fl = g.feature;  
  while (i<(int) fl.size()) {
    int n;
    char txt[10];
    nameOfFeature (txt, fl[i].id, n);
    for (int j=0; j<n; j++) {
      assert (fl[i+j].id==fl[i].id+j);
      printf ("%4s[%1d],", txt, j);
    }    
    i += n;
  }
  printf ("\n");
  XymMatrixVC& R = g.R;  
  for (int i=0; i<R.rows(); i++) {
    for (int j=0; j<R.cols(); j++) {
      if (R(i,j)==0) printf ("       ,");
      else printf ("%7.1f,",R(i,j));
    }
    printf("\n");
  }  
}


//*********** TmTreemap::Move
void TmTreemap::Move::clear ()
{
  subtree = oldAbove = above = NULL;
  whichChild = -1;
  cost = vmInf();
  join = false;        
}


void TmTreemap::Move::setOldAbove()
{
  if (subtree==NULL || subtree->parent==NULL) {
    oldAbove = NULL;
    whichChild = -1;
  }
  else{          
    TmNode* p = subtree->parent;        
    if (p->child[0]==subtree) {
      oldAbove = p->child[1];
      whichChild = 0;
    }
    else {
      oldAbove = p->child[0];          
      whichChild = 1;
    }        
  }        
}  


void TmTreemap::Move::tryIt ()    
{
  assert (!join && subtree->isFlag(TmNode::CAN_BE_MOVED));        
  setOldAbove ();
  subtree->moveTo (above, true);
  subtree->resetFlag (TmNode::CAN_BE_MOVED);  
}    


void TmTreemap::Move::doIt ()
{ 
  assert (subtree->isFlag (TmNode::CAN_BE_MOVED));  
  setOldAbove ();
  subtree->moveTo (above, false);
  if (join) subtree->tree->joinSubtree (subtree->parent);
}    


void TmTreemap::Move::undoIt ()
{ 
  assert (above->parent==subtree->parent && !join); 
  assert (!subtree->isFlag (TmNode::CAN_BE_MOVED));  
  subtree->moveTo (oldAbove, true);
  subtree->makeChildNr (whichChild);
  subtree->setFlag (TmNode::CAN_BE_MOVED);  
}


/****** TmTreemap::Optimizer *********/

TmTreemap::Optimizer::Optimizer ()
  :tree (NULL), optimizationQueue(), lcaIndex(-1), initialCost(0), unsuccessfulMoves(), 
   maxNrOfUnsuccessfulMoves(-1), nrOfMovesPerStep(-1), report()
{}


TmTreemap::Optimizer::Optimizer (TmTreemap* tree, int nrOfMovesPerStep, int maxNrOfUnsuccessfulMoves)
  :tree (NULL), optimizationQueue(), lcaIndex(-1), initialCost(0), unsuccessfulMoves(), 
   maxNrOfUnsuccessfulMoves(-1), nrOfMovesPerStep(-1), report()
{
  create (tree, nrOfMovesPerStep, maxNrOfUnsuccessfulMoves);  
}


TmNode* TmTreemap::Optimizer::nextNodeToBeOptimized ()
{
  TmNode* n = NULL;  
  while (!optimizationQueue.empty() && n==NULL) {    
    if (lcaIndex<0) { // Fetch new node
      lcaIndex = optimizationQueue.front();
      n = tree->getNode (lcaIndex);
      if (n!=NULL) {
        initialCost = n->worstCaseUpdateCost;
      }      
    }    
    else n = tree->getNode (lcaIndex); // Just get ptr to old node    
    if (n==NULL) { // Node does not exist any more
      optimizationQueue.pop_front();
      lcaIndex = -1;
    }
  }
  return n;  
}
  

void TmTreemap::Optimizer::create (TmTreemap* tree, int nrOfMovesPerStep, int maxNrOfUnsuccessfulMoves)
{
  this->tree = tree;
  optimizationQueue.clear();
  lcaIndex = -1;  
  initialCost = 0;
  unsuccessfulMoves.clear();
  this->nrOfMovesPerStep         = nrOfMovesPerStep;  
  this->maxNrOfUnsuccessfulMoves = maxNrOfUnsuccessfulMoves;
  report = string();  
}

      
string TmTreemap::Optimizer::getAndClearReport ()
{
  string s = report;
  report = string();
  return s;
}


int TmTreemap::Optimizer::memory () const
{
  return sizeof(Optimizer) + optimizationQueue.size()*sizeof(int) + unsuccessfulMoves.capacity()*sizeof(MoveIndices);
}


void TmTreemap::Optimizer::oneKLRun ()
{
  if (optimizationQueue.empty()) return;  
  TmNode* lca = nextNodeToBeOptimized ();
  if (lca==NULL) return;  

#ifndef NDEBUG
  char txt[100];
#endif
  XycVector<Move> moves;
  bool oldIsGaussianValidValid = tree->isGaussianValidValid;
  tree->isGaussianValidValid = false;

#ifndef NDEBUG
  sprintf (txt, "optimizing %d ", lca->index);
  report += txt;  
#endif
  
  double initialCost = lca->worstCaseUpdateCost;
  double bestCost = initialCost;
  TmNode* currentLca = lca;  
  Move move;
  bool didSomething = false;  
  while ((int) moves.size()<maxNrOfUnsuccessfulMoves) {
    // Do greedy moves preliminary even if they increase worstCaseUpdateCost
    tree->optimalKLStep (currentLca, bestCost-TmNode::costEps(), move);
    assert (!move.join || move.cost<bestCost+TmNode::costEps());
    if (move.isEmpty()) break;
    moves.push_back (move);    
    if (move.subtree->parent==currentLca && move.above->parent!=currentLca) {
      // We move one complete side of lca so lca is not the lca any more.
      currentLca = move.oldAbove;
    }
    if (move.cost<bestCost-TmNode::costEps()/2) {
      // Really execute everything
      // Undo all changes 
      for (int i=(int) moves.size()-2; i>=0; i--) { // The last one has not been done yet
        moves[i].undoIt ();
        moves[i].subtree ->setFlag (TmNode::CAN_BE_MOVED);
      }
      // Redo the changes that let to an improvement permanently
      tree->stat.optimizationStatistics (moves.size(), true);      
      for (int i=0; i<(int) moves.size(); i++) {
#ifndef NDEBUG
        sprintf (txt, "%d --> %d, ", moves[i].subtree->index, moves[i].above->index);
        report += txt;        
        if (moves[i].join) report += " joined ";
#endif
        moves[i].doIt ();
      }
      bestCost = moves.back().cost;
      moves.clear();
      didSomething = true;      
      break; 
    }
    else {
      // Just try whether this will lead to something
      move.tryIt ();
      move.subtree->resetFlag (TmNode::CAN_BE_MOVED);
    }
    currentLca->updateFeaturePassed();    
    assert (move.join || fabs(move.cost-currentLca->worstCaseUpdateCost)<TmNode::costEps());
//    assert (move.cost==currentLca->worstCaseUpdateCost);
  }
  tree->stat.optimizationStatistics (moves.size(), false);  
  for (int i=moves.size()-1; i>=0; i--) {
    moves[i].undoIt ();
    moves[i].subtree ->setFlag (TmNode::CAN_BE_MOVED);
  }  
  if (didSomething) {
    optimizationQueue.push_back (optimizationQueue.front());  
#ifndef NDEBUG
    sprintf (txt, " (%6.4fms >> %6.4fms) ", initialCost*1000, bestCost*1000);    
    report += txt;    
#endif
  }
  else {
    TmNode* n = tree->getNode (optimizationQueue.front());
    lcaIndex = -1;    
    if (n!=NULL) n->setFlag (TmNode::IS_OPTIMIZED);
#ifndef NDEBUG
    report += " optimal | ";    
#endif
    if (n!=NULL) tree->checkForSparsification (n);
  }
  optimizationQueue.pop_front();
  tree->updateFeaturePassed();
  tree->isGaussianValidValid = oldIsGaussianValidValid;  
#if ASSERT_LEVEL>=2
  tree->assertIt();
  if (tree->root!=NULL) tree->root->recursiveAssertFlag (TmNode::CAN_BE_MOVED);  
#endif
}


