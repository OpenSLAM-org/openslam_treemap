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
#ifndef TMNODE_H
#define TMNODE_H

/*!\file tmNode.h 
   \brief Class \c TmNode
   \author Udo Frese

  Contains the class \c TmNode representing a
  node in the treemap tree. 
*/

#include "tmTypes.h"
#include "tmGaussian.h"
#include "tmExtendedFeatureId.h"
#include <limits.h>


class TmTreemap;

//! A node in the treemap tree.
/*! 
 */
class TmNode 
{
 public:
  //! A \c TmNode makes only sense in a TmTreemap
  friend class TmTreemap;
  friend class MainWindow;  

  //! An uninitialised node
  TmNode ();

  //! Virtual destructor for overloading nodes.
  /*! Overloaded nodes are used to implement non-Gaussian
      distributions which are linearized, i.e. approximated by a
      Gaussian periodically.
  */
  virtual ~TmNode();


  //! Returns a polymorphic shallow copy of \c this
  /*! Must be overloaded by any derived class. Calls the copy
      constructor of the same class.
  */
  virtual TmNode* duplicate () const;  
  

  //! A consecutive index for identifying nodes
  int index;  
  
  //! Pointer to the tree in which the node resides.
  TmTreemap* tree;  
  
  //! Pointer to the parent node in the tree
  /*! \c NULL in the root node. */
  TmNode *parent;
  //! Pointer to both children
  /*! The children are index by 0 and 1 and there is no conceptual
      difference between child 0 and 1
  */
  TmNode *child[2];
  
  //! Cost for updating this node from its children
  double updateCost;  

  //! Largest possible cost for updating any descendant up to this node
  /*! This cost is used as a goal function for optimizing both
      partitioning and balancing of the tree. The cost for a path is
      defined as the sum of \c updateGaussianCost of the number of landmarks in \c
      featuresPassed in each node. The worst case cost is the largest
      number for all pathes from leaves up to this node.
  */
  double worstCaseUpdateCost;

  //! List of features passed to the parent node.
  /*! The list is sorted by ascending feature ids. Apart from the
      feature index it contains a counter that is added when
      information from both children on the same feature is
      integrated. It is compared to the overall counter \c
      tree->.. for that feature. If both are equal then every
      information about that feature has been collected and the
      feature can be marginalized out in this node.
  */
  TmExtendedFeatureList featurePassed;

  //! See \c gaussian.linearizationPointFeature
  /*! We have to store the linearizationPointFeature also in the node,
      because it must be updated not only when updating the Gaussians
      but also when updating the \c featuresPassed list.
   */
  int linearizationPointFeature;  


  //! The gaussian distribution stored in this node.
  /*! The Gaussian resulting from integrating, i.e. multiplying the
      Gaussians passed by the two children. For a leaf it is the
      original distribution being input to the treemap. This
      distribution is the integration, i.e.  product of all original
      distributions stored at leafs below this node but marginalizing
      landmarks on the way that need not be passed to the parent.

      The gaussian is stored in triangular form.

      The gaussian involves both features that are passed to the
      parent and features that are marginalized out at this node. Due
      to its triangular form the marginalized Gaussian is described by
      a submatrix involving only rows and columns \c
      [firstFeaturePassed..gaussian.n]. Since the columns of \c
      gaussian.R must be grouped into marginalized
      vs. not-marginalized, the \c gaussian.feature is not
      sorted. However \c featuresPassed is sorted and \c
      featuresPassed[i] corresponds to \c
      gaussian.feature[firstFeaturePassed+i] as well as column \c
      firstFeaturePassed+i.

      \warning The application may not change \c gaussian by itself but must
      call \c changeGaussian instead.
   */
  TmGaussian gaussian;

  //! Features [firstFeaturePassed..] in \c gaussian are passed to the parent.
  /*! See \c gaussian. */
  int firstFeaturePassed;

  //! Different status bits
  /*! All validity flags follow the flow of information in the
      tree. If something is invalid at a node then it is invalid at
      all ancestor nodes too.
  */
  enum NodeFlags {IS_FEATURE_PASSED_VALID=1, IS_GAUSSIAN_VALID=2, 
                  IS_OPTIMIZED=4, DONT_UPDATE_ESTIMATE=8,
                  CAN_BE_MOVED=16, CAN_BE_INTEGRATED=32};  

  /*! \var IS_FEATURE_PASSED_VALID 
   
      whether \c featuresPassed is up to date.
   */
  /*! \var IS_GAUSSIAN_VALID 
     
      whether \c gaussian is up to date implies that
      \c IS_FEATURE_PASSED_VALID is true. There is a special case in the
      optimization subalgorithm (\c TmTreemap::Optimizer::oneKLRun).
  */
  /*! \var IS_OPTIMIZED 

      is set when the HTP optimization algorithm found it 
      could not improve that node any more. It is reset whenever some node
      enters or leaves the node's two children's subtrees. If it is not set, the
      node is in the optimization queue \c TmTreemap::Optimizer::optimizationQueue
  */
  /*! \var CAN_BE_MOVED 
 
      Flag used in the HTP optimization subalgorithm
      for a Kernighan-Lin run. There, after a node is moved, \c CAN_BE_MOVED
      is set to \c false so the same node is not moved during that KL
      run any more. Outside the HTP algorithm it is always \c true.
   */
  /*! \var DONT_UPDATE_ESTIMATE 

      indicates that the estimate for the landmarks
      marginalized at this node should be not recomputed in the estimation step.
      In theory all estimates must be updated even after one Gaussian changes.
      However if one is not interested in certain estimates, e.g. because they
      correspond to remote parts of the map. One can select to ommit these
      from the computation. Since the update time per landmark is very small
      this makes only sense for extremely large map (say >=20000).

      If both children are \c DONT_UPDATE_ESTIMATE, the parent sets
      \c DONT_UPDATE_ESTIMATE too.
  */
  /*! \var CAN_BE_INTEGRATED 

      whether the algorithm is allowed to integrate
      the Gaussian represented by this node into another Gaussian and then
      get completely rid of this node. If this is done, the Gaussian is fixed
      forever and cannot be recomputed with new linearizations points. 
  */

      
  //! Different status bits or'ed
  int status;

  //! Returnst whether the flag passed is set
  bool isFlag (int flag) const {return (status & flag) !=0;}  

  //! Sets the flag(s) passed in \c flag to true.
  void setFlag (int flag) {status |= flag;}

  //! Resets the flag(s) passed in \c flag to false
  void resetFlag (int flag) {status &= ~flag;}

  //! Resets the flag(s) passed in \c flag from \c this up to the root
  /*! The routine stops when it encounters a node where all these flags
      have already been reset. Thus computation time is bounded by the
      number of nodes actually reset. This is only valid for flags that
      imply that the same flag is set at the parent.

      If \c flag includes \c IS_GAUSSIAN_VALID, \c tree->isEstimateValid
      is also set to \c false.
  */
  inline void resetFlagUpToRoot (int flag);  

  //! Resets the flags \c flag in all nodes below this node
  void resetFlagEverywhere (int flag);  

  //! Sets the flags \c flag in all nodes below this node
  /*! If \c stopIfFound, the routine stops recursion when it finds a
      node that is already set.  Then it can only be used for flags
      where being set at a node implies being set at the children.
  */
  void setFlagEverywhere (int flag, bool stopIfFound=true);  


  //! Indicates that this node has to be optimized later on.
  /*! If \c IS_OPTIMIZED is not set nothing is done. If it is set, it is
      reset and \c this added to the end of \c tree->optimizationQueue.
   */
  void setToBeOptimized ();

  //! Calls \c setToBeOptimized from \c this up to the root
  /*! This routine can be called even for leaves but does not set the
    leaves flag. */
  void setToBeOptimizedUpToRoot ();
  




  //! Recursively updates \c .featurePassed and assorted as far as necessary
  /*!
      Specifically updates \c .worstCaseUpdateCost, \c .featurePassed,
      \c .status, \c .linearizationPointFeature. Does not update \c .gaussian. In a
      leaf it calls \c updateGaussianOfALeaf().
   */
  virtual void updateFeaturePassed ();

  //! Recursively updates \c .gaussian as far as necessary
  /*!
      Implicitly updates \c .featurePassed. 
      For an inner node, the Gaussians passed by both children are 
      integrated (multiplied) and the result is transformed into
      triangular form and stored in \c .gaussian. \c .gaussian contains
      both the part passed to the parent and the part marginalized out
      in this node. For a leaf the .gaussian is retransformed into
      triangular form. This is necessary since it may have changed,
      which features are eliminated here.

      If both of childrens defines a \c .linearizationPointFeature,
      the gaussian of the second child is exactly rotated by the
      estimated (\c TmTreemap::feature) difference between both
      orientations before integrating and the resulting Gaussian uses
      the first \c .linearizationPointFeature as its own feature. This
      way the relative orientation is fixed in the computation above
      this node but the absolute orientation can still be changed. If
      only one child defines a \c .linearizationPointFeature, it is
      rotated according to the estimate of that feature. The resulting
      node then defines no \c .linearizationPointFeature and cannot be
      rotated further. This rotation is performed by \c
      TmTreemap::rotateGaussian.
  */
  void updateGaussian ();

  //! Recursively estimates all features marginalized out at or below this node.
  /*! The estimate (\c tree->feature) for all features in \c
      featuresPassed must already be computed. 
   */
  void estimate ();  

  //! Same as \c estimate but uses the optimized representation in \c TmGaussian::RCompressed 
  void estimateUsingRCompressed ();
  

  //! Changes the original distribution of a leaf to \c gaussian.
  /*! Nodes are invalidated accordingly. Note, the application must
      use this routine or \c addLeaf or \c addNonlinearLeaf when
      changing an original distribution because otherwise nodes will
      not be invalidated properly. */
  void changeGaussian (const TmGaussian& gaussian);

  
  //! Decrements \c TmFeature::totalCount() for all features involved at this leaf
  /*! When changing the Gaussian at a leaf, first call \c beforeChange, then change
      the Gaussian, then call \c afterChange. Can only be called for leaves.
  */
  void beforeChange ();  


  //! Invalidate everything after the gaussian has changed
  /*! Can only be called for leaves.
 
      Invalidates this node, all marginalization nodes of features
      involved and their ancestors. Also copies the linearization
      point into the estimate of the linearization point feature if
      there is no estimate available.

      Increments the \c TmFeature::totalCount() counter for each
      feature involved in this node. If the node has been changed,
      i.e. it has been in the tree before, \c beforeChange() must have
      been called so all counters of the features present before the
      change have been decremented.

      Must be called after inserting a new node (already done by \c
      TmTreemap::addLeaf and \c TmTreemap::addNonlinearLeaf) and after
      the application changed the Gaussian in a node (already done by
      \c changeGaussian.

      Sets the node in the optimization queue if according to \c
      featurePassed the set of features passed changed.
   */
  void afterChange ();  


  //! Computes the list of features \c il involved at the Gaussian here
  /*! This is \c gaussian.features but for an inner node it is computed
      from both childrens \c .featuresPassed, so it is valid even if
      the Gaussians are not valid.
  */
  void computeFeaturesInvolved (TmExtendedFeatureList& list) const;  


  //! Computes the height of the subtree at this node
  /*! Since the height is not maintained, the routine recursively
      descends into the tree. If \c max>=0 it descends only until
      \c max thereby limiting computation time.
  */
  int getHeight (int max=-1) const;

  //! Position of \c this relative to \c n
  /*! Returns 0 if this is left descendant of \c n, 1 if right and 2
    if equaln and 3 if both are unrelated.
   */
  int whichSideOf (const TmNode* n) const;  

  //! Computes the number of nodes below (including) this node
  int nrOfNodes() const;  

  //! Adds pointers to all nodes below (including) into a vector
  void nodesBelow (XycVector<TmNode*>& nodes);  
  

  //! Invalidates this node if there is too much linearization error.
  /*! Checks, whether the linearization error in this node concerning
      rotation is too much. The current implementation does nothing.
      A derived class can overload this function. It should check,
      whether the rotation parameter used for this node and the
      current estimate of \c linearizationPointFeature differ too much.
      If they do, the node's Gaussian should be invalidated.
  */
  virtual void checkLinearizationPointFeature ();  

  //! Recursively print indices for \c this and all descendatns
  /*! An example is ((3 4) 5 (7 8)). Returns the lenght of the string. */
  int sprintTreeIndex (char* s);  

  //! Recursively print indices for \c this and all descendatns to stdout
  /*! An example is ((3 4) 5 (7 8)). */
  void printTreeIndex ();  

  //! Makes \c *this child nr. \c childNr of \c parent by swapping children if necessary
  void makeChildNr (int childNr);  

  //! moves \c this and all nodes below to above \c above.
  /*! The parent of \c subtree is deleted and reused as a new
      node above \c above having \c above and \c subtree as its
      children.

      Invalidates \c IS_FEATURE_PASSED_VALID, \c IS_GAUSSIAN_VALID
      and sets all nodes along the path to \c setToBeOptimized.

      \c invalidateFeatureOnly is only for \c bestMoveSubtree. If it
      is true only the \c IS_FEATURE_PASSED_VALID is reset. The \c
      IS_GAUSSIAN_VALID flags are still flagged as valid. This is
      certainly not true and must be corrected later. This is helpful
      if the tree maybe immediately moved back again and the time for
      recomputing the Gaussians may can be saved.

      \c subtree and \c above may neither be identical nor siblings.
  */
  void moveTo (TmNode* above, bool invalidateFeatureOnly=false);  

  
  //! Returns whether this node is a leaf
  bool isLeaf() const
    {return child[0]==NULL;}

  //! Returns whether \c n is an ancestor of \c this
  bool isAncestor(TmNode* n) const;

  

  //! Returns the child index so \c this==parent->child[whichChild()]
  /*! Returns -1 for the root node.
   */
  int whichChild () 
    {
      if (parent==NULL) return -1;
      else if (this==parent->child[0]) return 0;
      else return 1;
    }  

  //! Returns the sibling of this node
  /*! Returns NULL if the node is the root. */
  TmNode* sibling () 
    {
      if (parent==NULL) return NULL;
      if (this==parent->child[0]) return parent->child[1];
      else return parent->child[0];
    }

  //! Returns the sibling and which child this node is
  void getSibling (TmNode*& sibling, int& whichChild)
    {
      if (parent==NULL) {
        sibling = NULL;
        whichChild = -1;
      }
      else if (this==parent->child[0]) {        
        sibling = parent->child[1];
        whichChild = 0;        
      }
      else {
        sibling = parent->child[0];
        whichChild = 1;
      }      
    }  
  

  //! Return how many children pass \c feature
  int degreeOfFeature (TmFeatureId feature) const;


  //! Returns the sum of number of rows of all leaves below \c n
  int rowsBelow () const;


  //! Checks some invariants on this node
  void assertIt () const;  

  //! Asserts, that this node is in the tree \c .tree it says.
  void assertInTree () const;  

  //! Checks some invariants on this node and all ancestors
  void recursiveAssertIt () const;  

  //! Recursively asserts that \c flag is set in all nodes below \c this
  void recursiveAssertFlag (int flag);  


  //! Where is feature \c id represented relative to \c this
  /*! Returns, whether \c id is represented somewhere below \c this 
       (\c below) and somewhere not below \c this (\c notBelow)
       Not efficient. Used in \c assertIt. 
  */
  void isRepresentedWithRespectToNode (TmFeatureId id, bool& below, bool& notBelow) const;  

  //! Internal recursive subroutine for \c isRepresentedWithRespectToN
  void recursiveIsRepresentedWithRespectToNode (const TmNode* n, bool isN2BelowN, TmFeatureId id, bool& below, bool& notBelow) const;

  //! Internal recursive subroutine for \c identifyFeatures
  /*! Changes all occurrences of \c from to \c to in Gaussians at leaves
      below \c this. Invalidates those leaves.
  */
  void recursivelyIdentifyFeature (int from, int to);  

  //! Returns the least common ancestor of \c a and \c b
  static TmNode* leastCommonAncestor (TmNode* a, TmNode* b);  

  //! Returns the largest worstcase cost for \c child[whichChild] leading to \c <=bound at \c this
  /*! If the other child has a too large update cost, \c bound may be exceeded regardless of
      the update cost at child. Then -1 is returned.
  */
  double boundForChild (int whichChild, double bound) {
    bound -= updateCost;
    if (child[1-whichChild]->worstCaseUpdateCost>=bound) return -1;
    else return bound;
  }  

  //! Recursive subfunction for \c findLeavesInvolving
  /*! Adds all leaves below \c this that involve \c id to \c node. 
      \c ctr is incremented by the counter in the leaves found.
   */
  void recursiveAddLeavesInvolving (TmFeatureId id, XycVector<TmNode*>& node, int& ctr);  

  //! Cost of updating the Gaussian of an \c n feature node
  /*! See \c TmTreemap::calibrateGaussianPerformance(). */
  static double updateGaussianCost (int n) 
    {return 1.543037E-6 + n*(1.154801E-6 + n* (44.716E-9 + n*1.799E-9));};

  //! Cost of updating the featurePassed list of an \c n feature node
  static double updateFeaturePassedCost (int n) 
    {n=n; return 3.8E-6;};
  
  
  //! Epsilon used for == checks on costs
  /*! A solution is only accepted as better if it is at least better by
    \c costEps(). */
  static double costEps () {return 1E-9;}

  //! Returns the storage space (Bytes) of this node without children
  virtual int memory() const;  

  //! Returns the storage space (Bytes) of this node with children
  int recursiveMemory() const;
  
  
 private:

  //! Auxiliary routine for \c updateGaussian
  /*! Adds all features from \c fl2 that are not contained in \c featurePassed
      to \c fl. If a feature is already contained it is not added twice. The
      counters are set to 0.
      This routine is used by \c updateGaussian to generate the list of features
      in a Gaussian. The counter must be initialized as 0 because they will be
      later increased by \c TmGaussian::multiply.
  */
  void addMarginalizedFeatures (TmExtendedFeatureList& fl, TmExtendedFeatureList& fl2);

  //! Auxiliary routine for \c updateFeaturePassed
  /*! Merges \c child[0]->featurePassed and \c child[1]->featurePassed removing all features
    that have reached a total count of \c TmFeature::totalCount() and setting TmFeature::marginalizationNode
    of those. The return value specifies the size of the union (i.e. number of features in \c child[0] or \c child[1])
    it counts also those which reached totalCount.    
  */
  int mergeFeatureLists ();  
  


  //! Subroutine for \c printTreeIndex
  void recursivePrintTreeIndex (); 
};

#endif
