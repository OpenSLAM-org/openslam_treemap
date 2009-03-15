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

#ifndef TMFEATURE_H
#define TMFEATURE_H


/*!\file tmFeature.h
   \brief Class \c TmFeature representing information about a single 1-D feature.

   \author Udo Frese
*/

#include "tmTypes.h"

class TmNode;

//! Information stored for a single feature in the treemap (globally)
/*! Features are 1-D i.e. a landmarks x or y coordinate or the robot's
  x, y, OR theta coordinate.
*/
class TmFeature 
{
 public:
  //! empty feature with null link
  TmFeature ()
    :est (tmNan), marginalizationNode (NULL), multiPurposeField(0)
    {}
  

  //! The estimate (output of the algorithm)
  float est;
  
  //! The node at which this feature is marginalized out.
  /*! Distributions is passed upwards in the tree. If at some node
    all distributions involving a feature are integrated,
    i.e. are descendants of that node, all information is there
    and the feature can be marginalized out at that node. The
    decision is based on the features counter at the node and \c
    totalCount. This node is the marginalization node. The
    algorithm may decide to choose an ancestor node, i.e. to
    marginalize higher in the tree. This option is both used to
    simplify bookkeeping and to prevent marginalization of a
    feature that is used as linearization point of a
    distribution.
    
    The treemap algorithm sets \c marginalizationNode when
    marginalizing out, precisely when updating the \c .featurePassed
    lists. So when the treemap is not fully updated, not all
    marginalization nodes are up to date. If it is not up to date, \c
    marginalizationNode is set to \c NULL. In this case the former
    marginalization node has already been invalidated and will be
    updated later. If \c treemap.isFeaturePassedValid() then also
    the marginalization nodes are valid.
  */
  TmNode* marginalizationNode;
  
  
  
  //! The sum of all counters for that feature in leafs
  /*! Every original distribution passed as input to the treemap
    and stored in a leaf includes a counter (set to 1) for
    each feature it involves. The sum of all these counters is
    stored in \c totalCount. They are added when distributions
    are integrated along the tree. If the counter in a node 
    equals \c totalCount then the feature is not involved in
    distributions outside the node and can be marginalized out
    at that node.
    
    If \c totalCount is 0<=0 that feature is not represented in
    the treemap. This can happen for features which the
    algorithm is allowed to permanently marginalize out
    automatically (e.g. robot poses). In this case the member is
    used to provide a linked list of unused features: \c
    totalCount is minus the index of the next unused feature in
    the list minus 1 or it is 0 at the end of the list.
  */
  int totalCount() const
    {
      if (multiPurposeField>=0) return multiPurposeField & 0xffffff;
      else return 0;
    }

  //! Returns totalCount() for a feature that is existing
  /*! This avoids the 'if' test in \c totalCount() and is (moderately)
      important for \c TmNode::mergeFeatureLists
   */
  int totalCountOfExistingFeature () 
  {
    return multiPurposeField & 0xffffff;
  }  
  
  //! Sets \c totalCount in \c multiPurposeField
  /*! The feature must be active. I.e. not an unused feature.
   */
  void setTotalCount (int totalCount)
    {
      assert (multiPurposeField>=0);
      multiPurposeField = (multiPurposeField & 0xff000000) + totalCount;
    }

  //! Add count to \c totalCount
  void addTotalCount (int count)
    {
      assert (multiPurposeField>=0);
      multiPurposeField += count;      
    }
  
  //! Next unused feature entry
  /*! The unused feature entries are linked by a list stored in \c
    multiPurposeField. Returns the index of the next unused
    feature or -1 if there is no next. This function is only
    defined, if it is unused.  (\c isFlag (IS_EMPTY)).

    We maintain linked lists for unused feature blocks of different
    size. In each block only the first feature is contained in the
    linked list. All other have -1 links.
  */
  int nextUnusedFeature () const
    {
      assert (multiPurposeField<0);
      return -2-multiPurposeField;
    }
  
  //! Sets the link to the next unused feature entry
  void setNextUnusedFeature (int feature)
    {
      assert (multiPurposeField<0);
      multiPurposeField = -2-feature;
    }
  
  
  //! Different flags concerning a feature
  /*
    The application can define its own flags in the bits set in
    \c USER_FLAGS. This is used by different drivers (\c tmSlamDriver2DL)
    to indicate whether a feature is an x, y or theta coordinate of
    a landmark or of a robot pose. bit 0-23 are used by \c totalCount
  */
  enum FeatureFlags
    {
      IS_EMPTY                = 0x80000000,
      CAN_BE_MARGINALIZED_OUT = 0x40000000,
      CAN_BE_SPARSIFIED       = 0x20000000,
      USER_FLAGS              = 0x0f000000
    };
  /*! \var IS_EMPTY 

    whether this feature entry is currently
    unused. If this is true, \c multiPurposeField<0 and all
    other flags as well as \c totalCount are undefined.
  */
  /*! \var CAN_BE_MARGINALIZED_OUT 

     whether the algorithm is allowed
    to permanently marginalize out this feature removing it
    completely from the representation and not providing
    estimates for that feature any more. Examples are old robot
    poses. The application can change this flag to \c true at
    any time. Once it is \c true it cannot be unset anymore.
  */
 /*! \var CAN_BE_SPARSIFIED 
     
    If this flag is set, the algorithm can
    marginalize it out of certain distribution even if it is involved
    in another distribution. This is equivalent to treating both
    occurences of the feature as independent features and thus is a
    conservative way of sparsification. Implies \c
    CAN_BE_MARGINALIZED_OUT.
 */




  //! Returns whether this feature is empty, i.e. unused
  bool isEmpty () const
    {return multiPurposeField<0;}
  
  
  //! Returns, whether a flag is set in \c multiPurposeField
  bool isFlag (int flag) const
    {
      return ((multiPurposeField & flag)!=0);
    }

  //! Returns whether all \c mustFlag flags are set and no \c mayNotFlag
  bool flagCondition (int mustFlag, int mayNotFlag) const
    {
      return ((multiPurposeField & mustFlag)==mustFlag) && ((multiPurposeField & mayNotFlag)==0);      
    }
  
  
  //! Sets alls flags in \c flag to \c flagTo
  /*! If the feature is switched to ISEMPTY, totalCount() is set
    to 0, if the feature is switched to ~ISEMPTY,
    nextUnusedFeature() is set to -1. This is the only way to
    switch between empty and not empty mode.
  */
  void setFlag (int flag, int flagTo)
    {
      if (multiPurposeField<0 && (flag & IS_EMPTY)!=0 && (flagTo & IS_EMPTY)==0) {
        // make not empty
        multiPurposeField = flagTo;
        est = tmNan;        
      }
      else if (multiPurposeField>=0 && (flag & IS_EMPTY)!=0 && (flagTo & IS_EMPTY)!=0) {
        // make empty
        multiPurposeField = -1;
        est = tmNan;        
      }
      else multiPurposeField = (multiPurposeField & ~flag) | flagTo; // update flags
    }

  //! Returns whether there should be an estimate for this feture
  /*! A feature is defined if it is not empty and totalCount()>0 */
  bool isDefined () const
    {return !isEmpty() && totalCount()>0;}

  //! Returns the user flag part
  int userFlag () const
    {
      return multiPurposeField & USER_FLAGS;
    }
  
  
  
  
  
 protected:
  //! An int containing \c totalCount(), \c isFlag() and \c nextUnusedFeature()
  /*! Don't use it directly. Use the provided access functions instead.
   */
  int multiPurposeField;

};
#endif
