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
#ifndef TMEXTENEDFEATUREID_H
#define TMEXTENEDFEATUREID_H
/*!\file tmExtendedFeatureId.h 
   \brief Contains \c TmExtendedFeatureId
   \author Udo Frese

   Contains \c TmExtendedFeatureId
*/


#include "tmTypes.h"

//! Used in maintaining the list of features represented at a node
/*! The \c TmExtendedFeatureId objects is stored in a \c
    TmExtendedFeatureList at a node. It contains a feature id together
    with the number of leaves below the node that involve that
    feature.  This information is used to determine when a feature can
    be marginalized out by comparing \c count to
    TmFeature::totalCount().
*/
class TmExtendedFeatureId 
{
 public:
  //! The id (index to \c TmTreemap::feature) of that feature
  int id;
  //! How many original distributions involving that feature have been multiplied into this distribution.
  int count;
  
  //! Default constructor
  TmExtendedFeatureId ()
    : id(-1), count(0)
    {}

  //! Std constructor
  TmExtendedFeatureId (int id, int count=1)
    :id(id), count(count)
    {}    

  //! sort by feature indices
  bool operator < (const TmExtendedFeatureId& b) const
    {return id<b.id;}

  //! Both id and count must be equal
  bool operator == (const TmExtendedFeatureId& b) const
    {
      return id==b.id && (count==b.count);
    }  
};

//! A list of features with corresponding counters
/*! The list is sorted by ascenting \c .id unless noted otherwise.
 */
typedef XycVector<TmExtendedFeatureId> TmExtendedFeatureList;

//! Returns whether \c id is contained in \c list.
bool isElement (const TmExtendedFeatureList& list, TmFeatureId id);

//! Returns whether \c id is contained in \c list.
bool isElement (const TmFeatureList& list, TmFeatureId id);

//! sorts \c list and combines repeated entries by adding \c .count
void sumUp (TmExtendedFeatureList& list);

//! Merges two sorted lists so the result is sorted and duplicate entries are added
/*! \c result may be identical to \c a or \c b
*/
void merge (TmExtendedFeatureList& result, const TmExtendedFeatureList& a, const TmExtendedFeatureList& b);

//! Returns whether the is no feature that is contained in both lists
/*! \c a and \c b must be sorted 
 */
bool isIntersectionEmpty (const TmExtendedFeatureList& a, const TmFeatureList& b);

//! Returns whether the is no feature that is contained in both lists
/*! \c a and \c b must be sorted 
 */
bool isIntersectionEmpty (const TmExtendedFeatureList& a, const TmExtendedFeatureList& b);

//! Returns the number of features both in \c a and \c b
/*! \c a and \c b may be unsorted. */
int nrOfIntersecting (const TmExtendedFeatureList& a, const TmExtendedFeatureList& b);


//! Function to convert a landmark id to a string name consisting of noncapital letters.
/*! The name is stored into \c txt (with 0 termination) and the length
    is returned. The landmark names are:
  
     a,b,c,...,z,aa,ab,...,az,ba,bb,bz,...,za,zb,...,zz,aaa,...
     
     The length is returned and it is never greater than 5.
*/
int featureToString (char* txt, TmFeatureId lm);

//! Converts a string naming a feature number into the feature
/*! See \c featureToString */
int stringToFeature (const char* txt);


//! Converts a string with blank separated feature numbers into a \c TmFeatureList
/*! The return value specifies, whether conversion was successfull. If
    \c doSort is \c true, the list is sorted otherwise features are
    returned as ordered in \c txt.
*/
bool stringToExtendedFeatureList (TmExtendedFeatureList& fl, const char* txt, bool doSort=true);

//! Sorts \c fl and removes duplicates
void listSortAndMakeUnique (TmExtendedFeatureList& result, const TmExtendedFeatureList& fl);

//! Prints a feature list to stdout
/*! For debugging */
void pfl (TmFeatureList& fl);

//! Prints an extended featurelist to stdout
/*! For debugging */
void pfl (TmExtendedFeatureList& fl);

#endif
