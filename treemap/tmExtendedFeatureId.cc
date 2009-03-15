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
//!\author Udo Frese

/*!\file tmExtendedFeatureId.cc Contains the implementation for class
  \c TmExtendedFeatureId
*/

#include <algorithm>
#include "tmExtendedFeatureId.h"

bool isElement (const TmExtendedFeatureList& list, TmFeatureId id)
{
  int lo=0, hi=(int) list.size()-1;
  if (hi<0) return false; // lmList empty
  if (id<list[lo].id) return false;
  if (list[hi].id<id) return false;    
  if (list[hi].id==id) return true;
  // Invariant: lm is contained in [lo..hi]
  while (lo<hi) {
    int mid = (lo+hi)/2; // it is important that 1/2==0
    TmFeatureId lmMid = list[mid].id;
    if (id<lmMid) hi = mid;
    else lo = mid+1;
  }
  return (lo>0 && list[lo-1].id==id);
}


bool isElement (const TmFeatureList& list, TmFeatureId id)
{
  for (int i=0; i<(int) list.size(); i++) if (list[i]==id) return true;
  return false;  
}


void sumUp (TmExtendedFeatureList& list)
{
  sort (list.begin(), list.end());
  int j=-1;  // j is the large entry of the already merged part of list
  for (int i=0; i<(int) list.size(); i++) {
    if (j>=0 && list[i].id==list[j].id) 
      list[j].count += list[i].count; // merge into list[j]
    else {
      j++; // add new entry to merged part
      list[j] = list[i];
    }    
  }
  list.resize (j+1);  
}


void merge (TmExtendedFeatureList& result, const TmExtendedFeatureList& a, const TmExtendedFeatureList& b)
{
    TmExtendedFeatureList tmp;
    int szA = a.size();
    int szB = b.size();
    tmp.reserve (szA+szB);
    int j=0;
    for (int i=0; i<szA; i++) {
        while (j<szB && b[j].id<a[i].id) {
            tmp.push_back(b[j]); // Only in 'b'
            j++;
        }
        if (j<szB && a[i].id==b[j].id) {
          tmp.push_back (TmExtendedFeatureId (a[i].id, a[i].count+b[j].count));
          j++;  // is contained in both
        }
        else tmp.push_back(a[i]); // only In 'a'
    }
    while (j<szB) {
        tmp.push_back(b[j]); // Only in b
        j++;
    }
    for (int i=0; i<(int) tmp.size(); i++) assert(tmp[i].id>=0);    
    result = tmp;
}



int featureToString (char* txt, TmFeatureId lm)
{
    if (lm<0) {
      txt[0]='%';      
      txt[1]='\0';
      return 1;
    }
    // Landmarks are shown as 'a', 'b', .., 'z', 'aa', .., 'az', 'ba'..'bz', 'za'..'zz'
    int len;
    if (lm>=26) { // more than one digit
        len = featureToString (txt, (lm-26)/26);
    }
    else len=0;
    txt[len  ]='a'+(lm%26);
    txt[len+1]='\0';
    return len+1;
}


int stringToFeature (const char* txt)
{
    if (txt[0]<'a' || txt[0]>'z') return -1;
    int lm=0;
    int i=0;
    while (txt[i]!='\0') {
        if (txt[i]<'a' || txt[i]>'z') return -1;
        if (i>0) lm = (lm+1)*26;
        lm = lm + (txt[i]-'a');
        i++;
    }
    return lm;
}


bool stringToExtendedFeatureList (TmExtendedFeatureList& fl, const char* txt, bool doSort)
{
    int j=0;
    fl.clear();
    while (txt[j]==' ') j++;
    while (txt[j]!='\0') {
        // Extract the text for one landmark into buffer ( to next ' ')
        char buffer[8];
        int k=0;
        while (k<7 && txt[j]!=' ' && txt[j]!='\0') {
            buffer[k] = txt[j];
            k++;
            j++;
        }
        buffer[k]='\0';
        if (k>=7) { // Error
            fl.clear();
            return false;
        }
        // Convert to a landmark and add to 'fl'
        int oneF = stringToFeature(buffer);
        if (oneF==-1) {
            fl.clear();
            return false;
        }
        else fl.push_back (TmExtendedFeatureId (oneF, 1));
        while (txt[j]==' ') j++;
    }
    if (doSort) listSortAndMakeUnique (fl, fl);
    return true;
}


void listSortAndMakeUnique (TmExtendedFeatureList& result, const TmExtendedFeatureList& fl)
{
  result = fl;  
  sort(result.begin(), result.end());
  // Sum up duplicate entries
  int j=0;
  for (int i=0; i<(int) result.size(); i++) {
    if (j>0 && result[i].id==result[j-1].id) result[j-1].count += result[i].count;
    else {
      result[j] = result[i];      
      j++;
    }
  }
  result.resize (j);  
}


void pfl (TmFeatureList& fl)
{
  printf("{");  
  for (int i=0; i<(int) fl.size(); i++) {
    if (i>0) printf(" ");
    char txt[10];
    featureToString (txt, fl[i]);
    printf ("%s", txt);
  }
  printf("}\n");  
}

void pfl (TmExtendedFeatureList& fl)
{
  printf("{");  
  for (int i=0; i<(int) fl.size(); i++) {
    if (i>0) printf(" ");
    char txt[10];
    featureToString (txt, fl[i].id);
    printf ("%s(%d)", txt, fl[i].count);
  }
  printf("}\n");  
}


bool isIntersectionEmpty (const TmExtendedFeatureList& a, const TmFeatureList& b)
{
  int szA = a.size();
  int szB = b.size();
  int j=0;
  for (int i=0; i<szA; i++) {
    while (j<szB && b[j]<a[i].id) j++;  // Only in 'b'
    if (j<szB && a[i].id==b[j]) return false; // is contained in both
    // else only in 'a'
  }
  // remainings 'j's are only in b
  return true;
}

bool isIntersectionEmpty (const TmExtendedFeatureList& a, const TmExtendedFeatureList& b)
{
  int szA = a.size();
  int szB = b.size();
  int j=0;
  for (int i=0; i<szA; i++) {
    while (j<szB && b[j].id<a[i].id) j++;  // Only in 'b'
    if (j<szB && a[i].id==b[j].id) return false; // is contained in both
    // else only in 'a'
  }
  // remainings 'j's are only in b
  return true;
}


int nrOfIntersecting (const TmExtendedFeatureList& a, const TmExtendedFeatureList& b)
{
  int szA = a.size();
  int szB = b.size();
  int j=0;
  int ctr=0;  
  for (int i=0; i<szA; i++) {
    for (j=0; j<i; j++) if (a[j].id==a[i].id) break;
    if (j==i) {      
      for (j=0; j<szB; j++) if (a[i].id==b[j].id) break;
      if (j<szB) ctr++;
    }    
  }
  return ctr;  
}
