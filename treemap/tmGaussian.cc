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
/*!\file tmGaussian.cc 
   \brief Implementation of \c TmGaussian
   \author Udo Frese

  Contains the implementation for class \c TmGaussian representing a
  multidimensional Gaussian. 
*/
#include "tmGaussian.h"
#include "tmExtendedFeatureId.h"
#include <xymlapack/xymLAPACK.h>

TmGaussian::TmGaussian()
  :isTriangular (false), R(), RCompressed(), feature(), linearizationPointFeature(-1), linearizationPoint(0)
{}


TmGaussian::TmGaussian (const TmExtendedFeatureList& feature, int rowsReserved)
  :isTriangular(false), R(), RCompressed(), feature(feature), linearizationPointFeature(-1), linearizationPoint(0)
{
  R.reserve (rowsReserved, feature.size()+1);
  R.create  (0, feature.size()+1);  
}


TmGaussian::TmGaussian (const char* features)
  :isTriangular (false), R(), RCompressed(), feature(), linearizationPointFeature(-1), linearizationPoint(0)
{
  TmExtendedFeatureList fl;  
  if (stringToExtendedFeatureList (fl, features, false)) create (fl, 0);
  else throw runtime_error ("Could not parse feature string.");
}


void TmGaussian::create (const TmExtendedFeatureList& feature, int rowsReserved)
{
  this->feature = feature;
  isTriangular = false;  
  R.reserve (rowsReserved, feature.size()+1);
  R.create  (0, feature.size()+1);
  RCompressed.clear();  
  linearizationPointFeature = -1;
  linearizationPoint = 0;  
}


void TmGaussian::create (const XymMatrixVC& R, const TmExtendedFeatureList& feature, bool isTriangular)
{
  this->feature = feature;
  this->isTriangular = isTriangular;
  this->R = R;
  RCompressed.clear();  
  linearizationPointFeature = -1;
  linearizationPoint = 0;  
}


TmGaussian::TmGaussian (const XymMatrixVC& R, const TmExtendedFeatureList& feature, bool isTriangular)
  :isTriangular(isTriangular), R(R), feature(feature), linearizationPointFeature(-1), linearizationPoint(0)

{
  assert (R.cols()==(int) feature.size()+1);
}


void TmGaussian::triangularize ()
{
  xymGEQR2 (R, R);
  isTriangular = true;  
}

void TmGaussian::triangularize (XymVector& workspace)
{
  xymGEQR2 (R, R, workspace);
  isTriangular = true;  
}


void TmGaussian::multiply (const TmGaussian& gaussian, int fromFeature)
{
  assert (R.isValid());  
  int n = gaussian.rows();
  if (n<0) return;
  // This is  a  rank deficient case that usually does not appear
  // because it means that the Gaussian is even indefinite after being
  // conditioned on all features passed to the parent. This in turn means
  // that the overall solution is indefinite.
  int base = R.rows();  
  R.appendRow (n, true);

  for (int j=fromFeature; j<gaussian.cols(); j++) {
    // search which to which column to copy gaussian.R.col(j)
    int dstJ=-1;
    if (j<(int) gaussian.feature.size()) {
      int srcFeature = gaussian.feature[j].id;      
      for (int i=0; i<(int) feature.size(); i++) 
        if (srcFeature==feature[i].id) {
          dstJ = i;
          break;
        }
      assert (dstJ>=0);
      feature[dstJ].count += gaussian.feature[j].count;    
    }
    else dstJ = R.cols()-1;
    

    // now add gaussian.R.col(j)[fromFeature..] to R.col(dstCol)
    if (!gaussian.RCompressed.empty()) {
      const float *srcP;
      int srcRowInc;      
      double *dstP, *dstPE;
      srcP  = &gaussian.RCompressedAt (fromFeature, j);
      srcRowInc = n-fromFeature-1;      
      R.loopCol (dstJ, dstP, dstPE);
      dstP += base;
      dstPE = dstP + j-fromFeature+1;
      while (dstP!=dstPE) {
        *dstP += *srcP;
        srcP  -= srcRowInc;
        srcRowInc--;        
        dstP++;
      }
    }
    else {      
      double *srcP, *srcE, *dstP, *dstPE;
      gaussian.R.loopCol (j, srcP, srcE);
      srcP += fromFeature;    
      R.loopCol (dstJ, dstP, dstPE);
      dstP += base;
      while (srcP!=srcE) {
        *dstP += *srcP;
        srcP++;
        dstP++;
      }    
    }    
  }
}

void TmGaussian::meanCompressed (float* x, int upToFeature)
{
  assert (isTriangular && !RCompressed.empty());
  if (upToFeature==0) return;  
  int n = feature.size()+1;  
  float* xP = x;
  float* rP = &RCompressedAt (upToFeature-1, n-1); // first entry of \c RCompressed used
  float* xDest  = x+n-upToFeature; // First x entry we will compute (feature[upToFeature-1])
  float* xDestE = x+n; // One after the last x entry we will compute (feature[0])

#ifdef USESSE
#error asdlkasjd
  for (float* xx=xDest; xx!=xDestE; xx++) *xx = 0;
  xDestE[0] = xDestE[1] = xDestE[2] = xDestE[3] = 0; // We have to clear 3 floats after  
  // We clear everything uncomputed so we don't have to care about adding only parts of an xmms registers  
  // eax = %0 = xDestE, ebx = %1 = xP, ecx = %2 = rp, edx = %3 = xDest
  asm volatile (
                "\n\t"
    
                "movl %0, %%eax                    \n\t"
                "movl %2, %%ecx                    \n\t"
                "movl %3, %%edx                    \n\t"
                "cmp %%eax, %%edx                  \n\t"
                "jae .end                          \n\t"
".align 16  \n\t"
".loopCol:  \n\t"
                "xorps %%xmm1,%%xmm1               \n\t"
                "movl %1, %%ebx                    \n\t"
".loopRow:  \n\t"
                "movups (%%ebx), %%xmm2            \n\t"   // fetch *xP
                "movups (%%ecx), %%xmm3            \n\t"   // fetch *rP
                "prefetcht0 32(%%ecx)                \n\t"
                "add $16, %%ebx                    \n\t"
                "add $16, %%ecx                    \n\t"      
                "mulps %%xmm3, %%xmm2              \n\t"   // multiply 
                "subps %%xmm2, %%xmm1              \n\t"   // accumulate
                "cmp %%edx, %%ebx                  \n\t"
                "jb .loopRow                       \n\t"
                "sub %%edx, %%ebx                  \n\t"   // xP has exceeded xDest
                "sub %%ebx, %%ecx                  \n\t"   // so correct rP as if we had exactly hit xDest
                

                "movaps %%xmm1, %%xmm2             \n\t"   // now sum up all 4 floats in xmm1
                "shufps $0x0E, %%xmm2, %%xmm2      \n\t"   // float 2-->0, float 3-->1
                "addps %%xmm2,  %%xmm1             \n\t"
                "movaps %%xmm1, %%xmm2             \n\t"
                "shufps $0x01, %%xmm2, %%xmm2      \n\t" // (2+3) from 1-->0
                "addps %%xmm2, %%xmm1              \n\t"
                "divss (%%ecx), %%xmm1             \n\t" 
                "movss %%xmm1, (%%edx)             \n\t"
                "add $4, %%ecx                     \n\t"
                "add $4, %%edx                     \n\t"
                "cmp %%eax, %%edx                  \n\t"
                "jb .loopCol                       \n\t"
".end: \n\t"
                : : "m" (xDestE), // %0 = xDestE
                  "m" (xP),       // %1 = xP
                  "m" (rP),       // %2 = rP
                  "m" (xDest)     // %3 = xDest
                : "%eax", "%ebx", "%ecx", "%edx", "memory"//, "%xmm1", "%xmm2", "%xmm3", gcc 4.1 does not know these any more, thats strange
                );  

#if ASSERT_LEVEL>=1
  // Check whether everything is correct
  while (xDest!=xDestE) {
    // Handle i-th line of R. Compute x_i (estimate for feature[i] stored in x[n-i-1]) such that (R*x)_i=0
    float sum = 0;
    xP = x;    
    while (xP!=xDest) {
      sum += *rP * *xP;
      rP++;
      xP++;
    }
    float result = -sum / *rP; 
    assert (fabs(result - *xDest)<1E-3);    
    rP++;
    xDest++;
  }
#endif
#else
  while (xDest!=xDestE) {
    // Handle i-th line of R. Compute x_i (estimate for feature[i] stored in x[n-i-1]) such that (R*x)_i=0
    float sum = 0;
    xP = x;    
    while (xP!=xDest) {
      sum += *rP * *xP;
      rP++;
      xP++;
    }
    *xDest = -sum / *rP; 
    rP++;
    xDest++;
  }
#endif
}


void TmGaussian::mean (XymVectorV& x, int upToFeature)
{
  assert (isTriangular);  
  if (upToFeature<0 || upToFeature>(int) feature.size()) upToFeature = feature.size();
  int n = x.size();
  // The routine can handle rank deficient matrices which routinely occur in the algorithm.
  // However for the mean to be well defined, the features conditioned on feature[upToFeature..]
  // must be of full rank. This means that at least the first \c upToFeature rows must be valid.
  assert (R.rows()>=upToFeature);  
  double* xP, *xPP, *xE, *rP;
  xE = x.base()+n;
  xPP = x.base()+upToFeature;
  for (int i=upToFeature-1; i>=0; i--) {
    // solve row i of Ax=0
    int rIncr;    
    R.loopRow (i, rP, rIncr);
    rP += i*rIncr;
    double divisor = *rP;
    rP += rIncr;
    xP = xPP;
    xPP--;

    double sum=0;    
    while (xP!=xE) {
      sum += *xP * *rP;
      xP++;
      rP+= rIncr;
    }
    sum += *rP; // virtual '1' entry appended to 'x'    
    assert (divisor!=0);    
    x[i] = -sum/divisor; // set 'x[i]' so the result in this row is 0
  }
  assert (xPP==x.base());
}


void TmGaussian::clear()
{
  feature.clear();
  R.clear();
  isTriangular = false;  
}


void TmGaussian::setLinearizationPoint (int linearizationPointFeature, double linearizationPoint)
{
  this->linearizationPointFeature = linearizationPointFeature;
  this->linearizationPoint        = linearizationPoint;  
}


void TmGaussian::assertIt () const
{
  assertFinite ();
  assert (RCompressed.empty() || RCompressed.size()==rCompressedSize(feature.size()+1));  
  if (R.isValid() && !RCompressed.empty()) {
    for (int i=0; i<R.cols(); i++) for (int j=0; j<R.cols(); j++) {      
      if (i<=j) {
        if (i<R.rows()) assert (fabs(R(i,j)-RCompressedAt(i,j))<1E-2);
        else assert (RCompressedAt(i,j)==0);
      }    
      else  {
        if (i<R.rows()) assert (R(i,j)==0);    
      }      
    }    
  }  
}


void TmGaussian::assertFinite () const
{
  assert (R.isFinite());  
}


void TmGaussian::computeMarginal (TmGaussian& result, int fromFeature)
{
  assert (isTriangular && 0<=fromFeature && fromFeature<(int) feature.size());  
  result.clear();
  result.isTriangular = true;
  int m = rows()-fromFeature;
  int n = cols()-fromFeature;
  assert (n>0 && m>=0);  
  result.R.create (m, n, false);
  R.extract (result.R, fromFeature, fromFeature);
  result.linearizationPoint = 0;
  result.linearizationPointFeature = -1;
  result.feature.clear();
  result.feature.reserve (n-1);
  for (int i=fromFeature; i<(int) feature.size(); i++) {
    if (feature[i].id==linearizationPointFeature) {      
      result.linearizationPoint = linearizationPoint;
      result.linearizationPointFeature = linearizationPointFeature;
    }
    result.feature.push_back (feature[i]);    
  }  
}


int TmGaussian::memory () const
{
  int mem = sizeof (TmGaussian);
  mem += R.memoryUsage();  // sizeof(XymMatrixVC) is included in \c sizeof(TmGaussian)
  mem += feature.capacity() * sizeof(TmExtendedFeatureId);
  mem += RCompressed.capacity() * sizeof(float);  
  return mem;  
}


void TmGaussian::transferFrom (TmGaussian& g2)
{
  isTriangular = g2.isTriangular;
  R.transferFrom (g2.R);
  feature.clear();
  feature.swap (g2.feature);
  RCompressed.clear();
  RCompressed.swap (g2.RCompressed);  
  linearizationPointFeature = g2.linearizationPointFeature;
  linearizationPoint = g2.linearizationPoint;    
}


void TmGaussian::compress ()
{
  assert (R.isValid() && isTriangular);
  int m = R.rows(), n = R.cols();
  // We extend the matrix with 0s to a full triangle
  RCompressed.clear();
  int rSize = rCompressedSize(n);  
  RCompressed.resizeCompactlyWithUndefindedData (rSize);  
  float* rc = RCompressed.begin();
  int incr = R.colOfs();
  int nm1Incr = incr*(n-1);  
  for (int i=n-1;i>=0;i--) {
    if (i<m) {
      // Reverse copy R(i,j) for (int j=n-1;j>=i;j--)
      double* pEnd = R.rowBase (i);
      double* p    = pEnd+nm1Incr;
      pEnd += incr*(i-1);
      while (p!=pEnd) {
        *rc = (float) *p;
        rc++;
        p-= incr;
      }
    }
    else {
      // Add n-i 0's
      float* rcEnd = rc + (n-i);
      while (rc!=rcEnd) {
        *rc = 0;
        rc++;
      }
    }        
//    for (int j=n-1;j>=i;j--) {
//      if (i>=m) RCompressed.push_back (0);
//      else RCompressed.push_back (R(i,j));    
  }  
  *rc = 0; rc++;
  *rc = 0; rc++;
  *rc = 0; rc++;
  // This is to avoid the SSE routines accessing NaN memory. This is questionable  
#if ASSERT_LEVEL>=2
  assertIt ();  
#endif
  R.clear();  
}


