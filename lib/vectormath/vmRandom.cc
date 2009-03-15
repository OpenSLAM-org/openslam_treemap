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

/*! \file vmRandom.cc
    \brief Implementation of \c VmRandom
    \author Udo Frese

    Implementation of the class \c VmRandom

*/

#include "vmRandom.h"
#include <stdio.h>
#include <stdlib.h>

VmRandom::VmRandom (int seed)
        : _iset(false), _gset(0), _idums(0)
{
    create (seed);
}


double VmRandom::gauss (double mean, double cov)
{
  return gauss()*sqrt(cov) + mean;
}


void VmRandom::gauss (VmVector3& x, const VmMatrix3x3& cov)
{
    VmMatrix3x3 L;
    vmCholesky (cov, L, 1E-9*vmTrace(cov));
    VmVector3 xRaw;
    xRaw[0] = gauss();
    xRaw[1] = gauss();
    xRaw[2] = gauss();
    vmMultiply (x, L, xRaw);
}


void VmRandom::gauss (VmVector2& x, const VmMatrix2x2& cov)
{
    VmMatrix2x2 L;
    vmCholesky (cov, L, 1E-9*vmTrace(cov));
    VmVector2 xRaw;
    xRaw[0] = gauss();
    xRaw[1] = gauss();
    vmMultiply (x, L, xRaw);
}


double VmRandom::gauss()
{
    return gasdev (&_ctr);
}

double VmRandom::random()
{
  return ran4(_idums, &_ctr);
}

double VmRandom::deterministicRandom (long int seed, long int i)
{
  return ran4(seed, &i);
}


double VmRandom::random (double low, double high)
{
    return low + random()*(high-low);
}


void VmRandom::create (long int seed)
{
    _seed = seed;
    _idums = seed;
    _ctr = 0;
    _iset = false;
    _gset = 0;
}


float VmRandom::gasdev(long int* idum)
{
    float fac,r,v1,v2;
    float ran1();
    
    if  (!_iset) {
        do {
          v1=2.0*ran4(_idums, idum)-1.0;
          v2=2.0*ran4(_idums, idum)-1.0;
          r=v1*v1+v2*v2;
        } while (r >= 1.0 || r==0);
        fac=sqrt(-2.0*log(r)/r);
        // we have produced two gaussians, so store one!
        _gset = v1*fac;
        _iset = true;
        return v2*fac;
    } else {
        // We still have a stored gaussian
        _iset = false;
        return _gset;
    }
}


float VmRandom::ran4 (long seed, long *idum)
{
    unsigned long irword, lword;
    // We have 24bit mantissa so we mask out all higher order bi
    // With this setting we exactly return the same numbers as the
    // original numerical recipes code but in a portable way. The
    // original code was namely wrongly optimized by gcc 4.0.
    static unsigned long jflmsk = 0x007fffff;
    
    irword=(*idum);
    lword=seed;
    psdes(&lword,&irword);
    ++(*idum);
    irword &= jflmsk;
    return float(irword)/float(jflmsk+1);    
}

//! Number of iterations of the pseudo DES algorithm used in the random number generator
#define NITER 4

void VmRandom::psdes(unsigned long *lword, unsigned long *irword)
{
    unsigned long i,ia,ib,iswap,itmph=0,itmpl=0;
    static unsigned long c1[NITER]={
        0xbaa96887L, 0x1e17d32cL, 0x03bcdc3cL, 0x0f33d1b2L};
    static unsigned long c2[NITER]={
        0x4b0f3b58L, 0xe874f0c3L, 0x6955c5a6L, 0x55a7ca46L};
    
    for (i=0;i<NITER;i++) {
        ia=(iswap=(*irword)) ^ c1[i];
        itmpl = ia & 0xffff;
        itmph = ia >> 16;
        ib=itmpl*itmpl+ ~(itmph*itmph);
        *irword=(*lword) ^ (((ia = (ib >> 16) |
                              ((ib & 0xffff) << 16)) ^ c2[i])+itmpl*itmph);
        *lword=iswap;
    }
}

#undef NITER


//! Internal class used to call a test for ran4 on program startup
class Checker 
{
public:
  Checker ();  
};

//! Auxiliary variable used to call \c Checker::Checker on library startup
Checker checker;

// This routine is called on program startup and verifies that ran4 works correctly
Checker::Checker ()
{
  long idum =1, seed = 1;
  float rv1 = VmRandom::ran4(seed, &idum);

  idum = 99;
  float rv2 = VmRandom::ran4(seed, &idum);

  seed = 99;  
  idum = 1;
  float rv3 = VmRandom::ran4(seed, &idum);

  idum = 99;
  float rv4 = VmRandom::ran4(seed, &idum);
  if (rv1!=0.2191203833f ||
      rv2!=0.8492462635f ||
      rv3!=0.3752903938f ||
      rv4!=0.4573339224f) {
    fprintf (stderr, "Random number generator VmRandom::ran4 is not working properly\n.");
    abort();    
  }
}
