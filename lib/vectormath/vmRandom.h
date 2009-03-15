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

/*! \file vmRandom.h
    \brief A high quality random number generator \c VmRandom.
    \author Udo Frese

    A high quality random number generator \c VmRandom.
*/
#ifndef VMRANDOM_H
#define VMRANDOM_H

#include "vectormath.h"

//! This class encapsulates a high quality random number generator
/*! The generator is based on a reimplementation of the numerical
    recipes ran4 and gausdev algorithm. The routine \c VmRandom::ran4
    avoids a binary cast from int to float to avoid problems with 
    gcc 4.0 but returns exactly the same numbers as the method proposed
    in Numerical recipes.

    \sa Numerical Recipes 3rd Edition: The Art of Scientific Computing 
        by William H. Press (Author), Saul A. Teukolsky (Author), 
        William T. Vetterling (Author), Brian P. Flannery (Author) 

        section ???
 */
class VmRandom
{
public:
    //! See \c create
    VmRandom(int seed=0);

    //! Initialise the random number generator with seed \c seed
    void create (long int seed=0);

    //! Returns a uniform [0..1] random deviate
    double random();

    //! Returns the i-th uniform [0..1] random number when starting with \c seed
    /*! The result is reproducible and does not depend on the random state. */
    static double deterministicRandom (long int seed, long int i);    

    //! Returns a uniform [lo..hi] random deviate
    double random (double low, double high);

    //! Returns a gaussian random deviate
    /*! mean 0, std. deviation 1 */
    double gauss();

    //! Returns a gaussian random deviate with given mean and covariance
    double gauss(double mean, double cov);    

    //! Return a 3D gaussian random deviate with covariance \c cov
    void gauss (VmVector3& x, const VmMatrix3x3& cov);

    //! Return a 3D gaussian random deviate with covariance \c cov
    void gauss (VmVector2& x, const VmMatrix2x2& cov);


protected:
    //! Random seed
    long int _seed;
    //! Counter of the random number used
    long int _ctr;

private:
    //! Internal routine for \c gauss()
    float gasdev(long int* idum);

    //! Internal routine for \c random()
    static float ran4(long seed, long *idum);

    //! Internal routine for \c random()
    static void psdes(unsigned long *lword, unsigned long *irword);

    //! Auxiliary variable for \c gasdev and \c ran4
    bool _iset;
    //! Auxiliary variable for \c gasdev and \c ran4
    float _gset;
    //! Auxiliary variable for \c gasdev and \c ran4
    long _idums;

    //! This class is used to check the random number generator on startup.
    friend class Checker;    
};

#endif  /* RANDOM_H */
