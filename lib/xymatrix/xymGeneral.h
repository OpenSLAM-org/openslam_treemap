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

/*!\file xymGeneral.h Contains general definitions (types, macros, auxiliary
   functions) for the XyMatrix package.
*/

#ifndef XYMGENERAL_H
#define XYMGENERAL_H



#include <stdlib.h>
//#include <algo.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdexcept>

using namespace std;


//! xymreal can be either double default or float (with \C XYMFLOAT) defined
#ifndef XYMFLOAT
typedef double xymreal;
#else
typedef float xymreal;
#endif

//! Whether the xymatrix package has been compiled with NDEBUG (and thus hopefully -O3)
bool xymIsCompiledWithOptimization();


#ifndef NDEBUG
//! Extened \c assert macro for the xyMatrix package
/*! If \c B does not hold an error message is printed and the program
    is aborted via an segmentation fault. The advantage of using a
    segmentation fault is, that no subroutine is called, so the
    debugger is automatically at the right place and the program can
    be continued in the debugger using "set execution position".
*/
#define XYMASSERT(B,MSG) if (!(B)) (*xymDummyP=fprintf(stderr, "%s [%s:%d]\n%s\n", MSG, __FILE__, __LINE__, __STRING(B)))
#else

#define XYMASSERT(B,MSG) {}
#endif

//! Macro for format assertion
#define XYMFORMAT(B) XYMASSERT(B,"Format error in matrices")

#if !defined(NORANGECHECKING) && !defined(NDEBUG)

//! Macro for range checking error
/*! Range checking can be switched off using \c NDEBUG or \c
    NORANGECHECKING.
 */
#define XYMRANGECHECK(B) XYMASSERT(B,"index range check error while accessing matrix/vector")

#else

#define XYMRANGECHECK(B) 

#endif

//! Internal pointer variable initialised to \c NULL for generating segmentation faults.
extern int* xymDummyP;

//! Auxiliary function to read numbers from streams.
/*! Reads the next word from \c in and asserts, that it is equal to
    \c word and throws a \c runtime_error otherwise.
    This is an auxiliary routine, that can be used by this class and
    derived classes to implement \c load.
*/
void xymExpect(istream& in, const string& word) throw (runtime_error);


#endif
