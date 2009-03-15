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

/*!\file xymVectorV.h Contains the class \c XymVectorV representing a
  vector view.
*/


#ifndef XYMVECTORV_H
#define XYMVECTORV_H

#include "xymGeneral.h"
#include <vector>
class XymVector;

//! Vector view
/*! A Vector view, that can cope with non unit stride vectors,
    like the row of a column major, the column of a row major
     or the diagonal of a matrix. Does not own the memory.
*/
class XymVectorV
{
public:
    //! Empty constructor: A non valid vector.
    XymVectorV ();

    //! Constructor from raw data
    /*! Constructs a view to the vector following \c *d
        with \c n entries having an offst of \c stride.
        So the \c i-th entry is \c d[i*stride].
    */
    XymVectorV (double* d, int n, int stride=1);

    //! Element access
    /*! Range checking is performed, if the code is neither
        compiled with \c NDEBUG nor with \c NORANGECHECKING.
    */
    double& operator() (int i) {
        XYMRANGECHECK ((0<=i) && (i<_n));
        return _d[i*_stride];
    }


    //! Element access
    /*! Range checking is performed, if the code is neither
        compiled with \c NDEBUG nor with \c NORANGECHECKING.
    */
    double operator() (int i) const {
        XYMRANGECHECK ((0<=i) && (i<_n));
        return _d[i*_stride];
    }

    //! Element access
    /*! Range checking is performed, if the code is neither
        compiled with \c NDEBUG nor with \c NORANGECHECKING.
    */
    double& operator[] (int i) {
        XYMRANGECHECK ((0<=i) && (i<_n));
        return _d[i*_stride];
    }


    //! Element access
    /*! Range checking is performed, if the code is neither
        compiled with \c NDEBUG nor with \c NORANGECHECKING.
    */
    double operator[] (int i) const {
        XYMRANGECHECK ((0<=i) && (i<_n));
        return _d[i*_stride];
    }

    //! Element access
    /*! Range checking is performed, if the code is neither
        compiled with \c NDEBUG nor with \c NORANGECHECKING.
    */
    double& get (int i) {
        XYMRANGECHECK ((0<=i) && (i<_n));
        return _d[i*_stride];
    }

    //! Element access
    /*! Range checking is performed, if the code is neither
        compiled with \c NDEBUG nor with \c NORANGECHECKING.
    */
    double get (int i) const {
        XYMRANGECHECK ((0<=i) && (i<_n));
        return _d[i*_stride];
    }

    //! Element access without range checking
    double& getNC (int i) {
        return _d[i*_stride];
    }

    //! Element access without range checking
    double getNC (int i) const {
        return _d[i*_stride];
    }

    //! Number of elements
    int size() const {
            return _n;
        }

    //! Whether the vector is initialized (has a certain size, even 0)
    bool isValid()  const {
        return _d!=NULL;
    }

    //! Whether the vector is empty (has size 0)
    bool isEmpty() const 
        {return _n==0;}
    
    
    //! Returns a view to the subvector starting at entry \c i
    //! with size \c n entries.
    XymVectorV subVector (int i, int n) {
        XYMRANGECHECK (0<=i && i+n<=_n);
        return XymVectorV (_d+i*_stride, n, _stride);
    }

    //! Returns a view to the subvector starting at entry \c i
    //! with size \c n entries.
    const XymVectorV subVector (int i, int n) const {
        XYMRANGECHECK (0<=i && i+n<=_n);
        return XymVectorV (_d+i*_stride, n, _stride);
    }


    //! Twodimensional vector
    typedef double Vector2[2];

    //! Threedimensional vector
    typedef double Vector3[3];
    
    //! Extracts a 2 vector starting at position \c i
    void extract (Vector2& v, int i) const {
        XYMRANGECHECK (0<=i && i+2<=_n);
        v[0] = _d[i];
        v[1] = _d[i+_stride];
    }
    

    //! Extracts a 3 vector starting at position \c i
    void extract (Vector3& v, int i) const {
        XYMRANGECHECK (0<=i && i+3<=_n);
        v[0] = _d[i];
        v[1] = _d[i+_stride];
        v[2] = _d[i+_stride+_stride];
    }
    
    //! Stores a 2 vector starting at position \c i
    void store (const Vector2& v, int i) {
        XYMRANGECHECK (0<=i && i+2<=_n);
        _d[i] = v[0] ;
        _d[i+_stride] = v[1];
    }
    
    //! Stores a 3 vector starting at position \c i
    void store (const Vector3& v, int i) {
        XYMRANGECHECK (0<=i && i+3<=_n);
        _d[i] = v[0];
        _d[i+_stride] = v[1];
        _d[i+_stride+_stride] = v[2];
    }

    //! Asserts, that \c this is valid and a \c n vector
    void assertFormat (int n) const {
        XYMFORMAT (_d!=NULL && _n==n);
    }
    
    //! Returns, whether \c this contains only finite numbers
    bool isFinite() const;

    //! Sets \c *this(i)=v(idx[i]).
    void extractFrom (const XymVectorV& v, const vector<int>& idx);

    //! Copies the data from \c v into \* this
    void copyFrom (const XymVectorV& v);

    //! Fills the vector with \c alpha
    void fill (double alpha);
    

    //! Print the vector
    void print() const;

    //! See the operator itself.
    friend ostream& operator<<(ostream& os, const XymVectorV& v);

    //! See the operator itself.
    friend istream& operator>>(ifstream& ifs, XymVectorV& v) throw (runtime_error);


    //************************************************************
    /*!\name Direct Access Methods

       It is often much faster to use direct access to the data via
       'double*' pointer.  As this vector implements a concretely
       defined memory model (sequentially with non-unit stride) it is
       allowed to directly access the elements via pointers returned
       by the following functions.  However we \em strongly \em advice
       to limit this type of access to the inner loop of computational
       important algorithms, to keep encapsulation as strong as
       possible.

       The location of the \c i-th element is \c base()+i*elmOfs().

       The parameter to pass to LAPACK are: \c base(), \c size(), \c
       elmOfs(). Most LAPACK functions can't take non-unit stride
       vectors, in which case a copy must be made, if \c elmOfs()!=1.

       It is most advisable to use Level 1 BLAS function contained in
       the file \c xyBLAS.h.
     */
     //@{

    //! Pointer to element 0. 
    double* base() const {
            return _d;
        }

    //! How much to add to the pointer to get to the next element.
    int elmOfs() const {
            return _stride;
        }

    //! Special function for BLAS/LAPACK
    /*! Returns a pointer to the integer member variable holding the
        stride of the vector. This pointer can directly be passed to
        FORTRAN.
    */
    int* st() const {
            return (int*) &_stride;
        }
    

    //! Function to construct efficient loops
    /*!
    Most linear algebra operations perform some elementwise
    operation on some vectors or rows/column of a matrix in
    the inner loop. To speed up that loop, one can implement
    it using pointers that are initialized by 'loop(...)'
    like in the following sheme:

    \code
    double *p, *pEnd;
    int inc;
    v.loop (p, pEnd, inc);
    for (;p!=pEnd;p+=incr) {
        // operations on a single vector element
    }
    \endcode
    */
    void loop(double*& pStart, double*& pEnd, int& increment) const
        {
            pStart    = _d;
            increment = _stride;
            pEnd      = pStart+_n*_stride;
        }
    
    //! Function to construct loops (See \c loop)
    void loop(double*& pStart, int& increment) const
        {
            pStart    = _d;
            increment = _stride;
        }
    
//@}

protected:
    //! Pointer to the first element \c (0)
    double* _d;

    //! Number of elements
    int _n;

    //! stride (How many increment to \c _d for one element)
    int _stride;
};

//! Prints \c v
/*! This function can be used in the debugger to view the
    content of a vector. */
void print (const XymVectorV& v);


//! Checks, whether 'v1' and 'v2' have the same length
inline void checkEqualLength (const XymVectorV& v1, const XymVectorV& v2)
{
    XYMASSERT (v1.size()==v2.size(), "Vectors must have same size");
}

//! stores the vector \c v into stream \c os
ostream& operator<<(ostream&, const XymVectorV& v);


//! loads the vector from the stream.
/*! \warning A \c XymVectorV can not change it's size, thus the
  vector specified in the stream must have the correct number of
  elements.  If the vectors size is a-priori unknown, use \c
  XymVector::operator>>
*/
istream& operator>>(ifstream&, XymVectorV& v) throw (runtime_error);

//! Computes a+b
/*! Note, c=a+b is less efficient than add \c (c, a, b) (\c xymOperations.h). */
XymVector const operator+ (const XymVectorV& a, const XymVectorV& b);

//! Computes a=a+b
void operator+= (XymVectorV& a, const XymVectorV& b);

//! Computes a-b
/*! Note, c=a-b is less efficient than sub \c (c, a, b) (\c xymOperations.h). */
XymVector const operator- (const XymVectorV& a, const XymVectorV& b);

//! Computes a=a-b
void operator-= (XymVectorV& a, const XymVectorV& b);

//! Computes a*labmda
/*! Note, c=a*lambda is less efficient than multiply \c (c, a, lambda) (\c xymOperations.h). */
XymVector const operator* (const XymVectorV& a, double lambda);

//! Computes lambda*a=a*lambda
/*! Note, c=a*lambda is less efficient than multiply \c (c, a, lambda) (\c xymOperations.h). */
XymVector const operator* (double lambda, const XymVectorV& a);

//! Computes a=a*lambda
void operator*= (XymVectorV& a, double lambda);

//! Dot product a*b
double operator* (const XymVectorV& a, const XymVectorV& b);

#endif
