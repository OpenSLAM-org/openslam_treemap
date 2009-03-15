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

/*!\file xymMatrixVC.h Contains the class \c XymMatrixVC representing a
  matrix view.
*/

#ifndef XYMMATRIXVC_H
#define XYMMATRIXVC_H

#include "xymVector.h"
class XymMatrixC;

//! \c n times \c m matrix (dense column major) without own memory
/*! It is a view to a matrix, thus not owning the memory.  The model
    used is compatible with LAPACK and allows to define arbitry
    rectangular submatrices of a matrix as a matrix view.  Indices are
    row, column and range from 0 to dimension-1.

   In general an \c n times \c m matrix is a matrix with \c n rows and
   c\ m column and the element \c (i,j) is the element with row \c i
   and column \c j. Indices range from \c 0 to \c n-1 or \c m-1.  In
   the declaration \c i always refers to a row index, \c j to a column
   index, \c n to a number of rows (usually a matrix dimension) and \c
   m to a number of columns.

   The class implements range checking if neither \c NORANGECHECKING
   nor \c NDEBUG are set.

   Routines not tested are marked by a warning. 
*/
class XymMatrixVC 
{
public:
    //! Conctruct an emptty (uninitialised) matrix view.
    XymMatrixVC();

    //! Conctruct matrix view from raw data.
    /*! Returns a Matrix View that acts on the numbers stored in \c *d
        (column major), with dimension \c n, \c m and leading
        dimension (column) increment \c ld. If \c ld isn't specified,
        \c ld=n.

        \warning This is a reference to the memory. No entry is
        copied.
    */
    XymMatrixVC(double* d, int n, int m, int ld=0);


    //! Returns element (i,j)
    /*! Range checking is performed, if the code is neither
        compiled with \c NDEBUG nor with \c NORANGECHECKING.
    */
    double& operator() (int i, int j) {
        XYMRANGECHECK ((0<=i) && (i<_rowDim) && (0<=j) && (j<_colDim));
        return _d[j*_leadingDimension+i];
    }

    //! Returns element (i,j)
    /*! Range checking is performed, if the code is neither
        compiled with \c NDEBUG nor with \c NORANGECHECKING.
    */
    const double& operator() (int i, int j) const {
        XYMRANGECHECK ((0<=i) && (i<_rowDim) && (0<=j) && (j<_colDim));
        return _d[j*_leadingDimension+i];
    }

    //! Returns element (i,j)
    /*! Range checking is performed, if the code is neither
        compiled with \c NDEBUG nor with \c NORANGECHECKING.
    */
    double& get (int i, int j) {
        XYMRANGECHECK ((0<=i) && (i<_rowDim) && (0<=j) && (j<_colDim));
        return _d[j*_leadingDimension+i];
    }

    //! Returns element (i,j)
    /*! Range checking is performed, if the code is neither
        compiled with \c NDEBUG nor with \c NORANGECHECKING.
    */
    const double& get (int i, int j) const {
        XYMRANGECHECK ((0<=i) && (i<_rowDim) && (0<=j) && (j<_colDim));
        return _d[j*_leadingDimension+i];
    }

    //! Returns element (i,j) without performing range checking
    double& getNC (int i, int j) {
        return _d[j*_leadingDimension+i];
    }

    //! Returns element (i,j) without performing range checking
    const double& getNC (int i, int j) const {
        return _d[j*_leadingDimension+i];
    }

    //! Returns the number of rows of the matrix.
    int rows() const {
        return _rowDim;
    };

    //! Returns the number of cols of the matrix
    int cols() const {
        return _colDim;
    };

    //! Returns, whether this matrix is valid (having a non NULL) pointer
    bool isValid() const{
        return _d!=NULL;
    };
    

    //! Returns, whether it is a square matrix.
    bool isSquare() const{
        return (_rowDim == _colDim);
    };


    //! Returns a view of a submatrix.
    /*! The submatrix starts at row \c i, column \c j with dimensions
        \c n times \c m.
     */
    const XymMatrixVC submatrix (int i, int j, int n, int m) const {
            XYMRANGECHECK ((0<=i) && (i+n<=_rowDim) && (0<=n));
            XYMRANGECHECK ((0<=j) && (j+m<=_colDim) && (0<=m));  
            return XymMatrixVC(_d+i+j*_leadingDimension, n, m, _leadingDimension);
        }

    //! Returns a view of a submatrix.
    /*! The submatrix starts at row \c i, column \c j with dimensions
        \c n times \c m.
     */
    XymMatrixVC submatrix (int i, int j, int n, int m) {
            XYMRANGECHECK ((0<=i) && (i+n<=_rowDim) && (0<=n));
            XYMRANGECHECK ((0<=j) && (j+m<=_colDim) && (0<=m));  
            return XymMatrixVC(_d+i+j*_leadingDimension, n, m, _leadingDimension);
        }
    

    //! Returns a view of column \c j
    XymVectorV col (int j) {
            XYMRANGECHECK ((0<=j) && (j<_colDim));
            return XymVectorV (_d+j*_leadingDimension, _rowDim, 1);
        }

    //! Returns a view of column \c j
    const XymVectorV col (int j) const {
            XYMRANGECHECK ((0<=j) && (j<_colDim));
            return XymVectorV (_d+j*_leadingDimension, _rowDim, 1);
        }
    
    //! Returns a view of row \c i
    XymVectorV row (int i) {
            XYMRANGECHECK ((0<=i) && (i<_colDim));
            return XymVectorV (_d+i, _rowDim, _leadingDimension);
        }

    //! Returns a view of row \c i
    const XymVectorV row (int i) const {
            XYMRANGECHECK ((0<=i) && (i<_colDim));
            return XymVectorV (_d+i, _rowDim, _leadingDimension);
        }


    //! Returns a view of the matrices (sub-) diagonal
    /*! If \c \c d is \c 0, the view refers to the diagonal, if \c <0
        to a subdiagonal or to a superdiagonal if \c >0.
    */
    XymVectorV diag (int d=0) {
            if (d<0) return XymVectorV (_d-d, min(_rowDim+d, _colDim), 1+_leadingDimension);
            else return XymVectorV (_d+d*_leadingDimension, min(_rowDim, _colDim-d), 1+_leadingDimension);
        }

    //! Returns a view of the matrices (sub-) diagonal
    /*! If \c \c d is \c 0, the view refers to the diagonal, if \c <0
        to a subdiagonal or to a superdiagonal if \c >0.
    */
    const XymVectorV diag (int d=0) const {
            if (d<0) return XymVectorV (_d-d, min(_rowDim+d, _colDim), 1+_leadingDimension);
            else return XymVectorV (_d+d*_leadingDimension, min(_rowDim, _colDim-d), 1+_leadingDimension);
        }

    //! Returns a view of row \c i
    /*! \warning It is possible to write \c A[i][j] for the element \c
        (i,j) of matrix \c A. However this is (at least without
        optimizer) \c em much slower than using \c A(i,j), since a
        temporary row-view object has to be created.
    */
    XymVectorV operator[] (int i) 
        {
            XYMRANGECHECK ((0<=i) && (i<_colDim));
            return XymVectorV (_d+i, _rowDim, _leadingDimension);
        }

    //! Returns a view of row \c i
    /*! \warning It is possible to write \c A[i][j] for the element \c
        (i,j) of matrix \c A. However this is (at least without
        optimizer) \c em much slower than using \c A(i,j), since a
        temporary row-view object has to be created.
    */
    const XymVectorV operator[] (int i) const
        {
            XYMRANGECHECK ((0<=i) && (i<_colDim));
            return XymVectorV (_d+i, _rowDim, _leadingDimension);
        }


    //! Type of a 2 times 2 matrix
    typedef double Matrix2x2[2][2];

    //! Type of a 2 times 3 matrix
    typedef double Matrix2x3[2][3];

    //! Type of a 3 times 2 matrix
    typedef double Matrix3x2[3][2];

    //! Type of a 3 times 3 matrix
    typedef double Matrix3x3[3][3];

    //! Type of a 2 vector
    typedef double Vector2[2];

    //! Type of a 3 vector
    typedef double Vector3[3];
    

    //! Stores matrix \c a into \c *this beginning with entry \c (i,j).
    /*! \warning Untested
    */
    void store (const XymMatrixVC& a, int i, int j);

    //! Stores 2*2 matrix \c a into \c *this beginning with entry \c (i,j).
    inline void store (const Matrix2x2& a, int i, int j);

    //! Stores 2*3 matrix \c a into \c *this beginning with entry \c (i,j).
    inline void store (const Matrix2x3& a, int i, int j);

    //! Stores 3*2 matrix \c a into \c *this beginning with entry \c (i,j).
    inline void store (const Matrix3x2& a, int i, int j);

    //! Stores 3*3 matrix \c a into \c *this beginning with entry \c (i,j).
    inline void store (const Matrix3x3& a, int i, int j);
    
    //! Stores a 3 col vector \c a into \c *this beginning with entry \c (i,j)
    inline void storeCol (const Vector3& a, int i, int j);

    //! Stores a 2 col vector \c a into \c *this beginning with entry \c (i,j)
    inline void storeCol (const Vector2& a, int i, int j);

    //! Stores a 3 row vector \c a into \c *this beginning with entry \c (i,j)
    inline void storeRow (const Vector3& a, int i, int j);

    //! Stores a 2 row vector \c a into \c *this beginning with entry \c (i,j)
    inline void storeRow (const Vector2& a, int i, int j);
    

    //! Extracts matrix \c a from \c *this beginning with entry \c (i,j)
    /*! \warning Untested
    */
    void extract (XymMatrixVC& a, int i, int j) const;

    //! Extracts 2*2 matrix \c a from \c *this beginning with entry \c (i,j)
    inline void extract (Matrix2x2& a, int i, int j) const;

    //! Extracts 2*3 matrix \c a from \c *this beginning with entry \c (i,j)
    inline void extract (Matrix2x3& a, int i, int j) const;

    //! Extracts 3*2 matrix \c a from \c *this beginning with entry \c (i,j)
    inline void extract (Matrix3x2& a, int i, int j) const;

    //! Extracts 3*3 matrix \c a from \c *this beginning with entry \c (i,j)
    inline void extract (Matrix3x3& a, int i, int j) const;

    //! Extracts 3 column vector \c a from \c *this beginning with entry \c (i,j)
    inline void extractCol (Vector3& a, int i, int j) const;

    //! Extracts 3 row vector \c a from \c *this beginning with entry \c (i,j)
    inline void extractRow (Vector3& a, int i, int j) const;
    

    //! Extracts elements form \c A using a row and a column index.
    /*! Sets \c *this(i, j)=A[idxI[i], idxJ[j]]
     */
    void extractFrom (const XymMatrixVC& A, const vector<int>& idxI, const vector<int>& idxJ);

    //! Extracts columns from \c A using an index.
    /*! Sets \c *this(i, j)=A[i, idxJ[j]]
     */
    void extractColsFrom (const XymMatrixVC& A, const vector<int>& idxJ);

    //! Extracts rows from \c idxI using an index.
    /*!
      Sets \c *this(i, j)=A[idxI[i], j]
    */
    void extractRowsFrom (const XymMatrixVC& A, const vector<int>& idxJ);

    //! Returns the transpose of this matrix
    /*! Note, that using transpose(A,B) (\c xymOperations.h) is faster. */
    XymMatrixC transpose () const;    


    //! Asserts, that \c this is valid and an \c m times \c n matrix
    inline void assertFormat (int m, int n) const;

    //! Returns, whether \c this contains only finite numbers
    bool isFinite() const;

    //! Copies all data from 'm'.
    void copyFrom (const XymMatrixVC& m); 

    //! Fills the matrix with \c alpha and \c beta on the diagonal
    void fill (double alpha, double beta);

    //! Prints the matrix
    void print () const;

    //! stores the vector in the stream
    friend ostream& operator<<(ostream& os, const XymMatrixVC& v);

    //! loads the vector from the stream.
    friend istream& operator>>(ifstream& ifs, XymMatrixVC& v) throw (runtime_error);



    //************************************************************
    /*!\name Direct Access Methods

       It is often much faster to use direct access to the data via
       'double*' pointer.  As this matrix implements a concretely
       defined memory model (column major, dense rectangle as LAPACK)
       it is allowed to directly access the elements via pointers
       returned by the following functions.  However we \em strongly
       \em advice to limit this type of access to the inner loop of
       computational important algorithms, to keep encapsulation as
       strong as possible.

       The pointer for element \c (i,j) is \c
       base()+i*rowOfs()+j*colOfs().


       The parameter to pass to LAPACK are:
       \c base(), \c colOfs(), \c rows(), \c cols()

       It is most advisable to use Level 2 or Level 3 BLAS functions
       contained in the file \c xyLAPACK.h
     */
     //@{


    //! Returns a pointer to the first element
    double* base() const {return _d;};

    //! Returns a pointer to the first element of the \c i -th row
    double* rowBase(int i) const {
        XYMRANGECHECK(0<=i && i<_rowDim);
        return _d+i;
    };

    //! Returns a pointer to the first element of the \c i -th col
    double* colBase(int j) const {
        XYMRANGECHECK(0<=j && j<_colDim);
        return _d+j*_leadingDimension;
    };
    

    //! Returns the pointer increment for the next row (i.e. traversing
    //! a column).
    /*! \c rowOfs() is always 1 for a column major matrix, however it
        is recommended to use 'rowOfs()' to remain portable to column
        major matrices. (Not yet implemented)
    */
    int rowOfs() const {
        return 1;
    };


    //! Returns the pointer increment for the next column (i.e. traversing
    //! a row)
    int colOfs() const {
        return _leadingDimension;
    };

    //! Special function for LAPACK/BLAS. 
    /*! Returns a pointer to the integer member variable holding the leading
        dimension of the matrix. This pointer can directly be passed to
        FORTRAN. 
    */
    int* ld() const {
            return (int*) &_leadingDimension;
        }
    

    //! Function to construct efficient loops
    /*! Most linear algebra operations perform some elementwise
        operation on some vectors or rows/column of a matrix in the inner
        loop. To speed up that loop, one can implement it using pointers
        that are initialized by 'loop(...)'  like in the following sheme:
     
        \code
        double *p, *pEnd;
        int inc;
        v.loopRow (j, p, pEnd, inc); // or loopCol (...)
        for (;p!=pEnd;p+=incr) {
           // operations on a single matrix element
        }
        \endcode
    */
    void loopRow (int i, double*& p, double*& pEnd, int& increment) const
        {
            p = rowBase (i);
            increment = _leadingDimension;
            pEnd = p+_colDim*_leadingDimension;
        }

    //! See \c loopRow
    void loopRow (int i, double*& p, int& increment) const
        {
            p = rowBase (i);
            increment = _leadingDimension;
        }

    //! See \c loopRow
    void loopCol (int j, double*& p, double*& pEnd, int& increment) const
        {
            p = colBase (j);
            increment = 1;
            pEnd = p+_rowDim;
        }

    void loopCol (int j, double*& p, double*& pEnd) const
        {
            p = colBase (j);
            pEnd = p+_rowDim;
        }

    //! See \c loopRow
    void loopCol (int j, double*& p, int& increment) const
        {
            p = colBase (j);
            increment = 1;
        }
//@}

protected:
    //! Pointer to entry (0,0)
    double* _d;

    //! Numbers of rows in the matrix
    int _rowDim; 
    //! Numbers of columns in the matrix
    int _colDim;

    //! How much to add to 'd' to get one column/row further
    int _leadingDimension;
};

//! Checks, whether \c a*b is valid.
inline void checkForMultiply (const XymMatrixVC& a, const XymMatrixVC& b);

//! Checks, whether \c a*b is valid and the result is compatible to \c res.
inline void checkForMultiply (const XymMatrixVC& res, const XymMatrixVC& a, const XymMatrixVC& b);

//! Checks, whether \c a*v is valid and the result is compatible to \c res.
inline void checkForMultiply (const XymMatrixVC& a, const XymVectorV& v);

//! Checks, whether \c a*v is valid
inline void checkForMultiply (const XymVectorV& res, const XymMatrixVC& a, const XymVectorV& v);

//! Checks, whether \c a and \c b have the same format
inline void checkSameFormat (const XymMatrixVC& a, const XymMatrixVC& b);

//! Checks, whether \c a is square
inline void checkSquare (const XymMatrixVC& a);

//! copies the matrix
void copyFrom (XymMatrixVC& res, const XymMatrixVC& a);


//! Prints the matrix
/*! This function can be used in the debugger to view the
    content of a matrix. */
void print (const XymMatrixVC& m);

//! Store matrix into stream
ostream& operator<<(ostream& os, const XymMatrixVC& v);

//! Load matrix from stream
/*!\warning AXymMatrixVC can not change it's size, thus the vector specified
  in the stream must have the correct number of elements.
  If the vectors size is a-priori unknown, use \c XymVector::operator>>
*/
istream& operator>>(ifstream& ifs, XymMatrixVC& v) throw (runtime_error);

//! Computes A+B
/*! Keep in mind, that C=A+B is less efficient than \c add (C, A, B) (\c xymOperations.h)
 */
XymMatrixC const operator + (const XymMatrixVC& A, const XymMatrixVC& B);

//! Computes A=A+B
/*! Keep in mind, that C=A+B is less efficient than \c add (C, A, B) (\c xymOperations.h)
 */
void operator += (XymMatrixVC& A, const XymMatrixVC& B);

//! Computes A-B
/*! Keep in mind, that C=A-B is less efficient than \c sub (C, A, B) (\c xymOperations.h)
 */
XymMatrixC const operator - (const XymMatrixVC& A, const XymMatrixVC& B);

//! Computes A=A-B
/*! Keep in mind, that C=A-B is less efficient than \c sub (C, A, B) (\c xymOperations.h)
 */
void operator -= (XymMatrixVC& A, const XymMatrixVC& B);

//! Computes A*B
/*! Keep in mind, that C=A*B is less efficient than \c multiply (C, A, B) (\c xymOperations.h)
 */
XymMatrixC const operator * (const XymMatrixVC& A, const XymMatrixVC& B);

//! Computes A*v
/*! Keep in mind, that v2=A*v is less efficient than \c multiply (v2, A, v) (\c xymOperations.h)
 */
XymVector const operator * (const XymMatrixVC& A, const XymVectorV& v);

//! Computes A*labmda
/*! Keep in mind, that B=A*labmda is less efficient than \c multiply (B, A, lambda) (\c xymOperations.h)
 */
XymMatrixC const operator * (const XymMatrixVC& A, double lambda);

//! Computes A=A*labmda
/*! Keep in mind, that B=A*labmda is less efficient than \c multiply (B, A, lambda) (\c xymOperations.h)
 */
void operator *= (XymMatrixVC& A, double lambda);

//! Computes lambda*A=A*lambda
/*! Keep in mind, that B=A*labmda is less efficient than \c multiply (B, A, lambda) (\c xymOperations.h)
 */
XymMatrixC const operator * (double lambda, const XymMatrixVC& A);



#include "xymMatrixVC_i.h"

#endif
