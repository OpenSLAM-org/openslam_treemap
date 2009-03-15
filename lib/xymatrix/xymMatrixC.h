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

/*!\file xymMatrixC.h Contains header for \c XymMatrixVC.
*/


#ifndef XYMMATRIXC_H
#define XYMMATRIXC_H

#include "xymMatrixVC.h"

//! Matrix owning its memory 
/*! Additional to the methods of \c XymMatrixVC there are methods for
    allocation and resizing.  In contrast to \c XymMatrixVC all
    constructor copy all data.
*/
class XymMatrixC : public XymMatrixVC
{
public:
    //! Empty matrix. (Uninitialised without memory)
    /*! \c isValid() will return \c false.
     */
    XymMatrixC();

    //! n*m Matrix containing 0 if \c initialize==true.
    /*! Both \c n and \c m may be 0. Such a matrix does not contain any
      data but is initialised. (\c isValid()==true)
    */
    XymMatrixC(int n, int m, bool initialize=true);

    //! Constructs itself as a copy of 'a'. (Copying all data)
    XymMatrixC (const XymMatrixVC& a);

    //! Like above for \c XymMatrixC. 
    /*! We need to explicitly implement this version, since otherwise
        the default copy constructor takes effect. */
    XymMatrixC (const XymMatrixC& a);
  
    //! Constructs a matrix with shape like 'a' filled with 'entry'
    XymMatrixC (const XymMatrixVC& a, double entry);

    //! Contruct matrix from raw data (Copying all data)
    /*! Constructs a new matrix of dimension \c n times \c m and
        copies the data from \c *d. \c columnMajor specifies, whether
        the matrix \c *d is a column-major matrix, \c ld specifies
        it's leading dimension increment.  The default value \c
        columnMajor=false can be used to fill matrices with values
        predefinied by a C \c double* const.
    */
    XymMatrixC (double* d, int n, int m, bool columnMajor=false, int ld=0);

    //! Copy of a flat double[2][2] matrix
    XymMatrixC (const Matrix2x2& A);    

    //! Copy of a flat double[2][3] matrix
    XymMatrixC (const Matrix2x3& A);    

    //! Copy of a flat double[3][2] matrix
    XymMatrixC (const Matrix3x2& A);    

    //! Copy of a flat double[3][3] matrix
    XymMatrixC (const Matrix3x3& A);    

    //! Constructs a 2*2 block matrix : {{a00,a01},{a10,a11}}
    XymMatrixC (const XymMatrixVC& a00, const XymMatrixVC& a01, const XymMatrixVC& a10, const XymMatrixVC& a11);

    //! Construct a block matrix
    /*! Constructs a general block matrix consisting of \c n times \c
        m blocks, where block \c (i,j) is \c a[i*n+j] \warning The
        memory format for this array of pointers is row major.  This
        is to make it possible to pass an array defined by

        \code
        XymMatrixVC* blocks = {&block1, &block2, ...}
        \endcode
    */
    XymMatrixC (const XymMatrixVC* a, int n, int m);

    //! Destructor. Frees the allocated data space.
    ~XymMatrixC();

    //! Creates an \c n times \c m matrix. 
    /*! If memory has been reserved, that is taken.  If the matrix is
        initialized, the former memory is freed.  If \c initialized is
        true, the matrix is set to 0.
    */
    void create (int n, int m, bool initialize=true);

    //! Creates a matrix of same format as A
    void createLike (const XymMatrixVC& A, bool initialize=true);

    //! Frees the matrix memory, returning to uninitialized state.
    void clear();

    //! Reserves memory for being able to store an \c n times \c m matrix.
    /*! If enough memory is reserved, nothing is done. Otherwise new
        memory is allocated and the matrix is copied.  \warning: This
        operation renders \em all views generated for this matrix
        invalid.  \c autoGrow determines the behaviour, when an
        growing operation exceeds the reserved size. If \c
        autoGrow==true (default for not reserved matrices), new memory
        will be allocated and the matrix copied (rendering all matrix
        views invalid).  If \c autoGrow==false an error will be
        thrown.
    */
    void reserve (int n, int m, bool autoGrow=false);

    //! Set the autoGrow flag
    /*! If \c autoGrow==true an automatic reallocation of the memory
        is performed if a call to \c append() or \c insert(...) exceed
        the reserved memory, otherwise this is an error. Note, that this
        process renders all views invalid.
    */
    void setAutoGrow (bool autoGrow=true);

    
    //! Append \c m rows (see \c appendRowAndColumn())
    void appendColumn (int m, bool initialize=true);

    //! Append \c n columns (see \c appendRowAndColumn())
    void appendRow    (int n, bool initialize=true);

    //! Appends \c n rows and / or \c m columns to the matrix
    /*! If \c initialize==true the appended entries are initializing
        with \c 0. Those operation can be performed without copying if
        an appropriate amount of memory has been reserved by calling
        \c reserve().
    */
    void appendRowAndColumn (int n, int m, bool initialize=true);

    //! Deletes the last \c m columns. (See \c deleteLastRowAndColumn())
    void deleteLastColumn (int m);

    //! Deletes the last \c n rows. (See \c deleteLastRowAndColumn())
    void deleteLastRow (int n);

    //! Deletes the last \c n rows and or \c m columns without freeing the memory.
    /*! The memory of the columns / rows deleted remains reserved so
        new columns and rows can be appended without reallocating the
        memory.
    */
    void deleteLastRowAndColumn(int n, int m);

    //! Insert \c m column before column \c j. (See \c insertRowAndColumn())
    void insertColumn (int j, int m, bool initialize=true);

    //! Insert \c n rows before row \c i. (See \c insertRowAndColumn())
    void insertRow    (int i, int n, bool initialize=true);

    //! Inserts \c n rows and  \c m columns before row \c i and column \c j,
    /*! If possible it uses reserved space. But necessarily the
        routine has to copy the right and/or lower part of the
        matrix. (Thus, even if enough memory is reserved all views to
        the right lower part refer to the same indices but as the
        matrix is expanded not to the same elements as before) If \c
        initialize is \c true, the new entries are initialized with 0.
    */
    void insertRowAndColumn (int i, int j, int n, int m, bool initialize=true);

    //! Deletes \c n rows strating with row \c i
    void deleteRow (int i, int n);

    //! Deletes \c m column starting with column \c j
    void deleteColumn (int j, int m);

    //! Deletes \c n rows and \c m columns starting with row \c i and column \c j
    /*! The right lower part of the matrix has to be copied, so 
        all views to the right lower part refer to the same indices but as the
        matrix is shrunken not to the same elements as before)
    */
    void deleteRowAndColumn (int i, int j, int n, int m);

    //! Sets the matrix to be the block matrix consisting of '{{a00,a01},{a10,a11}}'
    void block (const XymMatrixVC& a00, const XymMatrixVC& a01, const XymMatrixVC& a10, const XymMatrixVC& a11);

    //! Sets the matrix to be a block matrix.
    /*! The block matrix has \c n times \c m blocks, where block \c
        (i,j) is \c 'a[i*m+j]' \warning: The array of blocks is
        defined in a rowmajor way.
    */
    void block (const XymMatrixVC* a, int n, int m);

    //! Sets \c *this(i, j)=A[idxI[i], idxJ[j]]
    /*! If \c this is uninitialized, it is initialize to the right format.
     */
    void extractFrom (const XymMatrixVC& A, const vector<int>& idxI, const vector<int>& idxJ);

    //! Returns the meory usage (in bytes) of this matrix
    /*! Does not include \c sizeof(*this). */
    int memoryUsage() const;

    //! Copies all data from \c m.
    void copyFrom (const XymMatrixVC& m);
    
    //! The whole data is copied. If necessary new memory is allocated.
    XymMatrixC& operator = (const XymMatrixVC& A);

    //! The whole data is copied. If necessary new memory is allocated.
    XymMatrixC& operator = (const XymMatrixC& A);

    //! Transfer all data from \c m
    /*! The data is moved from \c m to \c this by moving the pointer
        without copying. As a consequence \c m is invalid afterward
        containing no data. This function is for instance usefull to
        move computation results from a local variable to somewhere
        else, where the loca variable would be destroyed anyway.
     */
    void transferFrom (XymMatrixC& m);

protected:
    //! Pointer to the beginning of the allocated memory block.
    double* _memBase;

    //! Number of rows reserved.
    int _rowReserved;

    //!  Number of columns reserved.
    int _colReserved;
#
    //! Whether the matrix shall reallocate its memory if necessary.
    bool _autoGrow;
};

//! copies the matrix. 
/*! 'res' can be uninitialised
 */
void copyFrom (XymMatrixC& res, const XymMatrixVC& a);

#endif  /* XYMMATRIXC_H */
