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

/*!\file xymMatrixC.cc Contains implementation for \c XymMatrixVC.
*/

#include "xymMatrixC.h"
#include "xymOperations.h"
#include <vector>


XymMatrixC::XymMatrixC()
        :XymMatrixVC(), _memBase(NULL), _rowReserved(0), _colReserved(0), _autoGrow(true)
{}


XymMatrixC::XymMatrixC(int n, int m, bool initialize)
        :XymMatrixVC (new double[n*m], n, m), _rowReserved(n), _colReserved(m), _autoGrow(true)
{
    _memBase = _d;
    if (initialize) {
        double *d = _d, *dEnd = d+n*m;
        for (;d!=dEnd; d++) *d = 0;
    }
}

    
XymMatrixC::XymMatrixC (const XymMatrixVC& a00, const XymMatrixVC& a01, const XymMatrixVC& a10, const XymMatrixVC& a11)
        :XymMatrixVC(), _memBase(NULL), _rowReserved(0), _colReserved(0), _autoGrow(true)    
{
    block (a00, a01, a10, a11);
}


XymMatrixC::XymMatrixC(const XymMatrixC& a)
        :XymMatrixVC(), _memBase(NULL), _rowReserved(0), _colReserved(0), _autoGrow(true)
{
    copyFrom(a);
}


XymMatrixC::XymMatrixC(const XymMatrixVC& a)
        :XymMatrixVC(), _memBase(NULL), _rowReserved(0), _colReserved(0), _autoGrow(true)
{
    copyFrom(a);
}
 

XymMatrixC::XymMatrixC (const XymMatrixVC* a, int n, int m)
        :XymMatrixVC(), _memBase(NULL), _rowReserved(0), _colReserved(0), _autoGrow(true)    
{
    block (a, n, m);
}


XymMatrixC::XymMatrixC (const Matrix2x2& A)
        :XymMatrixVC(), _memBase(NULL), _rowReserved(0), _colReserved(0), _autoGrow(true)
{
  create (2, 2);
  _d[0] = A[0][0];
  _d[1] = A[1][0];
  _d[2] = A[0][1];
  _d[3] = A[1][1];  
}


XymMatrixC::XymMatrixC (const Matrix2x3& A)
        :XymMatrixVC(), _memBase(NULL), _rowReserved(0), _colReserved(0), _autoGrow(true)
{
  create (3, 2);
  _d[0] = A[0][0];
  _d[1] = A[1][0];
  _d[2] = A[0][1];
  _d[3] = A[1][1];  
  _d[4] = A[0][2];
  _d[5] = A[1][2];  
}


XymMatrixC::XymMatrixC (const Matrix3x2& A)
        :XymMatrixVC(), _memBase(NULL), _rowReserved(0), _colReserved(0), _autoGrow(true)
{
  create (2, 3);
  _d[0] = A[0][0];
  _d[1] = A[1][0];
  _d[2] = A[2][0];  
  _d[3] = A[0][1];
  _d[4] = A[1][1];  
  _d[5] = A[2][1];  
}
    

XymMatrixC::XymMatrixC (const Matrix3x3& A)
        :XymMatrixVC(), _memBase(NULL), _rowReserved(0), _colReserved(0), _autoGrow(true)
{
  create (3, 3);
  _d[0] = A[0][0];
  _d[1] = A[1][0];
  _d[2] = A[2][0];  
  _d[3] = A[0][1];
  _d[4] = A[1][1];  
  _d[5] = A[2][1];  
  _d[6] = A[0][2];
  _d[7] = A[1][2];  
  _d[8] = A[2][2];  
}


XymMatrixC::~XymMatrixC()
{
    clear();
}


void XymMatrixC::create (int n, int m, bool initialize)
{
    if (_memBase!=NULL && (_rowReserved<n || _colReserved<m)) {
        delete[] _memBase;
        _memBase = NULL;
    }
    if (_memBase==NULL) {
        _memBase = new double[n*m];
        _rowReserved = n;
        _colReserved = m;
    }
    _d = _memBase;
    _rowDim = n;
    _colDim = m;
    _leadingDimension = _rowReserved;
    _autoGrow = true;
    if (initialize) {
        double *d = _d, *dEnd = d+_leadingDimension*m;
        for (;d!=dEnd; d++) *d = 0;
    }
}


void XymMatrixC::createLike (const XymMatrixVC& A, bool initialize)
{
    create (A.rows(), A.cols(), initialize);
}


void XymMatrixC::block (const XymMatrixVC& a00, const XymMatrixVC& a01, const XymMatrixVC& a10, const XymMatrixVC& a11)
{
    XymMatrixVC a[] = {a00, a01, a10, a11};
    block (a, 2, 2);
}


void XymMatrixC::block (const XymMatrixVC* a, int n, int m)
{
    int cols=0;
    for (int j=0; j<m; j++) {
        cols += a[j].cols();
        for (int i=0; i<n; i++) {
            XYMASSERT (a[i*m+j].cols()==a[i].cols(), "all block of a block column must have the same column size");
        }
    }
    int rows=0;
    for (int i=0; i<n; i++) {
        rows += a[i*m].rows();
        for (int j=0; j<n; j++) {
            XYMASSERT (a[i*m+j].rows()==a[i*m].rows(), "all block of a block row must have the same row size");
        }
    }
    create (rows, cols);
    int storeI=0, storeJ;
    for (int i=0; i<n; i++) {
        storeJ=0;
        for (int j=0; j<n; j++) {
            store (a[i*m+j], storeI, storeJ);
            storeJ += a[j].cols();
        }
        XYMASSERT (storeJ==cols, "sum of rows must equal overall matrix rows");
        storeI+=a[i*m].rows();
    }
    XYMASSERT (storeI==rows, "sum of cols must equal overall matrix cols");
}


void XymMatrixC::clear()
{
    if (_memBase!=NULL) delete[] _memBase;
    _memBase = _d = NULL;
    _colDim = _colReserved = 0;
    _rowDim = _rowReserved = 0;
    _autoGrow = true;
    _leadingDimension = 0;
}


void XymMatrixC::reserve (int n, int m, bool autoGrow)
{
    if (n>_rowReserved || m>_colReserved) {
        XYMASSERT (n>=_rowDim && m>=_colDim, "argument too small");
        XymMatrixVC oldM(*this);
        double* oMB = _memBase;
        _memBase = _d = new double[n*m];
        _leadingDimension = n;
        _rowReserved = n;
        _colReserved = m;
        _autoGrow    = autoGrow;
        for (int i=0; i<_rowDim; i++) for (int j=0; j<_colDim; j++)
            (*this)(i,j) = oldM(i,j);
        if (oMB!=NULL) delete[] oMB; // oldM gets invalid
    }
}

    
void XymMatrixC::setAutoGrow (bool autoGrow)
{
    _autoGrow = autoGrow;
}


void XymMatrixC::appendColumn (int m, bool initialize)
{
    appendRowAndColumn (0, m, initialize);
}


void XymMatrixC::appendRow    (int n, bool initialize)
{
    appendRowAndColumn (n, 0, initialize);
}


void XymMatrixC::appendRowAndColumn (int n, int m, bool initialize)
{
    XYMASSERT (0<=m && 0<=n, "n,m must be >=0");
    if (_autoGrow) reserve (_rowDim+n, _colDim+m, _autoGrow);
    XYMASSERT (_rowDim+n<=_rowReserved && _colDim+m<=_colReserved, "reserved space exceeded");
    int oR = _rowDim;
    int oC = _colDim;
    _rowDim += n;
    _colDim += m;
    if (initialize) {
      double *d, *dEnd;
      for (int j=0; j<_colDim; j++) {
        loopCol (j, d, dEnd);
        if (j<oC) d+=oR;        
        for (;d!=dEnd; d++) *d = 0;
      }
    }
}


void XymMatrixC::deleteLastColumn (int m)
{
    deleteLastRowAndColumn (0, m);
}


void XymMatrixC::deleteLastRow (int n)
{
    deleteLastRowAndColumn (n, 0);
}


void XymMatrixC::deleteLastRowAndColumn(int n, int m)
{
    XYMASSERT ((0<=n) && (n<=_rowDim), "index error");
    XYMASSERT ((0<=m) && (m<=_colDim), "index error");
    _rowDim -= n;
    _colDim -= m;
}


void XymMatrixC::insertColumn (int j, int m, bool initialize)
{
    insertRowAndColumn (0, j, 0, m, initialize);
}


void XymMatrixC::insertRow    (int i, int n, bool initialize)
{
    insertRowAndColumn (i, 0, n, 0, initialize);
}


void XymMatrixC::insertRowAndColumn (int insI, int insJ, int n, int m, bool initialize)
{
    XYMASSERT (0<=m && 0<=n, "index error");
    if (_autoGrow) reserve (_rowDim+n, _colDim+m, _autoGrow);
    XYMASSERT (_rowDim+n<=_rowReserved && _colDim+m<=_colReserved, "reserve space exceeded");
    int oR = _rowDim;
    int oC = _colDim;
    _rowDim += n;
    _colDim += m;
    for (int i=0; i<insI; i++) for (int j=oC-1; j>=insJ; j--) (*this)(i, j+m) = (*this)(i, j);
    for (int i=oR-1; i>=insI; i--) for (int j=0; j<insJ; j++) (*this)(i+n, j) = (*this)(i, j);
    for (int i=oR-1; i>=insI; i--) for (int j=oC-1; j>=insJ; j--) (*this)(i+n, j+m) = (*this)(i, j);
    if (initialize) {
        for (int i=0; i<_rowDim; i++) for (int j=insJ; j<insJ+m; j++) (*this)(i, j) = 0;
        for (int i=insI; i<insI+n; i++) for (int j=0; j<_colDim; j++) (*this)(i, j) = 0;
    }
}


void XymMatrixC::deleteRow (int i, int n)
{
    deleteRowAndColumn (i, 0, n, 0);
}


void XymMatrixC::deleteColumn (int j, int m)
{
    deleteRowAndColumn (0, j, 0, m);
}


void XymMatrixC::deleteRowAndColumn (int delI, int delJ, int n, int m)
{
    XYMASSERT ((0<=delI) && (0<=n) && (delI+n<=_rowDim), "index error (delI)");
    XYMASSERT ((0<=delJ) && (0<=m) && (delJ+m<=_colDim), "index error (delJ)");
    int nn = _rowDim-delI-n; // how many rows below the deleted rows
    int mm = _colDim-delJ-m; // how many cols right the deleted cols
    for (int i=delI; i<delI+nn; i++) for (int j=0; j<delJ; j++) (*this)(i,j) = (*this)(i+n,j);
    for (int i=0; i<delI; i++) for (int j=delJ; j<delJ+mm; j++) (*this)(i,j) = (*this)(i,j+m);
    for (int i=delI; i<delI+nn; i++) for (int j=delJ; j<delJ+mm; j++) (*this)(i,j) = (*this)(i+n,j+m);
    _colDim -= n;
    _rowDim -= m;
}


XymMatrixC::XymMatrixC (const XymMatrixVC& a, double entry)
        :XymMatrixVC(), _memBase(NULL), _rowReserved(0), _colReserved(0), _autoGrow(true)
{
    create (a.rows(), a.cols());
    for (int j=0; j<cols(); j++) {
        double* dst = &getNC (0,j);
        for (int i=0; i<rows(); i++) {
            *dst = entry;
            dst ++;
        }
    }
}

XymMatrixC::XymMatrixC(double* d, int n, int m, bool columnMajor, int ld)
        :XymMatrixVC(), _memBase(NULL), _rowReserved(0), _colReserved(0), _autoGrow(true)
{
    create (n, m);
    if (columnMajor && ld<n) ld = n;
    if (!columnMajor && ld<m) ld = m;
    for (int j=0; j<m; j++) {
        double* src;
        int srcIncr;
        if (columnMajor) {
            src     = d+ld*j;
            srcIncr = 1;
        }
        else {
            src     = d+j;
            srcIncr = ld;
        }
        double* dst = &getNC (0, j);
        for (int i=0; i<n; i++) {
            double xx = *src;
            *dst = xx;
            
            src+=srcIncr;
            dst++;
        }
    }
}

    
int XymMatrixC::memoryUsage() const
{
  return _rowReserved*_colReserved*sizeof(double);
}


void XymMatrixC::copyFrom (const XymMatrixVC& m)
{
    ::copyFrom (*this, m);
}


void copyFrom (XymMatrixC& res, const XymMatrixVC& a)
{
    if (!a.isValid())  res.clear();
    else {
        if (!res.isValid() || res.rows()!=a.rows() || res.cols()!=a.cols()) 
            res.create (a.rows(), a.cols());
        copyFrom ((XymMatrixVC&)res, a);
    }
}


void XymMatrixC::transferFrom (XymMatrixC& m)
{
    clear();
    if (m.isValid()) {
        _memBase = m._memBase;
        _rowReserved = m._rowReserved;
        _colReserved = m._colReserved;
        _d = m._d;
        _rowDim = m._rowDim;
        _colDim = m._colDim;
        _leadingDimension = m._leadingDimension;
        m._memBase = m._d = NULL;
        m._colDim = m._colReserved = 0;
        m._rowDim = m._rowReserved = 0;
        m._autoGrow = true;
        m._leadingDimension = 0;
    }
}


void XymMatrixC::extractFrom (const XymMatrixVC& A, const vector<int>& idxI, const vector<int>& idxJ)
{
    if (!isValid()) create (idxI.size(), idxJ.size(), false);
    XymMatrixVC::extractFrom (A, idxI, idxJ);
}



XymMatrixC& XymMatrixC::operator = (const XymMatrixVC& A)
{
    copyFrom (A);
    return *this;
}


XymMatrixC& XymMatrixC::operator = (const XymMatrixC& A)
{
    copyFrom (A);
    return *this;
}


istream& operator>>(ifstream& ifs, XymMatrixC& v) throw (runtime_error)
{
    vector<double> dBuf;
    xymExpect (ifs, "{");
    string buffer;
    ifs >> buffer;
    int colDim = -1;
    int rowCtr=0;
    while (buffer!="}") {
        rowCtr++;
/*
        if (buffer!="{")  ERR_THROW(GenFileIsCorrupted("(istream)", 
                                                        "missing  { found " + buffer + " instead", "operator>>(ifstream&, XymMatrixVC&)"));
*/
        if (buffer!="{") throw runtime_error ("XymMatrixVC::operator>>: missing  { found " + buffer + " instead");
        int colCtr=0;
        ifs >> buffer;
        while (buffer!="}") {
            double d;
            d = atof(buffer.c_str());
            colCtr++;
            if (colDim>=0 && colCtr>colDim) throw runtime_error ("XymMatrixVC::operator>>: too many columns in one row");
            dBuf.push_back (d);
            ifs >> buffer;
        }
        if (colDim>=0 && colCtr<colDim) throw runtime_error ("XymMatrixVC::operator>>: too few columns in one row");
        if (colDim==-1) colDim = colCtr;
        ifs >> buffer;
    }
    // Now allocate memory and copy the data from 'dBuf'
    v.clear();
    v.appendRowAndColumn (rowCtr, colDim);
    int k=0;
    for (int i=0; i<v.rows(); i++) for (int j=0; j<v.cols(); j++) {
        v(i,j) = dBuf[k];
        k++;
    }
    XYMASSERT (k==(int) dBuf.size() && k==colDim*rowCtr, "internal error");
    return ifs;
}


XymMatrixC const operator + (const XymMatrixVC& A, const XymMatrixVC& B)
{
  XymMatrixC result;
  add (result, A, B);
  return result;
}


void operator += (XymMatrixVC& A, const XymMatrixVC& B)
{
  add (A, B);
}


XymMatrixC const operator - (const XymMatrixVC& A, const XymMatrixVC& B)
{
  XymMatrixVC result;
  sub (result, A, B);
  return result;
}

void operator -= (XymMatrixVC& A, const XymMatrixVC& B)
{
  sub (A, B);
}


XymMatrixC const operator * (const XymMatrixVC& A, const XymMatrixVC& B)
{
  XymMatrixC result;
  multiply (result, A, B);
  return result;  
}


XymVector const operator * (const XymMatrixVC& A, const XymVectorV& v)
{
  XymVector result;
  multiply (result, A, v);
  return result;
}


XymMatrixC const operator * (const XymMatrixVC& A, double lambda)
{
  XymMatrixC result;
  multiply (result, A, lambda);
  return result;
}


void operator *= (XymMatrixVC& A, double lambda)
{
  multiply (A, A, lambda);
}


XymMatrixC const operator * (double lambda, const XymMatrixVC& A)
{
  XymMatrixC result;
  multiply (result, A, lambda);
  return result;
}
  
