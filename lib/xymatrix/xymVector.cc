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

/*!\file xymVector.cc Contains the implementation for class \c
  XymVector representing a vector.
*/

#include "xymVector.h"

XymVector::XymVector()
        :XymVectorV(), _memBase(NULL), _nReserved(0), _autoGrow(true)
{}


XymVector::XymVector (int n, bool initialize)
        :XymVectorV(new double[n], n, 1), _nReserved(n), _autoGrow(true)
{
    _memBase = _d;
    if (initialize) for (int i=0; i<n; i++) _d[i]=0;
}


XymVector::XymVector (int n, int nReserve, bool initialize, bool autoGrow)
        :XymVectorV (new double[nReserve], n, 1), _nReserved(nReserve), _autoGrow(autoGrow)
{
    XYMASSERT (n<=nReserve, "reserved space exceeded");
    _memBase = _d;
    if (initialize) for (int i=0; i<n; i++) _d[i]=0;
}


XymVector::XymVector(double* d, int n, int stride)
        :XymVectorV (new double[n], n, 1), _nReserved(n), _autoGrow(true)
{
    _memBase = _d;
    double* p = d;
    for (int i=0; i<n; i++) {
        _d[i] = *p;
        p+=stride;
    }
}

    
XymVector::XymVector(const XymVectorV& v)
        :XymVectorV (new double[v.size()], v.size()), _nReserved(v.size()), _autoGrow(true)
{
    _memBase = _d;
    for (int i=0; i<_n; i++) _d[i] = v[i];
}


XymVector::XymVector(const XymVector& v)
        :XymVectorV (new double[v.size()], v.size()), _nReserved(v.size()), _autoGrow(true)
{
    _memBase = _d;
    for (int i=0; i<_n; i++) _d[i] = v[i];
}


XymVector::~XymVector()
{
    clear();
}


void XymVector::setAutoGrow (bool autoGrow)
{
    _autoGrow = autoGrow;
}


void XymVector::reserve (int n, bool autoGrow)
{
    XYMASSERT (n>=_n, "too little space reserved");
    if (n>_nReserved) {
        double* d  = new double[n];
        for (int i=0; i<_n; i++) d[i] = _d[i];
        if (_memBase!=NULL) delete[] _memBase;
        _memBase = d;
        _nReserved = n;
        _d = d;
        _stride = 1;
    }
    _autoGrow = autoGrow;
}


    
void XymVector::append (int n, bool initialize)
{
    XYMASSERT (n>=0, "n must be >=0");
    if (_autoGrow) reserve (n+_n, _autoGrow);
    XYMASSERT (_n+n<=_nReserved, "reserved space exceeded");
    if (initialize) for (int i=0; i<n; i++) _d[i]=0;
    _n += n;
}

    
void XymVector::deleteLast (int n)
{
    XYMASSERT (0<=n && n<=_n, "index error");
    _n -= n;
}


void XymVector::insert (int i, int n, bool initialize)
{
    XYMASSERT (n>=0, "n must be >=0");
    if (_autoGrow) reserve (n+_n, _autoGrow);
    XYMASSERT (_n+n<=_nReserved, "reserved space exceeded");
    for (int k=_n-1; k>=i; k--) _d[k+n] = _d[k];
    if (initialize) for (int k=i; k<i+n; k++) _d[k]=0;
    _n += n;
}


void XymVector::remove (int i, int n)
{
    XYMASSERT ((0<=i) && (i+n)<=_n, "illegal index or length");
    int start = i+n;
    int rem   = _n-i-n;
    for (int k=0; k<rem; k++) _d[i+k] = _d[start+k];
    _n -= n;
}


void XymVector::resize (int n, bool initialize)
{
    if (n<_n) _n = n; // We do not change '_nReserved'
    else if (n>_n) {
        if (_autoGrow) reserve (n, _autoGrow);
        if (initialize) for (int i=n; i<_n; i++) _d[i]=0;
        _n = n;
    }
}


int XymVector::memoryUsage() const
{
    return _nReserved*sizeof(double);
}


void XymVector::clear()
{
    if (_memBase!=NULL) delete[] _memBase;
    _memBase = _d = NULL;
    _n = _nReserved = _stride = 0;
    _autoGrow = true;
}


void XymVector::create(int n, bool init, int reserve)
{
    clear();
    if (reserve==0) reserve=n;
    XYMASSERT(reserve>=n, "reserve must be 0 or >= n");
    
    _memBase = _d = new double[reserve];
    _n = n;
    _nReserved = reserve;
    _stride = 1;
    _autoGrow = true;
    if (init) for (int i=0; i<n; i++) (*this)[i] = 0;
}


void XymVector::extractFrom (const XymVectorV& v, const vector<int>& idx)
{
    if (!isValid()) create (idx.size(), false);
    XymVectorV::extractFrom (v, idx);
}


void XymVector::copyFrom (const XymVectorV& v)
{
    if (!v.isValid()) clear();
    else {
        if (!isValid() || size()!=v.size()) create (v.size(), false, 0);
        XymVectorV::copyFrom (v);
    }
}


void XymVector::transferFrom (XymVector& v)
{
    clear();
    if (v.isValid()) {
        _memBase = v._memBase;
        _nReserved = v._nReserved;
        _autoGrow = v._autoGrow;
        _d = v._d;
        _n = v._n;
        _stride = v._stride;
        v._memBase = v._d = NULL;
        v._n = v._nReserved = v._stride = 0;
        v._autoGrow = true;
    }
}


XymVector& XymVector::operator = (const XymVector& v)
{
    copyFrom (v);
    return *this;
}


XymVector& XymVector::operator = (const XymVectorV& v)
{
    copyFrom (v);
    return *this;
}

istream& operator>>(ifstream& ifs, XymVector& v) throw (runtime_error)
{

    xymExpect (ifs, "{");
    v.clear();
    string buffer;
    ifs >> buffer;
    if (buffer!="}") {
        do {
            double d;
            if (sscanf(buffer.c_str(), "%lf", &d)<1) return ifs;
            v.append(1, false);
            v[v.size()-1] = d;
            ifs >> buffer;
            if (buffer==",") ifs >> buffer;
            else break;
        } while (true); // break
    }
    return ifs;
}

