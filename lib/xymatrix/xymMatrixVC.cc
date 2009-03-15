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

/*!\file xymMatrixVC.h Contains the implementation of class \c
  XymMatrixVC representing a matrix view.
*/

#include "xymMatrixVC.h"
#include "xymMatrixC.h"
#include "xymOperations.h"
#include <math.h>

XymMatrixVC::XymMatrixVC()
        :_d(NULL),  _rowDim(0), _colDim(0), _leadingDimension(0)
{}


XymMatrixVC::XymMatrixVC (double* d, int n, int m, int ld)
        :_d(d), _rowDim(n), _colDim(m), _leadingDimension(max(ld,n))
{}


void XymMatrixVC::store (const XymMatrixVC& a, int i, int j)
{
    submatrix (i, j, a.rows(), a.cols()).copyFrom(a);
}


void XymMatrixVC::extract (XymMatrixVC& a, int i, int j) const
{
    a.copyFrom (submatrix (i, j, a.rows(), a.cols()));
}


void XymMatrixVC::extractFrom (const XymMatrixVC& A, const vector<int>& idxI, const vector<int>& idxJ)
{
    int n = rows();
    int m = cols();
    XYMFORMAT (n==(int) idxI.size());
    XYMFORMAT (m==(int) idxJ.size());
    for (int i=0; i<n; i++) XYMRANGECHECK(0<=idxI[i] && idxI[i]<A.rows());
    for (int i=0; i<m; i++) XYMRANGECHECK(0<=idxJ[i] && idxJ[i]<A.cols());
    for (int j=0; j<m; j++) {
        int idJ = idxJ[j];
        double* dP = colBase(j);
        const double* sP = A.colBase(idJ);
        const int* idxP  = &idxI[0];
        const int* idxE  = idxP + n;
        while (idxP!=idxE) {
            *dP = sP[*idxP];
            dP ++;
            idxP++;
        }
    }
}


void XymMatrixVC::extractColsFrom (const XymMatrixVC& A, const vector<int>& idxJ)
{
    int m = cols();
    XYMFORMAT (m==(int) idxJ.size());
    for (int i=0; i<m; i++) XYMRANGECHECK(0<=idxJ[i] && idxJ[i]<A.cols());
    for (int j=0; j<m; j++) {
        int idJ = idxJ[j];
        double* dP = colBase(j);
        double *sP, *sE;
        int sIncr;
        A.loopCol (idJ, sP, sE, sIncr);

        while (sP!=sE) {
            *dP = *sP;
            dP++;
            sP+= sIncr;
        }
    }
}


void XymMatrixVC::extractRowsFrom (const XymMatrixVC& A, const vector<int>& idxI)
{
    int n = rows();
    int m = cols();
    XYMFORMAT (n==(int) idxI.size());
    for (int i=0; i<n; i++) XYMRANGECHECK(0<=idxI[i] && idxI[i]<A.rows());
    for (int j=0; j<m; j++) {
        double* dP = colBase(j);
        const double* sP = A.colBase(j);
        const int* idxP  = &idxI[0];
        const int* idxE  = idxP + n;
        while (idxP!=idxE) {
            *dP = sP[*idxP];
            dP ++;
            idxP++;
        }
    }
}


bool XymMatrixVC::isFinite() const
{
    for (int j=0; j<cols(); j++) {
        double *p, *pE;
        int incr;
        loopCol (j, p, pE, incr);
        for (;p!=pE;p+=incr) if (!finite(*p)) return false;
    }
    return true;
}


XymMatrixC XymMatrixVC::transpose () const
{
  XymMatrixC result;
  ::transpose (result, *this);
  return result;  
}


void XymMatrixVC::copyFrom (const XymMatrixVC& m)
{
    ::copyFrom (*this, m);
}

void copyFrom (XymMatrixVC& res, const XymMatrixVC& a)
{
    checkSameFormat (res, a);
    for (int j=0; j<res.cols(); j++) {
        double *rP, *rPE, *aP;
        int rIncr, aIncr;
        res.loopCol (j, rP, rPE, rIncr);
        a.loopCol (j, aP, aIncr);
        for (;rP!=rPE;rP+=rIncr,aP+=aIncr) *rP = *aP;
    }
}


void XymMatrixVC::fill (double alpha, double beta)
{
    for (int j=0; j<cols(); j++) {
        double *rP, *rPE;
        int rIncr;
        loopCol (j, rP, rPE, rIncr);
        for (;rP!=rPE;rP+=rIncr) *rP = alpha;
        (*this)(j,j) = beta;
    }
}


ostream& operator<<(ostream& os, const XymMatrixVC& m)
{
    os << " { " << endl;
    for (int i=0; i<m.rows(); i++) {
        os << " { ";
        for (int j=0; j<m.cols(); j++) os << m(i,j) << " ";
        os << " } " << endl;
    }
    os << " } " << endl;
    return os;
}


istream& operator>>(ifstream& ifs, XymMatrixVC& v) throw (runtime_error)
{
    xymExpect (ifs, "{");
    string buffer;
    ifs >> buffer;
    int rowCtr=0;
    while (buffer!="}") {
        rowCtr++;
        if (rowCtr>v.rows()) throw runtime_error ("XymMatrixVC::operator>> too many rows in matrix");
        if (buffer!="{")  throw runtime_error ("XymMatrixVC::operator>>  missing  { found " + buffer + " instead");
        int colCtr=0;
        ifs >> buffer;
        while (buffer!="}") {
            double d;
            d = atof(buffer.c_str());
            colCtr++;
            if (colCtr>v.cols()) throw runtime_error ("XymMatrixVC::operator>> too many columns in matrix");
            v(rowCtr-1, colCtr-1) = d;
            ifs >> buffer;
        }
        if (colCtr<v.cols()) throw runtime_error ("XymMatrixVC::operator>> too few columns in matrix");
        ifs >> buffer;
    }
    if (rowCtr<v.rows()) throw runtime_error ("XymMatrixVC::operator>> too few rows in matrix");
    return ifs;
}


void XymMatrixVC::print () const
{
    double max=1E-32;
    for (int i=0; i<rows(); i++) for (int j=0; j<cols(); j++)
        if (max<fabs((*this)(i,j))) max=fabs((*this)(i,j));
    max*=10;
    double scale = exp(log(1000.0)*floor(log(max)/log(1000.0)));
    
    printf("%3.1e*{\n", scale);
    for (int i=0; i<rows(); i++) {
        printf ("{ ");
        for (int j=0; j<cols(); j++) {
            printf("%+7.3f", (*this)(i,j)/scale);
            if (j<cols()-1) printf(", ");
        }
        if (i<rows()-1) printf(" },\n");
        else printf(" }\n");
    }
    printf("}\n");
}


void print (const XymMatrixVC& M)
{
    M.print();
}

