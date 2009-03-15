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

#ifndef XYMMATRIXVC_I_H
#define XYMMATRIXVC_I_H

//!\author Udo Frese

/*!\file xymMatrixVC_i.h Contains the inline implementation of class
  \c XymMatrixVC representing a matrix view.
*/

void checkForMultiply (const XymMatrixVC& a, const XymMatrixVC& b)
{
    XYMASSERT (a.cols()==b.rows(), "format mismatch for multiplication");
}

            
void checkForMultiply (const XymMatrixVC& res, const XymMatrixVC& a, const XymMatrixVC& b)
{
    // Format of matrices to be multiplied does not match
    XYMASSERT(a.cols()==b.rows(), "format mismatch for matrix matrix multiplication");
    XYMASSERT(res.rows()==a.rows(), "format mismatch for matrix matrix multiplication");
    XYMASSERT(res.cols()==b.cols(), "format mismatch for matrix matrix multiplication");
}

void checkForMultiply (const XymVectorV& res, const XymMatrixVC& a, const XymVectorV& v)
{
    XYMASSERT(res.size()==a.rows(), "format mismatch for matrix vector multiplication");
    XYMASSERT(a.cols()==v.size(), "format mismatch for matrix vector multiplication");
}

void checkForMultiply (const XymMatrixVC& a, const XymVectorV& v)
{
    XYMASSERT(v.size()==a.rows(), "format mismatch for matrix vector multiplication");
    XYMASSERT(a.cols()==v.size(), "format mismatch for matrix vector multiplication");
}


void checkSameFormat (const XymMatrixVC& a, const XymMatrixVC& b)
{
    XYMASSERT(a.cols()==b.cols(), "matrices do not have same format");
    XYMASSERT(a.rows()==b.rows(), "matrices do not have same format");
}


void checkSquare (const XymMatrixVC& a)
{
    XYMASSERT (a.rows()==a.cols(), "matrix is not square");
}

void XymMatrixVC::store (const Matrix2x2& a, int i, int j)
{
    XYMRANGECHECK (0<=i && (i+1)<_rowDim);
    XYMRANGECHECK (0<=j && (j+1)<_colDim);    
    double* d = &get(i, j);
    int cO = colOfs();
    int rO = rowOfs();
    d[0]         = a[0][0];
    d[cO]        = a[0][1];
    d[rO]        = a[1][0];
    d[rO+cO]     = a[1][1];
}


void XymMatrixVC::store (const Matrix2x3& a, int i, int j)
{
    XYMRANGECHECK (0<=i && (i+1)<_rowDim);
    XYMRANGECHECK (0<=j && (j+2)<_colDim);    
    double* d = &get(i, j);
    int cO = colOfs();
    int rO = rowOfs();
    d[0]         = a[0][0];
    d[cO]        = a[0][1];
    d[2*cO]      =  a[0][2];
    d[rO]        = a[1][0];
    d[rO+cO]     = a[1][1];
    d[rO+2*cO]   = a[1][2];
}


void XymMatrixVC::store (const Matrix3x2& a, int i, int j)
{
    XYMRANGECHECK (0<=i && (i+2)<_rowDim);
    XYMRANGECHECK (0<=j && (j+1)<_colDim);    
    double* d = &get(i, j);
    int cO = colOfs();
    int rO = rowOfs();
    d[0]         = a[0][0];
    d[cO]        = a[0][1];
    d[rO]        = a[1][0];
    d[rO+cO]     = a[1][1];
    d[2*rO]      = a[2][0];
    d[2*rO+cO]   = a[2][1];
}


void XymMatrixVC::store (const Matrix3x3& a, int i, int j)
{
    XYMRANGECHECK (0<=i && (i+2)<_rowDim);
    XYMRANGECHECK (0<=j && (j+2)<_colDim);    
    double* d = &get(i, j);
    int cO = colOfs();
    int rO = rowOfs();
    d[0]         = a[0][0];
    d[cO]        = a[0][1];
    d[2*cO]      = a[0][2];
    d[rO]        = a[1][0];
    d[rO+cO]     = a[1][1];
    d[rO+2*cO]   = a[1][2];
    d[2*rO]      = a[2][0];
    d[2*rO+cO]   = a[2][1];
    d[2*rO+2*cO] = a[2][2];
}
    

void XymMatrixVC::storeCol (const Vector3& a, int i, int j) 
{
    double* d = &get(i, j);
    d[0] = a[0];
    d[1] = a[1];
    d[2] = a[2];
}



void XymMatrixVC::storeCol (const Vector2& a, int i, int j) 
{
    double* d = &get(i, j);
    d[0] = a[0];
    d[1] = a[1];
}


void XymMatrixVC::storeRow (const Vector3& a, int i, int j) 
{
    double* d = &get(i, j);
    d[0] = a[0];
    d[_leadingDimension] = a[1];
    d[2*_leadingDimension] = a[2];
}


void XymMatrixVC::storeRow (const Vector2& a, int i, int j) 
{
    double* d = &get(i, j);
    d[0] = a[0];
    d[_leadingDimension] = a[1];
}


void XymMatrixVC::extract (Matrix2x2& a, int i, int j) const
{
    XYMRANGECHECK (0<=i && (i+1)<_rowDim);
    XYMRANGECHECK (0<=j && (j+1)<_colDim);
    const double* d = &get(i, j);
    int cO = colOfs();
    int rO = rowOfs();
    a[0][0] = d[0];
    a[0][1] = d[cO];
    a[1][0] = d[rO];
    a[1][1] = d[rO+cO];
}


void XymMatrixVC::extract (Matrix2x3& a, int i, int j) const
{
    XYMRANGECHECK (0<=i && (i+1)<_rowDim);
    XYMRANGECHECK (0<=j && (j+2)<_colDim);
    const double* d = &get(i, j);
    int cO = colOfs();
    int rO = rowOfs();
    a[0][0] = d[0];
    a[0][1] = d[cO];
    a[0][2] = d[2*cO];
    a[1][0] = d[rO];
    a[1][1] = d[rO+cO];
    a[1][2] = d[rO+2*cO];
}


void XymMatrixVC::extract (Matrix3x2& a, int i, int j) const
{
    XYMRANGECHECK (0<=i && (i+2)<_rowDim);
    XYMRANGECHECK (0<=j && (j+1)<_colDim);
    const double* d = &get(i, j);
    int cO = colOfs();
    int rO = rowOfs();
    a[0][0] = d[0];
    a[0][1] = d[cO];
    a[1][0] = d[rO];
    a[1][1] = d[rO+cO];
    a[2][0] = d[2*rO];
    a[2][1] = d[2*rO+cO];
}


void XymMatrixVC::extract (Matrix3x3& a, int i, int j) const
{
    XYMRANGECHECK (0<=i && (i+2)<_rowDim);
    XYMRANGECHECK (0<=j && (j+2)<_colDim);
    const double* d = &get(i, j);
    int cO = colOfs();
    int rO = rowOfs();
    a[0][0] = d[0];
    a[0][1] = d[cO];
    a[0][2] = d[2*cO];
    a[1][0] = d[rO];
    a[1][1] = d[rO+cO];
    a[1][2] = d[rO+2*cO];
    a[2][0] = d[2*rO];
    a[2][1] = d[2*rO+cO];
    a[2][2] = d[2*rO+2*cO];
}


void XymMatrixVC::extractCol (Vector3& a, int i, int j) const
{
    const double* d = &get(i, j);
    a[0] = d[0];
    a[1] = d[1];
    a[2] = d[2];
}


void XymMatrixVC::extractRow (Vector3& a, int i, int j) const
{
    const double* d = &get(i, j);
    a[0] = d[0];
    a[1] = d[_leadingDimension];
    a[2] = d[2*_leadingDimension];
}


void XymMatrixVC::assertFormat (int m, int n) const
{
    XYMFORMAT (_d!=NULL && _rowDim==m && _colDim==n);
}


#endif  
