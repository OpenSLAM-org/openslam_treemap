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

#include "xymOperations.h"

#include <math.h>

double dot (const XymVectorV& v1, const XymVectorV& v2)
{
    checkEqualLength (v1, v2);
    double* vd1 = v1.base();
    double* vd2 = v2.base();
    int incr1  = v1.elmOfs();
    int incr2  = v2.elmOfs();
    double sum=0;
    for (int n = v1.size();n>0;n--) {
        sum += *vd1 * *vd2;
        vd1 += incr1;
        vd2 += incr2;
    }
    return sum;
}


double norm (const XymVectorV& v)
{
    return sqrt(dot(v,v));
}

double norminf (const XymVectorV& v){
	double norm=0;
	double *p, *pEnd;
	int inc;
	
	v.loop(p, pEnd, inc);
	for(;p!=pEnd; p+=inc){
		if(norm < fabs(*p)) norm = fabs(*p);
	}
	return norm;
}

double normE (const XymMatrixVC& A)
{
    double sum=0;
    int n = A.rows();
    int m = A.cols();
    for (int j=0; j<m; j++) sum += xymDot (A.colBase(j), A.colBase(j), n);
    return sqrt(sum);
}


void sub (XymVectorV& res, const XymVectorV& v1, const XymVectorV& v2)
{
    checkEqualLength (v1, v2);
    checkEqualLength (res, v1);
    double* vd1 = v1.base();
    double* vd2 = v2.base();
    double* vRes = res.base();
    int incr1  = v1.elmOfs();
    int incr2  = v2.elmOfs();
    int incrR  = res.elmOfs();
    for (int n = v1.size();n>0;n--) {
        *vRes = *vd1 - *vd2;
        vRes += incrR;
        vd1 += incr1;
        vd2 += incr2;
    }
}


void add (XymVectorV& res, const XymVectorV& v1, const XymVectorV& v2)
{
    checkEqualLength (v1, v2);
    checkEqualLength (res, v1);
    double* vd1 = v1.base();
    double* vd2 = v2.base();
    double* vRes = res.base();
    int incr1  = v1.elmOfs();
    int incr2  = v2.elmOfs();
    int incrR  = res.elmOfs();
    for (int n = v1.size();n>0;n--) {
        *vRes = *vd1 + *vd2;
        vRes += incrR;
        vd1 += incr1;
        vd2 += incr2;
    }
}


void add (XymVector& res, const XymVectorV& v1, const XymVectorV& v2)
{
    if (!res.isValid()) res.create (v1.size(), false);
    add ((XymVectorV&) res, v1, v2);
}


void sub (XymVector& res, const XymVectorV& v1, const XymVectorV& v2)
{
    if (!res.isValid()) res.create (v1.size(), false);
    sub ((XymVectorV&) res, v1, v2);
}


void scale (XymVectorV& res, const XymVectorV& u, double scaleFactor)
{
    double *s, *sEnd, *d;
    int sOfs, dOfs;
    checkEqualLength (res, u);
    u.loop (s, sEnd, sOfs);
    res.loop (d, dOfs);
    while (s!=sEnd) {
        *d = scaleFactor * *s;
        d += dOfs;
        s += sOfs;
    }
}


void scale (XymVector& res, const XymVectorV& u, double scaleFactor)
{
    if (!res.isValid()) res.create (u.size(), false);
    scale ((XymVectorV&) res, u, scaleFactor);
}


void scaleAdd (XymVectorV& res, const XymVectorV& u, double scale)
{
    double *s, *sEnd, *d;
    int sOfs, dOfs;
    checkEqualLength (res, u);
    u.loop (s, sEnd, sOfs);
    res.loop (d, dOfs);
    while (s!=sEnd) {
        *d += scale * *s;
        d += dOfs;
        s += sOfs;
    }
}


void linear (XymVectorV& res, const XymVectorV& u, double lambda, const XymVectorV& v, double mu)
{
    double *uP, *uPEnd, *d, *vP;
    int uOfs, dOfs, vOfs;
    checkEqualLength (res, u);
    checkEqualLength (res, v);
    u.loop (uP, uPEnd, uOfs);
    v.loop (vP, vOfs);
    res.loop (d, dOfs);
    while (uP!=uPEnd) {
        *d = lambda * *uP + mu * *vP;
        d += dOfs;
        uP += uOfs;
        vP += vOfs;
    }
}


void linear (XymVector& res, const XymVectorV& u, double lambda, const XymVectorV& v, double mu)
{
    if (!res.isValid()) res.create (u.size());
    linear ((XymVectorV&)res, u, lambda, v, mu);
}


void multiply (XymMatrixVC& res, const XymMatrixVC& a, const XymMatrixVC& b, double factor)
{
    checkForMultiply (res, a, b);
    for (int i=0; i<res.rows(); i++) for (int j=0; j<res.cols(); j++) {
        double *aP, *aPE, *bP;
        int aIncr, bIncr;
        a.loopRow (i, aP, aPE, aIncr);
        b.loopCol (j, bP, bIncr);
        double sum=0;
        for (;aP!=aPE;aP+=aIncr,bP+=bIncr) sum += *aP * *bP;
        res(i,j) = sum*factor;
    }
}


void scale (XymMatrixVC& res, const XymMatrixVC& A, double lambda)
{
    checkSameFormat (res, A);
    int n = res.cols();
    for (int j=0; j<n; j++) {
        double *aP, *aPE, *rP;
        int aOfs, rOfs;
        A.loopCol (j, aP, aPE, aOfs);
        res.loopCol (j, rP, rOfs);
        while (aP!=aPE) {
            *rP = lambda* *aP;
            aP += aOfs;
            rP += rOfs;
        }
    }
}


void scale (XymMatrixC& res,  const XymMatrixVC& A, double lambda)
{
    if (!res.isValid()) res.createLike (A, false);
    scale ((XymMatrixVC&) res,  A, lambda);
}




void multiplyAdd (XymMatrixVC& res, const XymMatrixVC& a, const XymMatrixVC& b, double factor)
{
    checkForMultiply (res, a, b);
    for (int i=0; i<res.rows(); i++) for (int j=0; j<res.cols(); j++) {
        double *aP, *aPE, *bP;
        int aIncr, bIncr;
        a.loopRow (i, aP, aPE, aIncr);
        b.loopCol (j, bP, bIncr);
        double sum=0;
        for (;aP!=aPE;aP+=aIncr,bP+=bIncr) sum += *aP * *bP;
        res(i,j) += sum*factor;
    }
}


void add (XymMatrixVC& res, const XymMatrixVC& a)
{
    checkSameFormat (res, a);
    for (int j=0; j<res.cols(); j++) {
        double *rP, *rPE, *aP;
        int rIncr, aIncr;
        res.loopCol (j, rP, rPE, rIncr);
        a.loopCol (j, aP, aIncr);
        for (;rP!=rPE;rP+=rIncr,aP+=aIncr) *rP += *aP;
    }
}

void multiply (XymMatrixVC &a, double lambda)
{
    for (int j=0; j < a.cols(); j++) {
        double *aP, *aPE;
        int aIncr;
        a.loopCol (j, aP, aPE, aIncr);
        for (;aP!=aPE;aP+=aIncr) *aP *= lambda;
    }
}

void multiply (XymMatrixVC& res, const XymMatrixVC& a, double lambda)
{
    checkSameFormat (res, a);
    for (int j=0; j<res.cols(); j++) {
        double *rP, *rPE, *aP;
        int rIncr, aIncr;
        res.loopCol (j, rP, rPE, rIncr);
        a.loopCol (j, aP, aIncr);
        for (;rP!=rPE;rP+=rIncr,aP+=aIncr) *rP = *aP * lambda;
    }
}


void multiply (XymVectorV& res, const XymMatrixVC& a, const XymVectorV& v, double factor)
{
    checkForMultiply (res, a, v);
    for (int i=0; i<res.size(); i++) {
        double *aP, *aPE, *vP;
        int aIncr, vIncr;
        a.loopRow (i, aP, aPE, aIncr);
        v.loop (vP, vIncr);
        double sum=0;
        for (;aP!=aPE;aP+=aIncr,vP+=vIncr) sum += *aP * *vP;
        res(i) = sum*factor;
    }
}


void multiply (XymVector& res, const XymMatrixVC& a, const XymVectorV& v, double factor)
{
    if (!res.isValid()) res.create(a.rows(), false);
    multiply ((XymVectorV&) res, a, v, factor);
}


void multiplyAdd (XymVectorV& res, const XymMatrixVC& a, const XymVectorV& v, double factor)
{
    checkForMultiply (res, a, v);
    for (int i=0; i<res.size(); i++) {
        double *aP, *aPE, *vP;
        int aIncr, vIncr;
        a.loopRow (i, aP, aPE, aIncr);
        v.loop (vP, vIncr);
        double sum=0;
        for (;aP!=aPE;aP+=aIncr,vP+=vIncr) sum += *aP * *vP;
        res(i) += sum*factor;
    }
}


void multiplyT (XymVectorV& res, const XymMatrixVC& a, const XymVectorV& v, double factor)
{
    XYMASSERT(res.size()==a.cols(),"res and a formats do not match");
    XYMASSERT(a.rows()==v.size(),"a and v formats do not match");
    for (int j=0; j<res.size(); j++) {
        double *aP, *aPE, *vP;
        int aIncr, vIncr;
        a.loopCol (j, aP, aPE, aIncr);
        v.loop (vP, vIncr);
        double sum=0;
        for (;aP!=aPE;aP+=aIncr,vP+=vIncr) sum += *aP * *vP;
        res(j) = sum*factor;
    }
}


void multiplyT (XymVector& res, const XymMatrixVC& a, const XymVectorV& v, double factor)
{
    if (!res.isValid()) res.create(a.cols(), false);
    multiplyT ((XymVectorV&) res, a, v, factor);
}


void multiplyTAdd (XymVectorV& res, const XymMatrixVC& a, const XymVectorV& v, double factor)
{
    XYMASSERT(res.size()==a.cols(),"res and a formats do not match");
    XYMASSERT(a.rows()==v.size(),"a and v formats do not match");
    for (int j=0; j<res.size(); j++) {
        double *aP, *aPE, *vP;
        int aIncr, vIncr;
        a.loopCol (j, aP, aPE, aIncr);
        v.loop (vP, vIncr);
        double sum=0;
        for (;aP!=aPE;aP+=aIncr,vP+=vIncr) sum += *aP * *vP;
        res(j) += sum*factor;
    }
}


void add (XymMatrixVC& res, const XymMatrixVC& a, const XymMatrixVC& b)
{
    checkSameFormat (res, a);
    checkSameFormat (res, b);
    for (int j=0; j<res.cols(); j++) {
        double *rP, *rPE, *aP, *bP;
        int rIncr, aIncr, bIncr;
        res.loopCol (j, rP, rPE, rIncr);
        a.loopCol (j, aP, aIncr);
        b.loopCol (j, bP, bIncr);
        for (;rP!=rPE;rP+=rIncr,aP+=aIncr,bP+=bIncr) *rP = *aP + *bP;
    }
}


void add (XymMatrixC& res, const XymMatrixVC& a, const XymMatrixVC& b)
{
    if (!res.isValid()) res.createLike (a, false);
    add ((XymMatrixVC&) res, a, b);
}


void sub (XymMatrixVC& res, const XymMatrixVC& a)
{
    checkSameFormat (res, a);
    for (int j=0; j<res.cols(); j++) {
        double *rP, *rPE, *aP;
        int rIncr, aIncr;
        res.loopCol (j, rP, rPE, rIncr);
        a.loopCol (j, aP, aIncr);
        for (;rP!=rPE;rP+=rIncr,aP+=aIncr) *rP -= *aP;
    }
}

void sub (XymMatrixVC& res, const XymMatrixVC& a, const XymMatrixVC& b)
{
    checkSameFormat (res, a);
    checkSameFormat (res, b);
    for (int j=0; j<res.cols(); j++) {
        double *rP, *rPE, *aP, *bP;
        int rIncr, aIncr, bIncr;
        res.loopCol (j, rP, rPE, rIncr);
        a.loopCol (j, aP, aIncr);
        b.loopCol (j, bP, bIncr);
        for (;rP!=rPE;rP+=rIncr,aP+=aIncr,bP+=bIncr) *rP = *aP - *bP;
    }
}


void sub (XymMatrixC& res, const XymMatrixVC& a, const XymMatrixVC& b)
{
    if (!res.isValid()) res.createLike (a, false);
    sub ((XymMatrixVC&) res, a, b);
}
    

void addvvt (XymMatrixVC& res, const XymVectorV& v, double scale)
{
    checkSquare (res);
    checkForMultiply (res, v);
    int n = v.size();
// res(i,j) += scale*v(i)*v(j)
    for (int i=0; i<n; i++) {
        double f = scale*v(i);
        double* resP, *resPEnd, *vP;
        int resIncr, vIncr;
        res.loopCol (i, resP, resPEnd, resIncr);
        v.loop (vP, vIncr);
        for (;resP!=resPEnd;resP+=resIncr,vP+=vIncr) *resP += f* *vP;
    }
}


void adduvt (XymMatrixVC& res, const XymVectorV& u, const XymVectorV& v, double scale)
{
    XYMFORMAT (res.rows()==u.size());
    XYMFORMAT (res.cols()==v.size());
    int n = v.size();
// res(i,j) += scale*u(i)*v(j)
    for (int j=0; j<n; j++) {
        double f = scale*v(j);
        double* resP, *resPEnd, *uP;
        int resIncr, uIncr;
        res.loopCol (j, resP, resPEnd, resIncr);
        u.loop (uP, uIncr);
        for (;resP!=resPEnd;resP+=resIncr,uP+=uIncr) *resP += f* *uP;
    }
}


XymMatrixC outer (const XymVectorV& u, const XymVectorV& v, double scale)
{
  XymMatrixC result (u.size(), v.size());
  adduvt (result, u, v, scale);
  return result;
}


void addJCJt (XymMatrixVC& res, const XymMatrixVC& J, const XymMatrixVC& C, double factor)
{
    checkSquare (res);
    checkSquare (C);
    checkForMultiply (res, J);
    checkForMultiply (J, C);
    int n = J.rows();
    int m = J.cols();
    // res(i,j) += J(i,k)*C(k,l)*J(j,l)
    
    // CODE using '(i,j)' operator
    // for (int j=0;j<n;j++) for (int k=0; k<m; k++) {
    //    double CklJjl = 0;
    //    for (int l=0; l<m; l++) CklJjl   += C(k,l)*J(j,l);
    //    for (int i=0; i<n; i++) res(i,j) += J(i,k)*CklJjl;
    // }
    for (int j=0;j<n;j++) for (int k=0;k<m;k++) {
        double CklJjl = 0;
        double* JP, *JPEnd, *CP;
        int jIncr, cIncr;
        J.loopRow (j, JP, JPEnd, jIncr);
        C.loopRow (k, CP, cIncr);
        for (;JP!=JPEnd;JP+=jIncr,CP+=cIncr) CklJjl += *JP * *CP;
        
        CklJjl *= factor;
        
        double *resP;
        int resIncr;
        J.loopCol (k, JP, jIncr);
        JPEnd = JP+j+1;
        res.loopCol (j, resP, resIncr);
        for (;JP!=JPEnd;JP+=jIncr,resP+=resIncr) *resP += *JP * CklJjl;
    }
    copyUpperToLower(res);
}


void linear (XymMatrixVC& res, const XymMatrixVC& a, double lambda)
{
    checkSameFormat (res, a);
    for (int j=0; j<res.cols(); j++) {
        double *rP, *rPE, *aP;
        int rIncr, aIncr;
        res.loopCol (j, rP, rPE, rIncr);
        a.loopCol (j, aP, aIncr);
        for (;rP!=rPE;rP+=rIncr,aP+=aIncr) *rP += *aP*lambda;
    }
}


double symVMV (const XymVectorV& u, const XymMatrixVC& M)
{
    int n = M.rows();
    XYMFORMAT (u.size()==n);
    checkSquare (M);
    int m = n-1;
    double sum=u[m] * M(m,m) * u[m];
    double* Mp = M.base();
    double* uP = u.base();
    XYMASSERT (u.elmOfs()==1, "stride of vector u must be 1");
    for (int j=0; j<m; j++) { // Exploiting symmetry
        sum += u[j]*((*uP) * (*Mp) + 2*xymDot (Mp+1, uP+1, m-j));
        Mp  += M.colOfs()+M.rowOfs();
        uP  ++;
    }
    return sum;
}


double VMV (const XymVectorV& u, const XymMatrixVC& M, const XymVectorV& v)
{
    int n = M.rows();
    XYMFORMAT (u.size()==n);
    XYMFORMAT (v.size()==n);
    checkSquare (M);
    double sum=0;
    double* Mp = M.base();
    double* uP = u.base();
    XYMASSERT (u.elmOfs()==1, "stride of vector u must be 1");
    for (int j=0; j<n; j++) { 
        sum += v[j]*xymDot (Mp, uP, n);
        Mp  += M.colOfs();
    }
    return sum;
}


double trace (const XymMatrixVC& A)
{
    checkSquare(A);
    double* p=A.base();
    int ofs = A.rowOfs()+A.colOfs();
    double sum=0;
    double* pEnd = p + A.rows()*ofs;
    for (;p!=pEnd;p+=ofs) sum += *p;
    return sum;
}


void multiply (XymMatrixC& res, const XymMatrixVC& a, const XymMatrixVC& b, double factor)
{
    if (!res.isValid()) res.create (a.rows(), b.cols(), false);
    multiply ((XymMatrixVC&) res, a, b, factor);
}


void transpose (XymMatrixC& result, const XymMatrixVC& A)
{
    int m=A.rows();
    int n=A.cols();
    result.create (n,m); // transposed format
    transpose( (XymMatrixVC&) result, A);
}	

void transpose (XymMatrixVC& result, const XymMatrixVC& A)
{
    int m=A.rows();
    int n=A.cols();
    XYMASSERT(result.rows()==n && result.cols()==m,
    	"Dimension mismatch");
    for (int j=0;j<n;j++) {
        // result(j,i) = A(i,j)
        double *pS, *pE, *pD;
        int ofsS, ofsD;
        A.loopCol (j, pS, pE, ofsS);
        result.loopRow (j, pD, ofsD);
        while (pS!=pE) {
            *pD = *pS;
            pD += ofsD;
            pS += ofsS;
        }
    }
}


void transpose (XymMatrixVC& A)
{
    checkSquare(A);
    int n=A.cols();
    for (int i=0;i<n;i++) {
        // result(j,i) = A(i,j)
        double *pS, *pE, *pD;
        int ofsS, ofsD;
        A.loopCol (i, pS, ofsS);
        pE = pS + ofsS*i;
        A.loopRow (i, pD, ofsD);
        while (pS!=pE) {
            // swap *pS and *pD
            double buf = *pD;
            *pD = *pS;
            *pS = buf;

            pD += ofsD;
            pS += ofsS;
        }
    }
}


void symmetrize (XymMatrixVC& A)
{
    checkSquare(A);
    int n=A.cols();
    for (int i=0;i<n;i++) {
        // result(j,i) = A(i,j)
        double *pS, *pE, *pD;
        int ofsS, ofsD;
        A.loopCol (i, pS, ofsS);
        pE = pS + ofsS*i;
        A.loopRow (i, pD, ofsD);
        while (pS!=pE) {
            *pS = *pD = (*pS + *pD)/2;

            pD += ofsD;
            pS += ofsS;
        }
    }
}

void copyUpperToLower (XymMatrixVC& A)
{
    checkSquare(A);
    int n=A.cols();
    for (int i=0;i<n;i++) {
        // result(j,i) = A(i,j)
        double *pS, *pE, *pD;
        int ofsS, ofsD;
        A.loopCol (i, pS, ofsS);
        pE = pS + ofsS*i;
        A.loopRow (i, pD, ofsD);
        while (pS!=pE) {
            *pD = *pS;

            pD += ofsD;
            pS += ofsS;
        }
    }
}

void zeroLower (XymMatrixVC& A)
{
  if (A.rows()==0) return;  
  int n = A.cols();
  int m = A.rows();  
  for (int i=0;i<n;i++) {
    double *pE, *pD;
    int ofsD;
    A.loopCol (i, pD, pE, ofsD);
    if (i+1<m) {      
      pD  +=  ofsD*(i+1);
      while (pD!=pE) {
        *pD = 0;
        pD += ofsD;
      }
    }    
  }
}


void zeroUpper (XymMatrixVC& A)
{
  if (A.rows()==0) return;  
  int n = A.cols();
  int m = A.rows();  
  for (int i=0;i<n;i++) {
    double *pE, *pD;
    int ofsD;
    A.loopCol (i, pD, ofsD);
    if (i<m) pE  = pD + i*ofsD;
    else pE = pD + m*ofsD;    
    while (pD!=pE) {
      *pD = 0;      
      pD += ofsD;
    }
  }
}


double symmetry (const XymMatrixVC& A)
{
    checkSquare (A);
    int n = A.rows();
    double maxDist = 0;
    for (int i=0; i<n; i++) for (int j=i+1; j<n; j++) {
        double dist = fabs(A(i,j)-A(j,i));
        if (dist>maxDist) maxDist = dist;
    }
    return maxDist;
}


void copyLowerToUpper (XymMatrixVC& A)
{
    checkSquare(A);
    int n=A.cols();
    for (int i=0;i<n;i++) {
        // result(j,i) = A(i,j)
        double *pS, *pE, *pD;
        int ofsS, ofsD;
        A.loopCol (i, pS, ofsS);
        pE = pS + ofsS*i;
        A.loopRow (i, pD, ofsD);
        while (pS!=pE) {
            *pS = *pD;

            pD += ofsD;
            pS += ofsS;
        }
    }
}


void addIndexed (XymVectorV& u, const XymVectorV& v, const vector<int>& idx, double factor)
{
    XYMFORMAT (v.size()==(int) idx.size());
    int n = idx.size();
#ifndef NDEBUG
    int m = u.size();
#endif
    const int* idxP = &idx[0];
    const int* idxPEnd = idxP + n;
    double *srcP;
    int srcInc;
    v.loop (srcP, srcInc);
    while (idxP!=idxPEnd) {
        int idx = *idxP;
        XYMASSERT (idx<m, "index must be insert vector u");
        if (idx>=0) u[idx] += factor * *srcP;
        srcP += srcInc;
        idxP ++;
    }
}


void addIndexed (XymMatrixVC& A, const XymMatrixVC& Q, const vector<int>& idxI, const vector<int>& idxJ, double factor)
{
    XYMFORMAT (Q.rows()==(int) idxI.size());
    XYMFORMAT (Q.cols()==(int) idxJ.size());
 
    int m=Q.rows();
    int n=Q.cols();
    for (int i=0; i<m; i++) {
      XYMASSERT(-1<=idxI[i] && idxI[i]<A.rows(), "index 'idxI' must be inside matrix A");
    }
    for (int i=0; i<n; i++) {
      XYMASSERT(-1<=idxJ[i] && idxJ[i]<A.cols(), "index 'idxJ' must be inside matrix A");
    }    
    for (int j=0; j<n; j++) {
        int jA = idxJ[j];
        if (jA>=0) {
            const int* idxP = &idxI[0];
            const int* idxPEnd = idxP + m;
            double* aCol = A.colBase(jA);
            double *srcP;
            int srcInc;
            Q.loopCol (j, srcP, srcInc);
            while (idxP!=idxPEnd) {
                int id = *idxP;
                if (id>=0) aCol[id] += factor * *srcP;
                srcP += srcInc;
                idxP ++;
            }
        }
    }
}


void addIndexedCol (XymMatrixVC& A, const XymMatrixVC& Q, const vector<int>& idxJ, double factor)
{
    XYMFORMAT (Q.cols()==(int) idxJ.size());
 
    int n=Q.cols();
    for (int i=0; i<n; i++) {
      XYMASSERT (-1<=idxJ[i] && idxJ[i]<A.cols(), "index 'idxI' must be inside matrix A");
      }    
    for (int j=0; j<n; j++) {
        int jA = idxJ[j];
        if (jA>=0) {
            double *srcP, *srcPEnd;
            double *dstP;
            int srcInc, dstInc;
            Q.loopCol (j, srcP, srcPEnd, srcInc);
            A.loopCol (jA, dstP, dstInc);
            while (srcP!=srcPEnd) {
                *dstP += factor * *srcP;
                srcP += srcInc;
                dstP += dstInc;
            }
        }
    }
}


void rotate2D (XymMatrixVC& result, const XymMatrixVC& A, double angle)
{
    // compute R*A*R^T
    XYMFORMAT (result.rows()==A.rows() && result.cols()==A.cols());
    XYMFORMAT (A.rows()%2==0 && A.cols()%2==0);
    double c = cos(angle), s = sin(angle), c2=c*c, cs=c*s, s2=s*s;
    for (int j=0; j<A.cols(); j+= 2) {
        double* srcP, * srcPEnd, *dstP;
        int inc, srcColOfs=A.colOfs(), dstColOfs=result.colOfs();
        A.loopCol (j, srcP, srcPEnd, inc);
        XYMASSERT (inc==1, "stride of columns of A must be 1");
        result.loopCol (j, dstP, inc);
        XYMASSERT (inc==1, "stride of columns of result must be 1");
        for (;srcP!=srcPEnd;srcP+=2,dstP+=2) {
            // Now process a 2*2 block
            double a00, a01, a10, a11, ad, ao;
            a00 = srcP[0];
            a01 = srcP[srcColOfs];
            a10 = srcP[1];
            a11 = srcP[1+srcColOfs];
            ao  = a01+a10;
            ad  = a00-a11;

            dstP[0]           = c2*a00 - cs*ao + s2*a11;
            dstP[dstColOfs]   = c2*a01 + cs*ad - s2*a10;
            dstP[1]           = c2*a10 + cs*ad - s2*a01;
            dstP[1+dstColOfs] = c2*a11 + cs*ao + s2*a00;
        }
    }
}


void rotate2D (XymMatrixC& result, const XymMatrixVC& A, double angle)
{
    if (!result.isValid()) result.createLike (A, false);
    rotate2D ((XymMatrixVC&) result, A, angle);
}


void rotate2DRows (XymMatrixVC& result, const XymMatrixVC& A, double angle)
{
    // compute R*A
    XYMFORMAT (result.rows()==A.rows() && result.cols()==A.cols());
    XYMFORMAT (A.rows()%2==0);
    double c = cos(angle), s = sin(angle);
    for (int j=0; j<A.cols(); j+= 1) {
        double* srcP, * srcPEnd, *dstP;
        int inc;
        A.loopCol (j, srcP, srcPEnd, inc);
        XYMASSERT (inc==1, "stride of columns of A must be 1");
        result.loopCol (j, dstP, inc);
        XYMASSERT (inc==1, "stride of columns of result must be 1");
        for (;srcP!=srcPEnd;srcP+=2,dstP+=2) {
            // Now process a 2*2 block
            double a00, a10;
            a00 = srcP[0];
            a10 = srcP[1];

            dstP[0] = c*a00 - s*a10;
            dstP[1] = c*a10 + s*a00;
        }
    }
}


void rotate2DRows (XymMatrixC& result, const XymMatrixVC& A, double angle)
{
    if (!result.isValid()) result.createLike (A, false);
    rotate2DRows ((XymMatrixVC&) result, A, angle);
}


void rotate2DCols (XymMatrixVC& result, const XymMatrixVC& A, double angle)
{
    // compute A*R^T
    XYMFORMAT (result.rows()==A.rows() && result.cols()==A.cols());
    XYMFORMAT (A.cols()%2==0);
    double c = cos(angle), s = sin(angle);
    for (int j=0; j<A.cols(); j+= 2) {
        double* srcP, * srcPEnd, *dstP;
        int inc, srcColOfs=A.colOfs(), dstColOfs=result.colOfs();
        A.loopCol (j, srcP, srcPEnd, inc);
        XYMASSERT (inc==1, "stride of columns of A must be 1");
        result.loopCol (j, dstP, inc);
        XYMASSERT (inc==1, "stride of columns of result must be 1");
        for (;srcP!=srcPEnd;srcP++,dstP++) {
            // Now process a 2*2 block
            double a00, a01;
            a00 = srcP[0];
            a01 = srcP[srcColOfs];

            dstP[0]           = c*a00 - s*a01;
            dstP[dstColOfs]   = c*a01 + s*a00;
        }
    }
}


void rotate2DCols (XymMatrixC& result, const XymMatrixVC& A, double angle)
{
    if (!result.isValid()) result.createLike (A, false);
    rotate2DCols ((XymMatrixVC&) result, A, angle);
}


void transform2D (XymVectorV& result, const XymVectorV& v, double dx, double dy, double angle)
{
    XYMFORMAT (result.size()==v.size());
    XYMFORMAT (v.size()%2==0);
    double c = cos(angle), s = sin(angle);
    double* srcP, * srcPEnd, *dstP;
    int srcInc, dstInc;
    v.loop (srcP, srcPEnd, srcInc);
    result.loop (dstP, dstInc);
    while (srcP!=srcPEnd) {
        double a0, a1;
        a0 = srcP[0];
        a1 = srcP[srcInc];
        
        dstP[0]      = a0*c - a1*s + dx;
        dstP[dstInc] = a0*s + a1*c + dy;
        dstP += 2*dstInc;
        srcP += 2*srcInc;
    }
}


void transform2D (XymVector& result, const XymVectorV& v, 
                  double dx, double dy, double angle)
{
    if (!result.isValid()) result.create (v.size(), false);
    transform2D ((XymVectorV&) result, v, dx, dy, angle);
}


void cholesky (XymMatrixC& U, const XymMatrixVC& A, double minDiagonal, bool setToMin)
{
    checkSquare(A);
    int n = A.rows();
    U.createLike (A, false);
    for (int i=0; i<n; i++) {
        for (int j=0; j<i; j++) U(i,j)=0;

        double sum, sqrtSum, sqrtSumInv;
        sum = A(i,i) - xymDot(U.colBase(i), U.colBase(i), i);
        if (sum<=minDiagonal) { // Error: Cholesky failed
            if (setToMin && minDiagonal>0) sum = minDiagonal;
            else {
                U.clear();
                return;
            }
        }
        sqrtSum = sqrt(sum);
        U(i,i)=sqrtSum;

        sqrtSumInv = 1/sqrtSum;
        for (int j=i+1; j<n; j++) {
            // Sum U(k,i)*U(k,j)
            double sum = A(i,j) - xymDot (U.colBase(i), U.colBase(j), i);
            U(i,j) = sqrtSumInv*sum;
        }
    }
}


void invertUT (XymMatrixC& result, const XymMatrixVC& U, bool refine)
{
    XYMASSERT (!refine, "refine not implemented");
    int n = U.rows();
    result.createLike (U, false);
    for (int i=0; i<n; i++) {
        result(i,i) = 1/U(i,i);
        for (int j=i+1; j<n; j++) {
            // Sum U(k,j)*result(k,i)
            double sum  = xymDot (U.colBase(j)+i, result.colBase(i)+i, j-i);
            result(j,i) = -sum/U(j,j);
        }
        for (int j=0;j<i;j++) result(j,i)=0;
    }
}


void multiplyLTL (XymMatrixC& result, const XymMatrixVC& L)
{
    checkSquare(L);
    int n = L.rows();
    result.create (n, n);
    for (int i=0; i<n; i++) {
        for (int j=0; j<=i; j++) {
            // result(i,j)=L(k,i)*L(k,j)
            double res = xymDot (L.colBase(i)+i, L.colBase(j)+i, n-i);
            result(i,j) = res;
            result(j,i) = res;
        }
    }
}

void multiplyJTJ (XymMatrixVC& result, const XymMatrixVC& J){
	checkSquare(result);
	checkForMultiply(J, result);
	
	for(int i=0; i<result.cols(); i++){
		for(int j=0; j<=i; j++){
			double *jp, *jtp, *jpE;
			int step;
			result(i,j) = 0;
			J.loopCol(i, jp, jpE, step);
			J.loopCol(j, jtp, step);
			for(;jp != jpE; jp+=step, jtp+=step){
				result(i,j) += *jp * *jtp;
			}
			result(j,i) = result(i,j);
		}
	}
}

void choleskyInvert (XymMatrixC& result, const XymMatrixVC& A, const XymMatrixVC& U, bool refine)
{
    if (U.isValid()) {
        XymMatrixC UI;
        invertUT (UI, U);
        multiplyLTL (result, UI);
    }
    else {
        cholesky (result, A);
        if (!result.isValid()) return;
        XymMatrixC UI;
        invertUT (UI, result);
        multiplyLTL (result, UI);
    }    
}


void inverseCholeskyFactor (XymMatrixC& result, const XymMatrixVC& A, bool refine)
{
  XymMatrixC U;
  cholesky (U, A);
  invertUT (result, U);  
}



void choleskySolve (XymVector& x, const XymMatrixVC& A, const XymVectorV& b, 
                    const XymMatrixVC& U, bool refine)
{
    if (!U.isValid()) {
        XymMatrixC UChol;
        cholesky (UChol, A);
        if (UChol.isValid()) choleskySolve (x, A, b, UChol);
        else x.clear();
        return;
    }
    
    if (!x.isValid()) x.create (A.cols());

    checkSquare (A);
    checkForMultiply (b, A, x);
    checkSquare (U);
    XYMFORMAT (U.rows()==A.rows());

    // Solve 'Uty=b' storing y in 'x'
    int n = A.cols();
    for (int i=0; i<n; i++) {
        // (b[i] - sum U(k,i)*x(k), k=0..i-1) / U(i,i)
        x[i] = (b[i] - xymDot (U.colBase(i), x.base(), i)) / U(i,i);
    }

    // Solve 'Ux=y'
    for (int i=n-1; i>=0; i--) {
        // (x[i] - sum U(i,k)*x(k), k=i+1..n-1)/U(i,i)
        x[i] = (x[i] - xymDot (U.rowBase(i)+(i+1)*U.colOfs(), U.colOfs(), x.base()+i+1, 1, n-i-1)) / U(i,i);
    }
}



void shermanMorisonUpdate (XymMatrixVC& AI, XymVectorV& u, double factor)
{
    checkSquare(AI);
    int n = AI.rows();
    XYMFORMAT (n==u.size());
    XymVector aIu (n);
    multiply (aIu, AI, u);
    double lambda = factor * dot(aIu, u);
    addvvt (AI, aIu, factor/(1+lambda));
}


void choleskyStep (XymMatrixVC& A, int col, XymVectorV& u, double minDiagonal)
{
    int n = A.cols();
    double diag = A(col,col);
    if (diag>=minDiagonal) {
        scale (u, A.col(col), 1/sqrt(diag));
        addvvt (A, u, -1);
    }
    else {
        for (int i=0; i<n; i++) {
#ifndef NDEBUG
            double uI=A(col, i);
#endif
            double d=A(i,i);
            if (d<minDiagonal) d=minDiagonal;
            // Assert, that the off-diagonal entries are small
            // too. This bound would follow from the S.P.D condition,
            // if the diagonal was 'minDiagonal.'
            XYMASSERT (fabs(uI*uI)<=minDiagonal*d, "matrix A is not SPD");
            u(i)=0;
        }
    }
    for (int i=0; i<n; i++) A(i,col) = A(col, i) = 0;
}


double xymDot (double* a, double* b, int n)
{
    double* end = a+n;
    double sum = 0;
    while (a!=end) {
        sum += *a * *b;
        a++;
        b++;
    }
    return sum;
}

double xymDot (double* a, int aOfs, double* b, int bOfs, int n)
{
    double* end = a+n*aOfs;
    double sum = 0;
    while (a!=end) {
        sum += *a * *b;
        a += aOfs;
        b += bOfs;
    }
    return sum;
}

