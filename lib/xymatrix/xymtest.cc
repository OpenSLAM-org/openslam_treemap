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

#include <xymatrix/xymMatrixC.h>
#include <xymatrix/xymVector.h>
#include <xymatrix/xymBLAS.h>
#include <math.h>


void invert (XymMatrixVC& aInv, const XymMatrixVC& a)
{
    checkSquare (aInv);
    checkSquare (a);
    if (a.cols()==1) {
        aInv(0,0) = 1/a(0,0);
    }
    else {
        // Submatrix formula (Numerical Recipes Chapter 2.7)
        int n = a.cols()/2;
        int m = a.cols()-n;
        // Submatrices of a
        XymMatrixVC P, Q, R, S;
        P = a.submatrix (0, 0, m, m);
        Q = a.submatrix (0, m, m, n);
        R = a.submatrix (m, 0, n, m);
        S = a.submatrix (m, m, n, n);
        // Submatrices of aInv
        XymMatrixVC PInv, QInv, RInv, SInv;
        PInv = aInv.submatrix (0, 0, m, m);
        QInv = aInv.submatrix (0, m, m, n);
        RInv = aInv.submatrix (m, 0, n, m);
        SInv = aInv.submatrix (m, m, n, n);
        // Subterms
        XymMatrixC Si (n, n), QSi (m,n), SiR(n,m), QSiR (m,m);
        invert (Si, S);
        multiply (SiR, Si, R);
        multiply (QSi, Q, Si);
        multiply (QSiR, QSi, R);
        XymMatrixC schurC(m, m);
        sub (schurC, P, QSiR);
        invert (PInv, schurC);
        
        multiply (QInv, PInv, QSi);
        multiply (QInv, QInv, -1);

        multiply (RInv, SiR, PInv);
        multiply (RInv, RInv, -1);
        
        multiply (SInv, RInv, QSi);
        sub (SInv, Si, SInv);

        // Check for correctness
        XymMatrixC id(a,0);
        multiply (id, a, aInv);
//        cout << " M " << id << endl;
    }
}


// Checks, whether M is SPD
void checkForSPD (const XymMatrixVC& M)
{
    checkSquare(M);
    int n = M.cols();
    XymMatrixC L(M); // copy
    XymVector p(n);
    for (int i=0;i<n;i++) for(int j=i;j<n;j++) {
        double sum=L(i,j);
        for(int k=i-1;k>=0;k--) sum -= L(i,k)*L(j,k);
        if (i==j) {
            assert (sum>0);
            p(i) = sqrt(sum);
        }
        else L(j,i)=sum/p(i);
    }
}


// Checks, whether a random submatrix of M is SPD
void checkSubMatrix (const XymMatrixVC& M)
{
    int n = M.cols(), i1, i2;
    if (n==0) return;
    i1=rand() % n;
    i2=rand() % n;
    if (i1>i2) {
        int park = i1;
        i1 = i2;
        i2 = park;
    }
    checkForSPD (M.submatrix(i1, i1, i2-i1+1, i2-i1+1));
}


// Checks, whether M and MInv are inverse to each other.
void checkInverse (const XymMatrixVC& M, const XymMatrixVC& MInv)
{
    XymMatrixC id;
    multiply (id, MInv, M);
    XymVectorV diag = id.diag();
    for (int i=0;i<diag.size();i++) diag(i)-=1;
    double sum=0;
    for (int i=0;i<id.rows();i++) for(int j=0;j<id.cols();j++) sum += fabs(id(i,j));
    assert (sum<1E-4); // id is identitymatrix
}

void checkAccess (XymMatrixVC& M)
{
    if (M.rows()==0 || M.cols()==0) return;
    // Checks whether the different access operators work properly
    for (int k=0; k<100; k++) {
        int i = rand()%M.rows();
        int j = rand()%M.cols();
        double a = M(i,j);
        // Element access
        assert (a==M[i][j]);
        assert (a==M.get(i,j));
        assert (a==M.getNC(i,j));

        // row/col/diag access
        assert (a==M.row(i)[j]);
        assert (a==M.col(j)[i]);
        assert (a==M.diag(i-j)[min(i,j)]);

        // Direct Access
        assert (a==M.base()[i*M.rowOfs()+j*M.colOfs()]);
        assert (a==M.rowBase (i)[j*M.colOfs()]);
        assert (a==M.colBase (j)[i*M.rowOfs()]);
        double*d;
        int dOfs;
        M.loopRow (i, d, dOfs);
        assert (a==d[j*dOfs]);
        M.loopCol (j, d, dOfs);
        assert (a==d[i*dOfs]);
        
        int top  = rand()%(i+1);
        int bottom = rand()%(M.rows()-i)+i;
        int left = rand()%(j+1);
        int right= rand()%(M.cols()-j)+j;
        XymMatrixVC SM = M.submatrix (top, left, bottom-top+1, right-left+1);
        assert (a==SM(i-top, j-left));
    }
}


void checkBlockMatrix (XymMatrixVC& M)
{
    // Checks the submatrix access, Block matrix composition and 'store'/'extract' functions.
}



void rowColumnTest()
{
    // row / column insertion / deletion example
    // We take a SPD matrix M an repeatedly insert and delete columns / rows
    // We maintain M^-1 by updating
    XymMatrixC M, MInv;
    int nMax=20;
    M.reserve (nMax, nMax);
    MInv.reserve (nMax, nMax);
    srand(1);
    for (int ctr=0; ctr<1000; ctr++) {
        int n=M.rows();
        checkSubMatrix (M);
        checkSubMatrix (MInv);
        checkInverse (M, MInv);
        checkAccess (M);
        int rC = rand()%10;
        if (n<nMax && rC<=1) { // 20% chance for inserting a row / column
            int i = rand()%(n+1);
            cout << "Inserting row/column " << i << endl;
            M.insertRowAndColumn (i, i, 1, 1);
            M(i,i)=1.0;
            MInv.insertRowAndColumn (i, i, 1, 1);
            MInv(i,i)=1.0;
        }
        else if (rC<=2 && n>0) { // 10% for deleting a row/column
            int i = rand()%n;
            cout << "Deleting row/column " << i << endl;
//            cout << M << MInv << endl;
            // P~^-1 = P - Q*S^-1*R, with M being A^-1 and MInv being A
            XymVector q(MInv.row(i));
            double si = -1/MInv(i,i);
            addvvt (MInv, q, si);
//            cout << "before delete " << M << MInv << endl;
            MInv.deleteRowAndColumn (i, i, 1, 1);
            M.deleteRowAndColumn (i, i, 1, 1);
//            cout << "after delete " << M << MInv << endl;
        }
        else if (n>=1) { // 70% for adding a positive definite submatrix
            int i1 = rand()%n, i2=rand()%n;
            cout << "Adding a 2*2 block with row/column " << i1 << " and " << i2 << endl;
            // sigma*(+e_i1-e_i2)^2
            //           double sigma = rand()%100/100.0; // 0..1
            double sigma = 1;
            XymVector u (n, true);
            u(i1) = sigma;
            u(i2) = -sigma;
            addvvt (M, u);
            // Now adapt the matrix via sherman-morrison formula
            XymVector MIu;
            multiply (MIu, MInv, u);
            double scale = -1/(1+dot(u, MIu));
            addvvt (MInv, MIu, scale);
        }
//        cout << endl << "M: " << endl << M << endl;
    };
} 

void polynomfitTest()
{
    // General Matrix handling : Polynom fit example
    double x[] = {1, 2, 3, 4, 5, 6};
    // 6*4 Vandermonde matrix for x_i = i. EVIL conditioned, be aware
    double d[] = {1, 1, 1, 1, 
                  1, 2, 4, 8,  
                  1, 3, 9, 27,  
                  1, 4, 16, 64,  
                  1, 5, 25, 125,  
                  1, 6, 36, 216}; 
    XymMatrixC v(d, 6, 4); // Takes a row major memory argument and copies

    cout << "Matrix V" << endl;
    cout << v;


    // Coefficient of a 3 degree polynom P
    double coef[] = {1, 0.3, 0.05, 0.007};
    
    // Compute b(i) = P(x[i])
    XymVector b (v.rows());
    for (int i=0; i<b.size(); i++) {
        double sum=0;
        double xPi = 1, xx=x[i];
        for (int k=0; k<4; k++) {
            sum += coef[k]*xPi;
            xPi *= xx;
        }
        b(i)  = sum;
    }

    XymVector b2;
    multiply (b2, v, XymVectorV(coef, 4));
    cout << "b=" << b << " = " << b2 << endl;
    

    // We want to make an least square fit of P(x[i]) to a 4th order polynomial.
    // For that we should compute VtVx=Vtb. Let's compute 'VtV'
    XymMatrixC VtV (v.cols(), v.cols());

    for (int i=0; i<v.cols(); i++) for (int j=0; j<v.cols(); j++) {
        double sum = 0;
        for (int k=0; k<v.rows(); k++) sum += v(k,i)*v(k,j);
        VtV(i,j)=sum;
    }


    
    // Now we compute (VtV)^-1 by partitioning. This is a silly way to do it, but
    // allows to check the submatrix functions
//    double d3[] = {1,1,0,1};
//    VtV.clear();
//    copyFrom(VtV, XymMatrixC (d3, 2, 2));
    XymMatrixC vtVInv (VtV, 0);
    invert (vtVInv, VtV);

    cout << " Matrices " << endl << VtV << endl << endl << vtVInv << endl;
    
    checkInverse (VtV, vtVInv);
    
     // Now compute result = (VtV)^-1Vtb
    XymVector result1, result;
    multiplyT (result1, v, b);
    multiply  (result, vtVInv, result1);

    // As we have fitted a fourth order polynomial 'coef[]' to a fourth
    // order polynomial, the originial polynomial should result;
    cout << "Input vector" << XymVectorV(coef, 4) << endl;
    cout << "Resulting vector " << result << endl;

    double sum=0;
    for (int i=0;i<4;i++) sum += fabs(coef[i]-result(i));
    assert(sum<1E-6);
}


void blasTest() 
{
    double v1B[3] = {1,2,3};
    double v2B[3] = {3,2,1};
    XymVector v1(v1B, 3);
    XymVector v2(v2B, 3);

    double scp = xymDOT (v1, v2);
    cout << "scalar product: " << scp << endl;
    assert (fabs(scp-10)<1E-6);

    XymVector vX(3), vY;
    xymCOPY (v1, vY);
    xymSWAP (vX, vY);
    for (int i=0;i<8;i++) xymROT (vX, vY, cos(M_PI/4), sin(M_PI/4));
    cout << "vX:" << vX << " vY: " << vY << endl;
    assert(xymAMAX(vX)==2);
    assert (xymNRM2(vY)<1E-6);
    xymAXPY (-1, v1, vX);
    assert (xymASUM(vX)<1E-6);
    
    
}


void choleskyTest()
{
    int n=200;
    XymMatrixC A(n, n);
    A(0,0) = 100;
    for (int i=0; i<n-1; i++) {
        A(i,i) +=1;
        A(i+1,i+1) +=1;
        A(i,i+1) -= 1;
        A(i+1,i) -= 1;
    }
    XymMatrixC AI;
    choleskyInvert (AI, A);
    XymMatrixC res;
    multiply (res, A, AI);
    double sum=0;
    for (int i=0; i<n; i++) for (int j=0; j<n; j++) {
        if (i==j) sum += fabs(res(i,j)-1);
        else sum += fabs(res(i,j));
    }
    cout << "Error of " << sum << " in inverse of " << n << "*" << n << " matrix " << endl;
    assert (sum<1E-3);
}


int main (void)
{
    /*
    cout << "Polynomfit test" << endl;
    polynomfitTest();
    cout << " **************************** " << endl;

    cout << "Row / Column inverse update test" << endl;
    rowColumnTest();
    cout << " **************************** " << endl;
    */

    cout << "Cholesky inverse test" << endl;
    choleskyTest();
    cout << " **************************** " << endl;

/*
    cout << "BLAS test" << endl;
    blasTest();
    cout << " **************************** " << endl;
*/
}
