//!\author Udo Frese

/*!\file xymVectorV.cc Contains the implementation for class \c
  XymVectorV representing a vector view.
*/

#include "xymVectorV.h"
#include "xymVector.h"
#include "xymOperations.h"
#include <math.h> 

XymVectorV::XymVectorV ()
        :_d(NULL), _n(0), _stride(0)
{}

XymVectorV::XymVectorV (double* d, int n, int stride)
        :_d(d), _n(n), _stride(stride)
{}


bool XymVectorV::isFinite() const
{
    double *p, *pE;
    int pInc;
    loop (p, pE, pInc);
    while (p!=pE) {
      // if isfinite is not found, maybe, you have to switch on C99
        if (!isfinite(*p)) return false;
        p+=pInc;
    }
    return true;
}


void XymVectorV::extractFrom (const XymVectorV& v, const vector<int>& idx)
{
    XYMFORMAT (size()==(int) idx.size());
    const int* idxP = &idx[0];
    double *p, *pE;
    int pInc;
    loop (p, pE, pInc);
    while (p!=pE) {
        *p = v[*idxP];
        p+=pInc;
        idxP++;
    }
}


void XymVectorV::copyFrom (const XymVectorV& v)
{
    XYMFORMAT (size()==v.size());
    double *src, *dst, *srcEnd;
    int srcOfs, dstOfs;
    loop (dst, dstOfs);
    v.loop (src, srcEnd, srcOfs);
    while (src!=srcEnd) {
        *dst = *src;
        src += srcOfs;
        dst += dstOfs;
    }
}


void XymVectorV::fill (double alpha)
{
    double *src, *srcEnd;
    int srcOfs;
    loop (src, srcEnd, srcOfs);
    while (src!=srcEnd) {
        *src = alpha;
        src += srcOfs;
    }    
}


void XymVectorV::print () const 
{
    double max=1E-32;
    for (int i=0; i<size(); i++) if (max<fabs((*this)(i))) max=fabs((*this)(i));
    max*=10;
    double scale = exp(log(1000.0)*floor(log(max)/log(1000.0)));
    
    printf("%3.1e*{ ", scale);
    for (int i=0; i<size(); i++) {
        printf("%+7.3f", (*this)(i)/scale);
        if (i<size()-1) printf(", ");
    }
    printf(" }\n");
}


ostream& operator<<(ostream& os, const XymVectorV& v)
{
    os << " { ";
    if (v._d!=NULL) for (int i=0; i<v._n; i++) {
        os << v._d[i];
        if (i<v._n-1) os << " ";
    }
    os << " } ";
    return os;
}


istream& operator>>(ifstream& ifs, XymVectorV& v) throw (runtime_error)
{
    xymExpect (ifs, "{");
    if (v._d==NULL) throw runtime_error ("XymMatrixV::operator>> v is empty");
    int ctr=0;
    string buffer;
    ifs >> buffer;
    if (buffer!="}") {
        do {
            ctr++;
            if (ctr>v._n) throw runtime_error ("XymMatrixV::operator>> too many entries in vector");
            v[ctr] = atof(buffer.c_str());
            ifs >> buffer;
            if (buffer==",") ifs >> buffer;
            else break;
        } while (true); // break
    }
    if (ctr<v._n) throw runtime_error ("XymMatrixV::operator>> too little entries in vector");
    return ifs;
}


void print (const XymVectorV& v)
{
    v.print();
}

XymVector const operator+ (const XymVectorV& a, const XymVectorV& b)
{
  XymVector result;
  add (result, a, b);
  return result;
}


void operator+= (XymVectorV& a, const XymVectorV& b)
{
  add (a, a, b);
}


XymVector const operator- (const XymVectorV& a, const XymVectorV& b)
{
  XymVector result;
  sub (result, a, b);
  return result;
}
  

void operator-= (XymVectorV& a, const XymVectorV& b)
{
  sub (a, a, b);
}


XymVector const operator* (const XymVectorV& a, double lambda)
{
  XymVector result;  
  scale (result, a, lambda);
  return result;
}


XymVector const operator* (double lambda, const XymVectorV& a)
{
  XymVector result;  
  scale (result, a, lambda);
  return result;
}


void operator*= (XymVectorV& a, double lambda)
{
  scale (a, a, lambda);
}


double operator* (const XymVectorV& a, const XymVectorV& b)
{
  return dot (a, b);
}

