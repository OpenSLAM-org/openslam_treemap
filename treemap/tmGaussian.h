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
#ifndef TMGAUSSIAN_H
#define TMGAUSSIAN_H

/*!\file tmGaussian.h 
   \brief Class TmGaussian
   \author Udo Frese

   Contains the class \c TmGaussian representing a
  multidimensional Gaussian. 
*/


#include "xymatrix/xymMatrixC.h"
#include "tmTypes.h"
#include "tmExtendedFeatureId.h"

//! An n-dimensional Gaussian. All statistical information is stored in these objecs.
/*! We use homogenous coordinates adding a 1 as last entry to the
    vector. This way information matrix, information vector and
    overall likelihood, i.e. the 2nd, 1st and 0th coefficients of the
    quadratic log-likelihood functions are stored in a single matrix.
    To improve numerical stability and make the overall computation
    more elegant, we store the Cholesky factors \f$ RR^T \f$ of that
    matrix.

    \f{equation}
          \newcommand{\mv}[1]{{\begin{pmatrix} #1 \end{pmatrix}}}
          \exp\left(-\frac{1}{2}x^TAx + 2x^Tb + \gamma\right)
        = \exp\left(-\frac{1}{2}
                    \mv{x\\1}^T
                    \mv{A&b\\b^T&\gamma}
                    \mv{x\\1}^T
              \right)
        = \exp\left(-\frac{1}{2}
                    \mv{x\\1}^T
                    R^TR
                    \mv{x\\1}
              \right)
    \f}

    Assume a measurement / constraint \f$y = f(x) + err\f$ giving rise to
    the Gaussian

    \f{equation}
          \newcommand{\mv}[1]{{\begin{pmatrix} #1 \end{pmatrix}}}
         \exp\left(-\frac{1}{2} (f(x)-y)^TC^{-1}(f(x)-y)\right)
    \f}

    After linearization \f$ f(x) \approx f(\breve{x}) + J (x-\breve{x}) \f$
    with J being the Jacobian of f at the linearization point \f$\breve{x}\f$.

    \f{equation}
          \newcommand{\mv}[1]{{\begin{pmatrix} #1 \end{pmatrix}}}
         \exp\left(-\frac{1}{2} \left(f(\breve{x})-y + J (x-\breve{x})\right)^TC^{-1}
         \left(f(\breve{x})-y + J (x-\breve{x})\right) \right)
    \f}

    If C is decomposed as \f$C=LL^T\f$ the result can be written as

    \f{equation}
          \newcommand{\mv}[1]{{\begin{pmatrix} #1 \end{pmatrix}}}
         \exp\left(-\frac{1}{2} \mv{x\\1}^T\mv{J&f(\breve{x})-y-J\breve{x}}^TL^{-T}
           \underbrace{L^{-1}\mv{J&f(\breve{x})-y-J\breve{x}}}_R\mv{x\\1}\right)
    \f}.

    Several independent measurements / links can be combined by
    stacking R.  Any R can be transformed into an equivalent Gaussian
    with a triangular R using QR decomposition. For such a
    triangularized representation mean and covariance can be computed.
    R can be (conceptually) divided into a 2*2 bloack matrix. Then the
    right lower block is again triangular and represents the Gaussian
    where the features corresponding to the first block column are
    marginalized out (CIB). The first block row in turn represents the
    conditional distribution of the marginalized features given the
    remaining ones. Thus the ordering of features matters.

    The object contains a list which column corresponds to which
    feature. It is not necessarily sorted since the ordering of
    features depends on which feature shall be marginalized out next.

    The object is self contained and may be used inside or outside the
    treemap data structure.
*/
class TmGaussian 
{
 public:
  /*! Constructs an uninitialised Gaussian. \c isValid will return false. */
  TmGaussian();

  /*! Constructs an empty Gaussians representing the features in \c
      feature.  The matrix R is initialized with 0 rows but memory for
      \c rowsReserved is reserved, so the matrix may later grow when
      multiply is called. 
      
      \warning The counters in \c feature[i] represent the number of
      "information bits" integrated into this Gaussian. You can set it
      to >0 if you manually supply the Gaussians rows
      afterwards. However if you multiply other Gaussians into this
      Gaussian, you must initialize with 0 because multiply will add
      counters itself.
  */
  TmGaussian (const TmExtendedFeatureList& feature, int rowsReserved=0);  

  /*! Builds the Gaussian from a matrix \c R and features list
    \c feature. */
  TmGaussian (const XymMatrixVC& R, const TmExtendedFeatureList& feature, bool isTriangular);

  //! Constructs an uniform Gaussian with features specified in \c features
  /*! This is mainly for creating debug and test code. */
  TmGaussian (const char* features);  

  //! Same as constructor.
  void create (const TmExtendedFeatureList& feature, int rowsReserved=0);

  //! Same as constructor
  void create (const XymMatrixVC& R, const TmExtendedFeatureList& feature, bool isTriangular);

  //! Computes the marginal of \c this for columns \c fromFeature to \c cols()
  /*! The routine marginalizes out the features corresponding to columns
      \c 0..fromFeature-1. The gaussian must be in triangular form. Then this
      boils down to throwing away the first \c fromFeature rows and columns. 

      If the linearization point is marginalized out, it is set to -1. I.e.
      no relinearization is possible any more.
  */
  void computeMarginal (TmGaussian& result, int fromFeature);  
  


  //! Sets \c linearizationPointFeature and \c linearizationPoint
  /*! Most distributions involved in SLAM are rotation and translation
      invariant.  After linearization they are only invariant with
      respect to linearized rotation. Thus the treemap algorithm
      provides a mechanism for nonlinearily rotating distributions
      thereby implicitly changing the orientation used for
      linearization by a certain angle. When setting \c
      linearizationPointFeature, the algorithm uses the angle defined
      in the current estimate of this feature to nonlinearily rotate
      the distribution thereby reducing the linerarization error
      generated by error in the orientation. \c linearizationPoint is
      the angle used in computing the current linearization.

      The rotation itself is a member of treemap \c TmTreemap::rotateGaussian ().
   */
  void setLinearizationPoint (int linearizationPointFeature, double linearizationPoint);  

  
  /*! Makes R triangular by performing a QR decomposition.  This
      doesn't change the represented Gaussian but makes the
      representation more compact. It further allows to determine mean
      or covariance and perform marginalization.

      This operation is only possible if there are no duplicate
      features in \c feature.
  */
  void triangularize ();  

  //! Overloaded
  /*! If a workspace is provided, it is passed to \c xymGEQR2 and
      allocation and deallocation can be avoided.
  */
  void triangularize (XymVector& workspace);  

  /*! Multiplies \c this Gaussian by a marginalized distribution from
      another Gaussian overwriting \c this. \c gaussian must be in
      triangular form. Features \c gaussian.feature[0..fromFeature-1]
      are marginalized out. The rest is integrated into \c this. This
      Gaussian must represent all features integrated.  Multiplication
      means to add log-likelihood thus to stack a part of \c
      gaussian.R into \c this->R. The resulting Gaussian is not
      triangular.

      This routine even works, if \c gaussian contains duplicate features
      in which case the corresponding columns are added up. The counters
      in \c gaussian.feature are added to the counter in \c this->feature.
      This corresponds to the fact that the information from gaussian has
      been integrated into \c this.
  */
  void multiply (const TmGaussian& gaussian, int fromFeature=0);  

  
  /*! Computes the Gaussians mean and stores it into \c x. If \c
      \c upToFeature>=0 the mean is conditioned on \c feature[i] being
      \c x[i] for all \c i>=upToFeature. The Gaussian must be in
      triangular form.
   */
  void mean (XymVectorV& x, int upToFeature=-1);  

  //! Highly optimized version of \c mean for \c RCompressed
  /*! Performs (partial) backsubstitution on the matrix R stored in
      compressed form as \c RCompressed. This routine is highly
      optimized since it consumes the most part of computation time
      for very large maps. Thus it has a rather specialized interface:

      Assume that R has n=feature.size()+1 columns. Then the task is
      to compute \f$ x_i \f$ such that \f$ (Rx)_i=0 \f$ for all
      i=0..upToFeature-1. However for performance reasons both \c
      RCompressed and \c x are stored reverse, i.e. \c x[i]
      corresponds to \f$ x_{n-i-1} \f$.  x[0] must be set to 1, which is
      the implicity 1 of homogenous vectors. Then follow the features
      \c features[upToFeature..] upon which the mean shall be
      conditioned in reverse order. Then there follows empty space
      which is filled by \c meanCompressed for the features \c
      features[0..upToFeature-1] in reverse order. After that there
      must be 4 float entries available which are filled with 0s. This
      is necesary to allow the SSE optimized version to access four floats
      regardless of where in the vector.
   */
  void meanCompressed (float* x, int upToFeature);  


  /*! Removes all information. \c isValid will return false. */
  void clear();

  //! Computes \c RCompressed from \c R and clears \c R.
  /*! See \c RCompressed. */
  void compress ();  

  //! Asserts internal consistency
  void assertIt () const;  

  /*! Asserts that all elements are finite. */
  void assertFinite () const;  

  /*! Memory consumption in bytes */
  int memory () const;  

  /*! Returns \c true if a Gaussian is represented by this object or
      \c false if it is still uninitialised, e.g. after clear() or
      TmGaussian(). */
  bool isValid() {return R.isValid();}

  //! Whether R is a \c triangular matrix
  /*! If \c true, all R(i,j)==0 for i>j. Triangular matrices are so to
      say digested information, because an arbitrary number of
      measurements on n features can be transformed into a triangular
      n*n matrix. It is also digested in the sense, that one can
      readily compute means or conditioned means by backsubstitution.
  */
  bool isTriangular;  
  
  //! Number of rows in the Gaussian
  /*! This is not necessarily the number of features involved.
      Even not if the gaussian is in triangular form since it
      may rather be trapezoidal.
  */
  int rows() const {
    if (RCompressed.empty()) return R.rows();
    else return feature.size()+1; // RCompressed is a full triangle    
  }

  //! Number of cols in the Gaussian
  int cols() const {
    return feature.size()+1;
  }



  //! Copies everything from \c g2
  /*! Does the same as *this = g2 but copies the data by moving
      pointers only. As a consequence the data in \c g2 is set
      to NULL. This routine is useful to copy from local variables
      which are going to be destroyed anyway.
  */
  void transferFrom (TmGaussian& g2);  
  
 
  /*! Matrix defining the Gaussian as 

    \f{equation}
          \newcommand{\mv}[1]{{\begin{pmatrix} #1 \end{pmatrix}}}
       \exp\left(-\frac{1}{2}
                 \mv{x\\1}^T
                 R^TR
                 \mv{x\\1}
          \right)
    \f}

    R can be rectangular (more rows than column). It can be normalized
    to a triangle / trapezoidal matrix by QR decomposition. Such a matrix
    can be used to compute the mean and covariance.
  */
  XymMatrixC R;  

  //! Compressed version of \c R in a format that optimizes storage and computation time
  /*! \c R consumes the largest part of \c TmTreemap 's overall
      memory. So it can be compressed by calling \c ??? into a less
      consuming format. \c RCompressed stores the nonzero
      (upper-triangular) part of R rowwise and in single precision
      roughly reducing storage space by a factor of 4.  Actually it
      is stored backwards (decreasing i, j) and padding 0 to any missing triangle
      if R has more cols than rows.
      
      R(i,j) = RCompressed[n-j-1 + (n-i-1)*(n-i)/2]

      For technical reasons RCompressed is further padded with 3 0s.
      This allows the SSE optimized version of \c meanCompressed to always access 
      four floats regardless of where in the vector.
   */
  XycVector<float> RCompressed;  

  //! Accessing row \c i, column \c j in R stored in a compressed way
  const float& RCompressedAt (int i, int j) const
    {
      int n = feature.size()+1;
      return RCompressed [n-j-1 + (n-i-1)*(n-i)/2];
    }

  //! Accessing row \c i, column \c j in R stored in a compressed way
  float& RCompressedAt (int i, int j) 
    {
      int n = feature.size()+1;
      return RCompressed [n-j-1 + (n-i-1)*(n-i)/2];
    }

  //! Returns the \c RCompressed index of R(i,j)
  /*! Only valid for i<=j. */
  int RCompressedIdx (int i, int j) const
    {
      int n = feature.size()+1;
      return n-j-1 + (n-i-1)*(n-i)/2;
    }

  //! Returns the size of \c RCompressed if \c R.cols()==nCols
  static int rCompressedSize (int nCols) {return nCols*(nCols + 1)/2 + 3; }
  


  /*! column i of R corresponds to \c feature[i]. column feature.size() always
     corresponds to the fixed 1 entry. The counter reflects "how many bits"
     of information about this feature have been collected in the Gaussian.
     It is set by the aplication in the input Gaussian (usually 1) and added
     up when multiplying Gaussians.
  */
  TmExtendedFeatureList feature;  

  //! Feature (global index) which determines the linearization point
  /*! The algorithm provides a mechanism for reducing the
      linearization error caused by errors in orientation estimates.

      If \c linearizationPointFeature is -1 the Gaussian is not
      rotated at all.  Either because the algorithm is used for pure
      linear computation or because the "true" nonlinear probability
      distribution is not rotation/translation invariant. This always
      holds for distributions that define the first robot pose as 0
      and certainly for absolute position measurements such as GPS.

      The framework prevents a feature from being marginalized out as
      long as it is used as a linearizationPointFeature.
  */
  int linearizationPointFeature;  

  //! The angle used when linearizing \c .gaussian.
  /*! During propagation in the treemap distributions are rotated by
      the estimate of \c linearizationPointFeature and minus \c
      linearizationPoint, so they are rotated to fit most closely to
      the curent estimate. If the difference is larger than a
      threshold the node is invalidated due to excessive linearization
      error.
   */
  double linearizationPoint;  
};



#endif
