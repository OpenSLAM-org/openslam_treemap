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
//!

/*!\file tmSlamDriver2DP.h 
   \brief Class TmSlamDriver2DP
   \author Udo Frese 

  Contains the class \c TmSlamDriver2DP
  representing treemap driver for 2D pose based SLAM (Lu & Milion style).
*/
#ifndef TMSLAMDRIVER2DP
#define TMSLAMDRIVER2DP

#include "tmTypes.h"
#include "tmTreemap.h"
#include "vectormath/vectormath.h"


//! Treemap driver for 2D pose relation based SLAM
/*! This class implements the most simple SLAM scenario. Random variables
    to be estimated are a set of 2D/3-DOF poses (x, y, theta). These poses
    are links by 3-DOF relative links which are the observations put into
    the system. No random variables are marginalized out.
*/
class TmSlamDriver2DP : public TmTreemap
{
 public:
  //! empty constructor
  TmSlamDriver2DP ();
  
  
  //! Initializes a new 2D pose relation based SLAM treemap
  /*! The system reserves memory of \c nrOfPosesReserved. It is okay to have more poses
      but this will trigger copying of some vectors involved.

      See \c TmTreemap for the parameters.
  */
  void create (int nrOfPosesReserved=0
, int nrOfMovesPerStep=4, int maxNrOfUnsuccessfulMoves=3);

  //! A single link, i.e. uncertain information on the relative pose of two poses
  /*!
      Such a link corresponds to the measurement, that the relation pose of \c poseA
      relative to \c poseB is \c d with uncertainty covariance \c dCov. The correspoding
      measurement equation is:

      \f{equation}
                \newcommand{\mv}[1]{{\begin{pmatrix} #1 \end{pmatrix}}}

      \mv{d_x\\d_y\\d_\theta} = 
      \mv{\cos p_\theta& \sin p_\theta& 0\\-\sin p_\theta& \cos p_\theta& 0\\0&0&1}
      \cdot
      \left(
      \mv{a_x\\a_y\\a_\theta} - \mv{b_x\\b_y\\b_\theta}
      \right)

      \f}
  */
  class Link 
  {
  public:
    //! The first pose involved
    /*! see \c poseB */
    int poseA;

    //! The second pose involved
    /*! The link gives information on pose nr. \c poseA relative to
        pose nr \c poseB.  \c poseB can be -1 refering to the world
        coordinate system. 
    */
    int poseB;    

    //! Position and Orientation of \c poseA in \c poseB coordinates
    VmVector3 d;

    //! The uncertainty (i.e. covariance) of \c d
    VmMatrix3x3 dCov;

    //! Empty link
    Link () :poseA(-1), poseB(-1) {vmNaN(d); vmNaN(dCov);} 
      
    //! A link between \c poseA and \c poseB with measurement \c d, \c dCov
    /*! A link represents the probabilistic information on the
      pose nr. \c poseB relative to pose nr. \c poseA. The information is,
      that it is \c d with covariance \c dCov.
      
      The exact meaning of \c d is documented in \c Link.
    */
      Link (int poseA, int poseB, const VmVector3& d, const VmMatrix3x3& dCov)
	:poseA(poseA), poseB(poseB)
	{
	  vmCopy (d, this->d);
	  vmCopy (dCov, this->dCov);	  
	}	

	//! Returns the pose with the larger index
	int largerPose () const
	{
	  if (poseA>poseB) return poseA;
	  else return poseB;
	}	
  };
  

  //! Adds uncertain information about the relative location of two poses to the system
  /*!
      If the link involves poses which are not there yet, they are added to the SLAM map.
      At least one of the poses must already be there in order to compute an initial
      estimate. The estimate is linearized using the measurement d and the orientation
      of the initial estimate as linearization.

      \c poseB can be -1, in which case the link defines the absolute pose of \c poseA
      with some uncertainty. This is necessary to define the coordinate origin (otherwise
      the system is singular) and it can also used to incorporate GPS/compass data. 

      See Nebot et al. concerning how to handle the correlation in GPS
      error.
          
   */
  void addLink (const Link& link);

  //! Returns the estimate of pose \c idx
  void poseEstimate (int idx, double& x, double& y, double &theta) const;

  //! Returns the estimate of pose \c idx
  void poseEstimate (int idx, VmVector3& pose) const
  {poseEstimate (idx, pose[0], pose[1], pose[2]);}  

  //! Returns the robot pose estimate
  /*! The pose with the largest index is interpreted as the current robot pose. */
  void robotEstimate (double& x, double& y, double &theta) const
  {poseEstimate (pose2Feature.size()-1, x, y, theta);}  
    
  
  //! Returns the robot pose estimate
  /*! The pose with the largest index is interpreted as the current robot pose. */
  void robotEstimate (VmVector3& pose) const  
  {poseEstimate (pose2Feature.size()-1, pose);}  


  //! Flags used to tag a 1D feature as x, y or theta
  /*! These flags are used together with the \c TmFeature::FeatureFlags flags. 
      All user flags must be part of \c TmFeature:USER_FLAGS.
   */
  enum FeatureFlag {POSEX=0x1000000, POSEY=0x2000000, POSETHETA=0x3000000};
  /*! \var POSEX     
 
      x coordinate of a robot pose (y, theta are next)
  */
  /*! \var POSEY     

      y coordinate of a robot pose (theta is next, x prev)
  */
  /*! \var POSETHETA 

      theta orientation of a robot pose (X, y are prev)
  */

  //! Overloaded \c TmTreemap function
  virtual void nameOfFeature (char* txt, int featureId, int& n) const;  

  //! Overloaded
  virtual void clear();  

  //! Overloaded \c TmTreemap function
  virtual int memory () const;  
  

// protected:
  //! \c pose2Feature[i] is the treemap-feature id of pose number \c i
  /*! If these pose has not been used yet it is -1 (or beyond the \c 
      i>=poseFeatureId.size())
  */
  XycVector<int> pose2Feature;

  //! Number of poses
  int p;  

  

 //! Computes a linearization of the measurement residual for a link
 /*! \c poseALinPoint and \c poseBLinPoint are the poses for a and b used as linearization
     point. \c d is the actual measurement of the relative pose. So in linear approximation
     at \c a=poseALinPoint and \c b=poseBLinPoint the \f$ chi^2 \f$ error is
 
     \f[
          |A \cdot (a,b,1)|_2^2
     \f]
 
     which is the representation used by treemap.
 
     This is equation (6) of Frese, Larsson, Duckett, 2004.

     If \c includeJacobianForB==true, the 3 columns corresponding to b are ommited. This is
     needed for links to the fixed ground frame -1.
  */
  static void linearizeLink (XymMatrixC& A, const VmVector3& poseALinPoint, 
			     const VmVector3& poseBLinPoint, const VmVector3& d,
			     bool includeJacobianForB=true); 

  //! reserves and initializes the pose random variables of a link
  /*! The routine assures, that for both \c poseA and \c poseB features are existing.
      At least one of them must have an initialized estimate and the other one is
      computed from that and from the link if necessary.
  */
  void setInitialEstimate (const Link& link);  

  //! Allocates feature entries for a new pose
  /*! A block of 3 features (x, y, theta) is reserved at \c TmTreemap and the first feature
      index is stored in \c pose2Feature So to say, this function registers the pose \c id at the 
      treemap kernel. */
  void allocatePose (int id);  
};

#endif
