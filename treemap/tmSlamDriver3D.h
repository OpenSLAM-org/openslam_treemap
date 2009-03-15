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
#ifndef TMSLAMDRIVER3D_H
#define TMSLAMDRIVER3D_H

/*!\file tmSlamDriver3D.h 
   \brief Class \c TmSlamDriver3D
   \author Udo Frese
   
    Contains the class \c TmSlamDriver3D, a driver for 3D
    landmark based SLAM with 3D landmarks and 3D/6DOF poses, but without
    any connection between successive poses.
*/

#include "tmTypes.h"
#include "tmTreemap.h"
#include "vectormath/vectormath.h"
#include <map>
#include <deque>

//! A pose in 3D space
extern VmMatrix4x4 globalPose;

//! Treemap driver for 3D landmark based SLAM without odometry
/*! The treemap algorithm itself is very generic mostly viewing the
    problem as least square estimation of some random
    variables. Especially it does not use the notion of landmarks and
    robot poses and makes no assumption on their
    parametrization. Neither does it assume a specific motion or
    sensor model. It also does not assume a specific sparsification
    policy.

    Thus for any specific scenario a driver class is needed that
    translates the measurements into appropriate Gaussians and manages
    the features corresponding to landmarks and / or robot poses.

    This class implements 3D feature based SLAM without odometry. This
    means, a measurement is a 3D vector giving the position of an
    identified feature relative to the freely floating robot. There is
    no information connecting successive poses, except from the
    features.

    Robot poses are marginalized out. No further approximation is
    necessary, because, as there is no connection between poses,
    they can be marginalized out without affecting sparsity.
*/
class TmSlamDriver3D : public TmTreemap
{
 public:
  //! Empty constructor
  TmSlamDriver3D ();
  
  //! Destructor
  ~TmSlamDriver3D ();

  //! Initializes a new 3D SLAM treemap
  /*! The robot pose is set to \c initialPose 
      
      See \c nrOfNonlinearLeaves.
      See \c TmTreemap::create for \c nrOfMovesPerStep and \c maxNrOfUnsuccessfulMoves.
   */
  void create (const VmMatrix4x4& initialPose, int nrOfMovesPerStep, int maxNrOfUnsuccessfulMoves, int nrOfNonlinearLeaves);
  

  //! Observation of a landmark from a robot pose
  /*! A \c LandmarkObservation gives the information, that the position
      of \c landmarkId relative to robot \c poseId is \c pose with
      an uncertainty of \c posCov.
  */
  class LandmarkObservation 
 {
  public:
    //! Which landmark has been observed (userId)
    int landmarkId;
    //! From which robot pose (userId)
    int poseId;    
    //! Observed position in robot coordinates
    VmVector3 pos;
    //! Observation covariance
    VmMatrix3x3 posCov;

    //! An empty landmark observation
    LandmarkObservation ()
      :landmarkId(-1), poseId(-1)
      {
        vmZero (pos);
        vmZero (posCov);
      }

      //! A landmark observation with all parameters given
      /*! Except \c poseId which is -1
       */
      LandmarkObservation (int id, const VmVector3& pos, const VmMatrix3x3& posCov)
        :landmarkId(id), poseId(-1)
        {
          vmCopy (pos, this->pos);
          vmCopy (posCov, this->posCov);          
        }
        

      //! Computes the jacobian of the measurement function f
      /*! The Jacobian with respect to the landmark position
          is returned in \c JPL. The jacobian with respect
          to a 6-DOF parametrization of the pose is returned
          in \c JPP. The Jacobian is evaluated at a position
          for the landmark \c pL0 and a robot pose \c T0.
      */
      static void jacobian (double JPL[3][3], double JPP[3][6],
                            const VmVector3& pL0, VmMatrix4x4& T0);

      //! Linearizes the observation and adds the resulting 3-DOF Gaussians to \c G
      /*! The observation equation is linearized at \c pL0 concerning the landmark position
          and \c T0 concerning the robot pose. Then a Gaussian is made from the Jacobian,
          the measurement \c pos and its covariance \c posCov and the result is added as
          3 rows to \c G. Column \c jPose..jPose+5 correspond to the pose and column \c jLm..jLm+2
          to the landmark.

          If either \c jPose or \c jLm is -1 the corresponding entries are not stored in G which
          is equivalent to condition on the landmark being equal to \c pL0 or the pose being equal
          to \c T0 respectively.
      */
      void addToGaussian (TmGaussian& G, int jPose, int jLm, const VmVector3& pL0, VmMatrix4x4& T0) const;      
  };

  //! Vector of landmark observations
  typedef XycVector<LandmarkObservation> LandmarkObservationList;  

  //! multi DOF random variable being a meaningful entity
  /*! For instance a landmark, a robot pose, etc. */
  class RandomVariable 
  {
  public:
    //! Id in \c TmTreemap 's feature array
    int featureId;   

    //! Id passed by the calling application
    int userId;

    //! Pointer to the overall tree class
    TmSlamDriver3D* tree;    

    //! Dimension
    /*! \c feature[featureId..featureId+dimension-1] are the 1D variables
        that correspond to this multi-dimension variable. */
    int dimension;

    //! Empty random variable
    RandomVariable () :featureId(-1), userId(-1), tree(NULL), dimension(0) {}

      //! Random variable with a given user id and dimension
      RandomVariable (int userId, int dimension) 
        :featureId(-1), userId(userId), tree(NULL), dimension(dimension) {}
      
	//! destructor
    virtual ~RandomVariable();

    //! Returns the memory consumption of \c this
    virtual int memory () const;    
  };
  
  

  //! A 3D point landmark as a random variable
  /*! A landmark is a random variable, i.e. a quantity estimated
      by treemap. Hence it is derived from \c RandomVariable,
      so both landmarks and poses can be kept in a single
      random variable container \c variablesByUserId.
  */
  class Landmark : public RandomVariable
  {
  public:
    //! An empty landmark with no id yet
    Landmark (): RandomVariable() {dimension = 3;}

      //! A Landmark with a given \c userId
    Landmark (int userId) : RandomVariable (userId, 3) {}

      //! destructor
    virtual ~Landmark ();      
        
      
    //! Read the estimate from the corresponding feature entries
    void getEstimate (VmVector3& p) const;
      
    //! Sets an initial estimate to \c p
    void setInitialEstimate (const VmVector3& p);    

    //! Overloaded \c RandomVariable function
    virtual int memory () const;    
  };  

  //! A map of \c RandomVariable index by their userId
  typedef std::map<int, RandomVariable*> IntRVMap;  

  //! Resets to no feature, no information
  virtual void clear();  

  //! Proceed to the next robot pose
  void step ();
  
  //! Incorporate a set of observations from a single pose
  /*! All */
  void observe (const LandmarkObservationList& obs); 

  //! Overloaded \c TmTreemap function
  /*! Sets \c estimateIncludes to reflect that all leaves are now
      incorporated into the estimate.
  */
  virtual void computeLinearEstimate ();
  
  //! Overloaded \c TmTreemap function
  /*! Relinearizes all leaves in \c nonlinearLeaf and calls
      \c computeLinearEstimate afterwards.
   */
  virtual void computeNonlinearEstimate ();

  //! Overloaded \c TmTreemap function
  virtual SlamStatistic slamStatistics () const;  

  //! Overloaded \c TmTreemap function
  virtual int memory () const;  

  //! Return the random variable corresponding to \c userId
  /*! If there is none with that id, \c NULL is returned
   */
  RandomVariable* variableFromUserId (int userId) const {
    IntRVMap::const_iterator it = variablesByUserId.find (userId);
    if (it==variablesByUserId.end()) return NULL;
    else return (*it).second;
  }  

  //! Return the landmark corresponding to \c userId
  /*! If there is none with that id, \c NULL is returned
      This is just the same as \c variablesByUserId followed
      by a dynamic cast to \c Landmark*
   */
  Landmark* landmarkFromUserId (int userId) const {
    return dynamic_cast<Landmark*> (variableFromUserId (userId));    
  }

  //! Different multi-DOF variables to be estimated indexed by \c .userId
  IntRVMap variablesByUserId;

  //! Whether we are now the pose defined by \c initialPose
  bool isFirstPose () const {return isFirstPoseX;}
  

protected:
  //! Statistic of the algorithm operation (nr of poses, etc.)
  SlamStatistic statistic;  

  //! Different multi-DOF variables to be estimated indexed by \c .featureId
  IntRVMap variablesByFeatureId;  


  //! A leaf containing a set of measurements from a single robot pose
  /*! This class is derived from \c TmNode and contains the original measurements
      that led to this node in the treemap. Thus it allows to update the linearization
      according to the current estimate.
   */
  class NonlinearLeaf:public TmNode 
    {
    public:
      //! Creates a default leave with \c tree set.
      NonlinearLeaf (TmSlamDriver3D* tree);

      //! Creates a leave with \c tree set and the passed observations.
      NonlinearLeaf (TmSlamDriver3D* tree, const LandmarkObservationList& obs);

      //! Creates a leave with \c tree set and the passed observations having a known \c pose
      NonlinearLeaf (TmSlamDriver3D* tree, const LandmarkObservationList& obs, const VmMatrix4x4& knownPose);

      //! overloaded function virtually calling the copy constructor
      virtual TmNode* duplicate() const;

      //! overloaded \c TmNode function
      virtual int memory() const;  

      //! Observations incorporated in this leaf
      LandmarkObservationList obs;

      //! Whether this leaf has already been incorporated into the estimate
      /*! The flag is set by \c computeLinearEstimate. It is used to
          decide whether it is better to linearize based on the
          estimate (\c hasBeenIncorporated=true) or based on the observation
          (\c hasBeenIncorporated=false). 
      */
      bool hasBeenIncorporated;

      //! Whether this observation set was taken from \c knownPose
      /*! Usually only the first set of observations corresponding to the
        initial pose has a \c knownPose. */
      bool hasKnownPose;

      //! If \c hasKnownPose the pose from which this observation set was taken
      VmMatrix4x4 knownPose;      
      

      //! Computes an initial estimate for the robot \c pose
      /*! Matches all observations  to the current
          estimate (if there is one) and computes the least square \c pose estimate
          therefrom (using \c vmLeastSquareTransform). If \c allLandmarks is \c false
          only landmarks in \c landmarks are considered for matching.

          There must be at least 3 overlapping landmarks in the observation
          and in \c landmarks.
      */
      void estimatePose (VmMatrix4x4& pose, bool allLandmarks=true, const XycVector<int>& landmarks=XycVector<int>());

      //! Sets an initial estimates for all involved features that have no estimate yet
      /*! by assuming that \c pose is the current observer pose. \c pose is used regardless
          if \c hasKnownPose.
      */
      void setInitialEstimate (const VmMatrix4x4& pose);      

      //! Compute \c .gaussian by linearizing \c obs
      /*! \c pose is used as a linearization point concerning the
          robot pose.  If \c useMeasurement is false, the current
          estimate is used for defining the linearization point
          regarding the landmarks positions. This is done when
          relinearizing a node. Otherwise (\c useMeasurement=false)
          the measurement itself provides the linearization point for
          the landmark positions. This is used when the leaf is
          created in the beginning because it provides a dramatically
          better linearization point when closing a loop.

          If \c hasKnownPose, the information is added that the current
          pose is exactly \c pose. This is used for the first observations
          to define the coordinate system. Otherwise the \c pose defines no
          statistical information but just the linearization point. Note, that
          if \c hasKnownPose it only makes sense to pass \c knownPose as \c pose.

          If \c index>=0, i.e. the node is already part of the treemap, the
          routine calls \c beforeChange and \c afterChange appropriately.
       */
      void linearize (VmMatrix4x4& pose, bool useMeasurement);

      //! Linearize base on the current landmark estimates
      /*! This routine should be called when relinearizing, i.e. when
          the leaf is already incorporated into the estimate. If not
          \c hasKnownPose the routine determines a pose estimate based
          on all observed landmarks otherwise it takes \c knownPose as
          linearization point. It uses that to linearize all measurements.

          If \c index>=0, i.e. the node is already part of the treemap, the
          routine calls \c beforeChange and \c afterChange appropriately.
      */        
      void linearize ();      
    };

  //! First pose defining the coordinate system
  /*! Only stored as a buffer until it is passed into the first leaf. */
  VmMatrix4x4 initialPose;  

  //! True if the current pose is the first one
  bool isFirstPoseX;

  //! We have 6 reserved feature indices for a dummy pose
  /*! The pose is only used intermediate in forming a Gaussian and marginalized
      out before passing the Gaussian to \c TmTreemap. But still we need indices
      that are not used by any other feature, so we reserve 6 DOF in the beginning.
  */
  int dummyPoseId;  
  
  //! Features that have been observed in the last pose
  /*! This list is used to compute an initial estimate for the robot
      pose by matching features that have been observed recently. 
      The stored numbers are user ids.
   */
  XycVector<int> recentlyObservedLandmarks;  

  //! Keep the observations at the last \c nrOfNonlinearLeaves poses nonlinearily.
  /*! Linearization errors are caused by different measurement being
      linearized at different linearization points. To reduce this problem
      a delayed state approach is persued. The last \c nrOfNonlinearLeaves are
      kept as \c TmSlamDriver2DL::NonlinearLeaf with their original measurements
      and are relinearized each round.
   */
  int nrOfNonlinearLeaves;

  //! List of nonlinear leaves in the tree.
  deque<int> nonlinearLeaf;



  

  //! Registers \c rv as a new random variable to be estimated
  /*! Allocates feature entries in \c feature, sets \c rv->featureId
      and adds \c rv to \c variablesByUserId and \c
      variablesByFeatureId. After calling \c register, \c this owns
      \c rv and will take care about freeing.
   */
  void registerRandomVariable (RandomVariable* rv);
  
  //! \c TmTreemap function overloaded
  /*! Removes \c id from all \c TmSlamDriver3D lists.
   */
  virtual void deleteFeature (TmFeatureId id);  

  //! Sets \c recentlyObservedLandmarks as the landmark ids in \c obs
  void setRecentlyObservedLandmarks (const LandmarkObservationList& obs);  
};


#endif
