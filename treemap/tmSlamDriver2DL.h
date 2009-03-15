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
#ifndef TMSLAMDRIVER2DL_H
#define TMSLAMDRIVER2DL_H

/*!\file tmSlamDriver2DL.h 
   \brief Class \c TmSlamDriver2DL
   \author Udo Frese

  Contains the class \c TmSlamDriver2DL
  representing treemap driver for 2D landmark based SLAM.
*/

#include "tmTypes.h"
#include "tmTreemap.h"
#include "vectormath/vectormath.h"
#include <deque>


//! Treemap driver for 2D landmark based SLAM
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

    This class implements the most classical SLAM scenario. Robot
    poses are planar, 3-DOF (x, y, theta) and marginalized out except
    for the current pose. Landmarks are points, 2-DOF (x, y) and kept
    as the variables to be estimated. Odometry measurements are 3-DOF
    and relate a poses with respect to its predecessor. Landmark
    measurements are 2-DOF and in turn relate the landmark with
    respect to the robot pose where the observation was taken.

    Robot poses are marginalized out as appropriate and it is allowed
    to sparsify a pose out, if it shares at least two landmarks with
    its predecessor pose.

    For technical reasons the implementation requires an upper bound
    on the landmark id, i.e. on the number of landmarks. It reserves
    memory and assign feature indices to robot poses accordingly.
*/
class TmSlamDriver2DL : public TmTreemap
{ 
 public:
  TmSlamDriver2DL ();

  //! Initializes a new 2D landmark SLAM treemap
  /*! \c n is the number of landmarks and required to allocate memory
      and to assign numbers to robot pose features. \c reservePoses
      is used to make the feature array large enough to contain
      \c reservePoses poses in addition.
      
      The robot pose is set to \c initialPose with covariance \c
      initialPoseCov. See \c TmTreemap::sparsificationDistance. 
      
      See \c nrOfNonlinearLeaves and \c sparsificationDistance.
      See \c TmTreemap::Optimizer::maxNrOfUnsuccessfulMoves.
   */
  void create (int nrOfLandmarks, int reservePoses, const VmVector3& initialPose, const VmMatrix3x3& initialPoseCov, 
               int nrOfMovesPerStep=4, int maxNrOfUnsuccessfulMoves=3, int nrOfNonlinearLeaves=3, double sparsificationDistance=4);


  //! A single absolute pose measurement
  class AbsolutePose 
    {
    public:
      //! Feature index for the pose that is measurement (\c feature)
      int poseFeature;

      //! Pose's x, y, theta 
      VmVector3 pose;

      //! Covariance of \c pos
      VmMatrix3x3 poseCov;

      //! Returns, whether this measurement is empty
      bool isEmpty () const {return poseFeature<0;}
      

      //! Empty constructor
      AbsolutePose ()
        : poseFeature (-1)
        {
          vmZero (pose);
          vmZero (poseCov);
        }

      //! Std constructor
      AbsolutePose (int poseFeature, const VmVector3& pose, const VmMatrix3x3& poseCov)
        :poseFeature (poseFeature)
        {
          vmCopy (pose, this->pose);
          vmCopy (poseCov, this->poseCov);
        }                
    };  
  

  //! A single odometry measurement
  class Odometry 
    {
    public:
      //! Feature index for the old robot pose (\c feature)
      /*! see \c newPoseFeature */
      int oldPoseFeature;
      //! Feature index for the new robot pose (\c feature)
      /*! If newPoseFeature==oldPoseFeature==-1, the odometry measurement is empty. 
          If \c newPoseFeature>=0 && oldPoseFeature==-1, the measurement defines the
          absolute pose. */
      int newPoseFeature;

      //! New pose x,y,theta in the old pose's coordinate frame
      VmVector3 pose;

      //! Covariance of \c pos
      VmMatrix3x3 poseCov;

      //! whether the measurement is empty
      /*! This can be used as a flag to indicate a not existing odometry
	  measurement in a variable.
      */
      bool isEmpty () const {return newPoseFeature<0 && oldPoseFeature<0;}
          

      //! Empty constructor
      Odometry ()
        :oldPoseFeature(-1), newPoseFeature(-1)
        {
          vmZero (pose);
          vmZero (poseCov);
        }

      //! Std constructor
      Odometry (int oldPoseFeature, int newPoseFeature, const VmVector3& pose, const VmMatrix3x3& poseCov)
        :oldPoseFeature (oldPoseFeature), newPoseFeature (newPoseFeature)
        {
          vmCopy (pose, this->pose);
          vmCopy (poseCov, this->poseCov);
        }      
    };  
  

  //! The robot has moved by \c odometry with covariance \c odometryCov
  /*!
      \c odometry describes the new robot pose \f$ q \f$ in the frame of
      the old pose \f$ p \f$. According to the equation

      \f{equation}
                \newcommand{\mv}[1]{{\begin{pmatrix} #1 \end{pmatrix}}}
      \mv{q_x\\q_y\\q_\theta} = 
      \mv{p_x\\p_y\\p_\theta} +
      \mv{\cos p_\theta& -\sin p_\theta& 0\\\sin p_\theta& \cos p_\theta& 0\\0&0&1}
      \cdot
      \mv{o_x\\o_y\\o_\theta}

      \f}

      with \f$o\f$ being \c odometry.

      Since integrating an odometry measurement does not affect the
      landmark estimate but only the robot pose, odometry measurements
      are accumulated and only passed to the treemap algorithm if an
      observation is made.
  */
  void step (const VmVector3& odometry, const VmMatrix3x3& odometryCov);

  //! Tells the SLAM system that the robot is now on level \c level
  /*! Levels can be used to extend 2D mapping to layered 2.5D mapping for
      instance in multi-storey buildings. The level is just stored and
      not used in any further processing. However it can be used by the
      application for visualization, data association and to limit the
      update to one level.

      Also allocates an entry in \c level.
   */
  void setLevel (int level);  

  //! Adds absolute information about the robot pose
  /*! Can be GPS (x, y) or compass (theta). You can use infinity in
      the covariance diagonal to indicate that no information about a
      specific coordinate is given. 

      See Nebot et al. concerning how to handle the correlation in GPS
      error.
  */
  void absolutePose (const VmVector3& pose, const VmMatrix3x3& poseCov);  

  //! A single landmark observation
  class Observation 
    {
    public:
      //! Id of the landmark observed
      /*! If the id is \c <0 the observation is ignored. This feature can
          be used to handle spurious measurements. 

	  As can be seen here, the algorithm requires data-association.
      */
      int id;
      //! position relative to the robot
      VmVector2 pos;
      //! covariance of position measurement
      VmMatrix2x2 posCov;

      //! Zero / empty observation
      Observation ()
        :id(-1)
        {
          vmZero (pos);
          vmZero (posCov);
        }      

	//! Std constructor
      Observation (int id, const VmVector2& pos, const VmMatrix2x2& posCov)
        :id(id)
        {
          vmCopy (pos, this->pos);
          vmCopy (posCov, this->posCov);          
        }
    };

  //! List of observations from a single robot pose.
  typedef XycVector<Observation> ObservationList;  
  
  //! Integrate a list of observations from the current robot pose
  /*! The current implementation does not support to call \c observe
      several times without \c step.
  */
  void observe (const ObservationList& observation);

  //! Returns the \f$ \chi^2 \f$ error of \c observation in the current estimate
  double chi2 (const ObservationList& observation);  

  //! Returns the robot pose estimate
  void robotEstimate (double& x, double& y, double &theta);
  
  //! Returns the robot pose estimate
  void robotEstimate (VmVector3& pose);  

  //! Returns the estimate for landmark \c id
  /*! The return value specifies, whether there is a landmark with
    this id. */
  bool landmarkEstimate (int id, double& x, double& y)
  {
    int fId = 2*id+landmarkBaseFeature;    
    if (feature.idx(fId) && feature[fId].isDefined()) {
      x = feature[fId].est;
      y = feature[fId+1].est;
      return true;      
    }
    else return false;    
  }

  //!  Returns the estimate for landmark \c id
  /*! The return value specifies, whether there is a landmark with
    this id. */
  bool landmarkEstimate (int id, VmVector2& est)
  {
    int fId = 2*id+landmarkBaseFeature;    
    if (feature[fId].isDefined()) {
      est[0] = feature[fId].est;
      est[1] = feature[fId+1].est;
      return true;      
    }
    else return false;    
  }

  //! A robot pose as a random variable
  /*! This class contains the information about a robot pose
      that is not contained in the \c TmTreemap::feature vector.
      It is maintained in parallel to the feature vector.
  */
  class Pose 
    {
    public:
      //! Consecutive number in the odometry chain 
      /*! (-1) for unused. */
      int poseNr;

      //! Index in \c feature[]
      /*! x,y,theta are \c feature[featureId+0],feature[featureId+1],feature[featureId+2]
       */
      int featureId;

      //! Link to the next pose.
      /*! Index to \c pose. */
      int nextPoseId;

      //! Link to the previous pose.
      /*! Index to \c pose. */
      int prevPoseId;

      //! Link to the next pose in the same level
      /*! The poses of one level form a list with \c level[z].PoseId
          being the startpointer and \c nextPoseInSameLevel==-1 indicating
          the end of the list. The indices correspond to \c pose[..]
      */
      int nextPoseInSameLevel;      
      

      //! Level to which the pose belongs
      int level;      

      //! Accumulated distance travelled.
      /*! Used to implement the sparsification policy, that only
        one pose every \c sparsificationDistance may be sparsified out. 
      */
      double distance;      
    };

  //! A landmark as a random variable
  /*! This class contains the information about a landmark
      that is not contained in the \c TmTreemap::feature vector.
      It is maintained in parallel to the feature vector.
  */
  class Landmark 
    {
    public:
      //! Index in \c feature[]
      /*! x,y,theta are \c feature[featureId+0],feature[featureId+1],feature[featureId+2]
       */
      int featureId;

      //! Level to which the pose belongs
      int level;

      //! \c landmark[nextLandmarkInSameLevel] is the next landmark in the same level
      /*! The landmarks of one level form a list with \c level[z].landmarkId
          being the startpointer and \c nextLandmarkInSameLevel==-1 indicating
          the end of the list. The indices correspond to \c landmark[..]
      */
      int nextLandmarkInSameLevel;      

      Landmark ()
        :featureId (-1), level(-1)
        {}      
    };
  
  //! The poses currently represented in \c feature
  /*! \c pose[i] corresponds to \c feature[poseBaseFeature(i)] and following (x, y, theta).
   */
  XycVector<Pose> pose;  

  //! The landmarks currently represented in \c feature
  /*! \c landmarks[i] corresponds to \c feature[landmarkFeature(i)] and following (x, y)
   */
  XycVector<Landmark> landmark;

  //! A single level (for multi-level maps)
  /*! In general the estimation is 2D and levels just provide an index. The list of level
      is needed to traverse all poses/landmarks of a specific level.
  */
  class Level 
    {
    public:
      //! \c landmark[firstLandmarkInLevel] and \c ..[nextLandmarkInSameLevel] build the list of landmarks in a level
      int firstLandmarkInLevel;
      //! \c pose[firstPoseInLevel] and \c ...[nextPoseInSameLevel] build the list of poses in a level
      int firstPoseInLevel;      
      Level ()
        :firstLandmarkInLevel (-1), firstPoseInLevel(-1)
        {}      
    };
  
  //! List of levels
  /*! Basically only needed to traverse the list of all poses/landmarks in a single level. */
  XycVector<Level> level;  
  

  //! Additional feature flags
  enum FeatureFlag {LANDMARKX=0x1000000, LANDMARKY=0x2000000, POSEX=0x3000000, POSEY=0x4000000, POSETHETA=0x5000000};
  /*! \var LANDMARKX 
  
      x coordinate of a landmark (y coordinate is next)
  */
  /*! \var LANDMARKY 
 
      y coordinate of a landmark (x coordinate is prev)
  */
  /*! \var POSEX     

      x coordinate of a robot pose (y, theta are next)
  */
  /*! \var POSEY     

      y coordinate of a robot pose (theta is next, x prev)
  */
  /*! \var POSETHETA 

      theta orientation of a robot pose (X, y are prev)
  */
  
  //! Returns the maximal number of landmarks reserved
  int maxNrOfLandmarks () const  {return nrOfLandmarks;}
  
  //! Returns the number of landmarks that are actually estimated
  int nrOfLandmarksDefined () const;  

  //! Overloaded \c TmTreemap function
  virtual void nameOfFeature (char* txt, int featureId, int& n) const;  

  //! Overloaded
  virtual void clear();  

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
  virtual int memory () const;  


  //! Overloaded \c TmTreemap function
  /*! This routine implements the actual sparsification policy. It
      is as follows: A pose is sparsified out, if it is at least
      \c sparsificationDistance before the current robot pose and
      both leaves representing it share at least 2 landmarks. The
      first condition ensures, that the original leaves derived from 
      the measurements have already been integrated as far as the
      HTP algorithm is able to do before we sacrifice information
      by sparsifying. The second criterion prevents sparsifying out
      poses in a way that leads to disintegration of the map.
  */
  virtual bool canBeSparsifiedOut (TmFeatureId id);  

  //! Overloaded \c TmTreemap function
  /*! Checks all features marginalized at \c n whether they can be
      sparsified out.
  */
  virtual void checkForSparsification (TmNode* n);  

  //! Only compute the estimate for the landmarks \c [from..to-1]
  /*! See \c TmTreemap::onlyUpdatesEstimateFor
   */
  void onlyUpdateEstimatesForLandmarks (int from, int to, bool setDontUpdateFlag=true);  

  //! Only compute the estimate for the landmarks and poses of level \c level
  void onlyUpdateLevel (int level, bool setDontUpdateFlag=true);

  //! Overloaded \c TmTreemap function.
  virtual void assertIt () const;      


  //! Returns the feature number of landmark \c id
  int landmarkFeature (int id) const
    {return landmarkBaseFeature + 2*id;}

  //! Overloaded \c TmTreemap function
  virtual SlamStatistic slamStatistics () const;  

 protected:
  //! Statistic of the algorithm operation (nr of poses, etc.)
  SlamStatistic statistic;  
  
  //! Maximal landmark index
  int nrOfLandmarks;

  //! Landmark \c i has feature id \c \c landmarkBaseFeature+2*i
  int landmarkBaseFeature;
  
  //! Feature index at which the block of poses begins
  int poseBaseFeature;  

  //! Feature (with 2 successors) in the treemap corresponding to the last robot pose integrated
  /*! see \c relativePose */
  int poseFeature;

  //! Index in \c pose corresponding to \c poseFeature
  int poseIndex;  

  //! Ommit an odometry measurement at most every \c sparsificationDistance steps
  /*!
    If \c sparsificationDistance is set to infinity, no sparsification
    is performed and the computation is exact up to linearization
    (just as EKF). When treemap combines several poses in one leaf it
    can marginalize out internal poses but still has to represent one
    pose also used in other leaves every time the robot leaves or
    enters that leaf. With a finite \c sparsificationDistance the
    driver once in a while ommits an odometry measurement. Thereby the
    sequence of poses is cut into pieces. Whenever combines all
    measurements of such a piece in a single leaf, it can marginalize
    out all involved poses. So in the end treemap does not need to
    represent poses except some recent ones.

    The odometry measurement is ommitted whenever the current and
    previous observations share at least two landmarks and a distance
    of \c sparsificationDistance has been travelled since ommitting an
    odometry measurement the last time.

    The first criterion prevents that the map disintegrates into
    unrelated parts by cutting odometry. The second criterion ensures
    that overall not too much information is lost. As larger \c
    sparsificationDistance as less information is lost but as more
    computation is needed. So \c sparsificationDistance should be
    chosen large enough, that the odometry error is comparable with a
    typical relocalization error. It should not be larger than the
    distance where the total number of landmarks observed is
    comparable with the number of landmarks observed from a single
    pose.
  */
  double sparsificationDistance;  

  //! A leaf that represents odometry and observations at a certain pose
  /*! This class is derived from \c TmNode and contains the original measurements
      that led to this node in the treemap. Thus it allows to update the linearization
      according to the current estimate.
   */
  class NonlinearLeaf:public TmNode 
    {
    public:
      //! Creates a default leave with \c tree set.
      NonlinearLeaf (TmSlamDriver2DL* tree);

      //! overloaded function virtually calling the copy constructor
      virtual TmNode* duplicate() const;      

      //! Returns the storage space (Bytes) of this node without children
      virtual int memory() const;  

      //! All measurements have been made at the pose corresponding to this feature.
      int poseFeature;      

      //! Original absolute measurement internal into this leaf
      /*! Maye empty. */
      AbsolutePose absolutePose;      

      //! Original odometry measurement integrated into this leaf
      /*! Maybe empty for leaves where the odometry sequence has been cut. */
      Odometry odometry;
      
      //! Original landmark measurements integrated into this leaf
      ObservationList observation;
      
      //! Compute \c .gaussian by linearizing the \c odometry and \c observation
      /*! If \c useMeasurement is false, the current estimate is used
          as linearization point. This is done when relinearizing a
          node. Otherwise the measurement itself provides part of the
          linearization point. This is used when the leaf is created
          in the beginning.

          If \c tree!=NULL, i.e. the node is already part of the treemap, the
          routine calls \c beforeChange and \c afterChange appropriately.
       */
      void linearize (bool useMeasurement);      


    protected:
      //! Adds a 2-DOF landmark measurement to \c gaussian by appending 2 rows
      /*! The routine assumes that the pose corresponds to rows cP..cP+2 and the landmark
        to columns cL..cL+1. Uses the current \c feature estimate for linearization.
        
        If \c useMeasurement is true, the routine uses the measurement
        \c pos itself for part of the linearization. This is better when
        linearizing a measurement the first time because the prior
        estimate may be subject to error accumulation. When
        relinearizing the estimate already incorporates the measurement and using the estimate
        is better because it results in consistent linearization points.
      */
      void addLandmarkMeasurement (int cP, int cL, const VmVector2& pos, const VmMatrix2x2& posCov, bool useMeasurement);  
      
      //! Adds a 3-DOF odometry measurement to \c gaussian by appending 3 rows
      /*! The routine assumes, that columns \c c0..c0+2 correspond to the old pose (x,y,theta)
        and columns \c c1..c1+2  correspond to the new pose (x, y, theta). Uses the currrent
        \c feature estimate for linearization. 
        
        For first linearization of an odometry measurement the estimate for the new pose is
        initialized by the odometry measurement itself. So we always use the measurement and
        do not need a \c useMeasurement flag as with \c addLandmarkMeasurement.        
      */
      void addOdometry (int cOld, int cNew, const VmVector3& odo, const VmMatrix3x3& odoCov);

      //! Adds absolute information on the robot pose
      void addAbsolute (int cP, const VmVector3& pos, const VmMatrix3x3& posCov);      
    };  
  

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

  //! The current estimate includes \c nonlinearLeaf[0..estimateIncludes]
  /*! Before computing an estimate the nonlinear leaves are updated using
      the current estimate as linearization point. This is only good for
      those leaves that are already incorporated into the estimate. New leaves
      use the measurements itself as linearization point, which can be much
      better. Thus they should not be relinearized.
  */
  int estimateIncludes;    

  //! Current robot pose relative to \c poseFeature.
  /*! We pass both odometry and landmarks in a single Gaussian to the treemap.
      Further we accumulate steps without observation.
  */
  VmVector3 relativePose;

  //! Covariance of \c relativePose
  VmMatrix3x3 relativePoseCov;  
  
  //! The distance travelled since poseFeature
  double relativeDistance;  

  //! reset \c relativePose and \c relativePoseCov
  void clearRelativePose ();  

  //! Level on which the robot currently resides
  /*! The level is stored in all robot poses (\c Pose) and 
      landmarks (\c Landmark). If a landmark is observed from several
      levels the first one decides.
  */
  int currentLevel;  
  

  //! Counts the number of landmarks both in \c a and \c b
  static int sharedLandmarks (const ObservationList& a, const ObservationList& b);  

  //! Counts the number of landmarks both in \c a and \c b
  int sharedLandmarks (const TmExtendedFeatureList& a, const TmExtendedFeatureList& b);  

  //! Adds a new pose both to the odometry chain
  /*! Returns the feature id of that pose. Initializes the flags. Sets
      pose number and accumulated distance, where \c distance is the
      distance travelled since adding the last pose. If there is no
      pose, the distance is 0 and the \c Pose::poseNr=-1. 
      Set the \c CAN_BE_SPARSIFIED flag, if this is okay according to
      \c latestSparsificationDistance.
  */
  int newPose (double distance);

  //! \c TmTreemap function overloaded to maintain \c pose.
  virtual void deleteFeature (TmFeatureId id);  

  //! Overloaded \c TmTreemap function
  /*! Uses \c pose to reset the \c CAN_BE_SPARSIFIED flag in all
      poses within \c sparsificationDistance range.
  */
  virtual void hasBeenSparsifiedOut (TmFeatureId id);  


  //! Sets the flags for the 3 features of a pose.
  /*! fId..fId+2 are all set to \c flags and \c POSEX, \c POSEY, or \c POSETHETA respectively
   */
  void setFlagsForPose (int fId, int flags);
  
  //! Sets the initial estimate of \c fId..fId+2 to \c pose
  void setInitialEstimateForPose (int fId, const VmVector3& pose);  

  //! Initializes landmark \c id
  /*! Sets \c landmark[id] and the landmark's feature flags
   */
  void initLandmark (int id);  

  //! Sets the flags for the 2 features of a landmark
  /*! fId..fId+1 are all set to \c flags and \c LANDMARKX or \c LANDMARKY respectively
   */
  void setFlagsForLandmark (int fID, int flags);  
};


#endif
