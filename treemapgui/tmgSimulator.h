#ifndef TMG_SIMULATOR_H
#define TMG_SIMULATOR_H

//!\author Udo Frese

/*!\file tmgSimulator.h Contains the class \c TmgSimulator that can simulate robot
   motion and landmark perception in a multi-story building specified by a bitmap.
*/

#include <vector>
#include "vectormath/vectormath.h"
#include "vectormath/vmRandom.h"
#include <string>
using namespace std;

//! Simulates robot motion and landmark perception in a multi-story building specified by a bitmap.
/*
   In the bitmap the robot moves along a red path. Walls are black and
   prevent the robot from observing landmarks. Landmarks are
   blue. Elevators green. The simulator can use the same bitmap for
   many stories so it is able to efficiently simulate a million
   landmarks map.

   The simulator provides very rudimentary support for generating 3D measurements. In this case the
   landmarks have a Z coordinate derived from their index in some deterministic random scheme. The
   robot makes the primary planar motion as in the 2D case overlayed with a random 6-DOF transform
   to make the motion full 3D. Note, that in the 3D case currently no odometry is supported.
*/
class TmgSimulator 
{
 public:
  //! Empty constructor
  TmgSimulator();
  
  //! Load a building. See \c load
  TmgSimulator(const char* configFilename);
  
  //! Desctructor
  ~TmgSimulator();  

  //! Reset to empty
  void clear();  

  //! Load a building from file \c configFilename.
  /*! \c configFilename is the name of the text file specifying simulation parameter,
      the different stories and how the robot shall move between stories. Each individual
      story is defined by a bitmap (see \c Story). The same bitmap can be used for several
      stories, so the simulator can efficiently simulate very large buildings.

      The configuration file consists of a sequence of lines with the following statements:
      \code
          # comment
          SCALE <length>
          Sets the scale used for the bitmaps. One pixel corresponds to a square of length
          \c length. The same scale is used for all bitmaps.
          
          STORY <nr> <scale><filename>
          Defines story number <nr> with bitmap <filename>

          ELEVATETO <nr>                 
          When reaching an elevator go to story number <nr>

          ASSUMECOVARIANCESCALE <scaleOdo> <scaleLandmark> 
          Give the SLAM algorith a covariance matrix that is scaled by \c scaleOdo (odometry) or
          \c scaleLandmark (landmark observations) relative to the true covariance used for computing
          artifical noise in the simulation.

          ODOMETRY <bias> <sigma> <radius>
          Noise in the odometry: A bias of \c bias (m/m) and Gaussian noise of
          \c sigma (m/sqrt(m)). The \c radius of the robot determines the factor between wheel noise
          and resulting orientation noise.

          SENSOR <fieldOfView> <maxDistance> <observationProb>
          Defines the landmarks sensor characteristic. A landmark is detected within a pie
          of \c +/-fieldOfView/2 (deg) and a range of \c maxDistance (m). It is observed with
          a probability of \c observationProb.

          SENSORNOISE <angleSigma0> <angleSigma1> <angleSigma2> <distSigma0> <distSigma1> <distSigma2>  <distBias>
          Landmark sensor noise. The angle is disturbed by a Gaussian with std deviation
          \c angleSigma0 + dist*angleSigma1 + dist^2*angleSigma2. (deg, deg/m, deg/m^2) 
          The distance is disturbed by a Gaussian
          with std dev \c distSigma0 + dist*distSigma1 + dist^2*distSigma2 (m, 1, 1/m). 
          Additional there is a bias of \c distBias*angle (m/deg) that is coupled with the angle. This sort of
          bias is one of the few errors that can create a really large orientation error during forward mapping.
          All other errors are pretty much reduced by observing several landmarks several times.

          CANSEETHROUGHELEVATOR
          Activates, that when the robot goes through an elevator it simultaneously sees landmarks in both stories.
          This makes a difference on how information is passed between stories. If the flag is activated different
          stories are connected by commonly observed landmarks otherwise only by odometry.

          RANDOMMOTION3D <rotation> <translation>
          The robot path defined by the simulation map is mostly 2D (except for passing between stories). To create
          a very basic 3D motion the 2D pose is perturbed by a random rotation of <rotation> (deg) and a random
          translation of <translation> before observing landmarks. This is equivalent to the camera jerking on the
          robot.

          STORYHEIGHT <height>
          Height of one story. Used for computing 3D observations that go across stories.
          

          ACTIVATENOISE <randomSeed>
          Switches artificial noise on with \c randomSeed to initialize the random number generator.
      \endcode
  */
  void load (const char* configFilename);

  //! Defines a single story of the simulated building.
  /*! An \c Story object can be referenced several times creating
      copies of the story. Thus it does not contain the story
      number. It is specified by a bitmap where the robot trajectory
      is marked red, wall black, landmarks blue and elevators green.
   */
  class Story 
    {
    public:
      Story ();
      
      //! Load a story (see \c load)
      Story (char* ppmFilename);      

      //! Load a story       
      /*! Pixel are converted into \c data. Landmarks are numbered consecutively. 
          \c ppmFilename must be an .ppm file.
      */
      void load (char* ppmFilename);      
      
      //! Filename from which the bitmap was loaded
      string ppmFilename;      

      //! Number of landmarks in this story
      int nrOfLandmarks;

      //! Number of path pixel
      int nrOfPath;      

      //! Image dimensions
      int width, height;

      //! Physical dimension of a pixel (m/pixel)
      double scale;  

      
      //! Special entries
      /*! \enum FREE a floor pixel (all other) */
      /*! \enum WALL a wall pixel (black). Block the sight for landmarks. */
      /*! \enum PATH a path pixel (red). The robot follows this pass. May have 90 deg crossings. */
      /*! \enum ELEVATOR an elevator pixel (green). When the robot reaches this position it goes to a different story.
                An elevator is also a waypoint setting a flag to execute a special operation (i.e. a snapshot).*/
      /*! \enum START the start pixel (cyan) Here the robot starts. It also serves as an elevator. */
      /*! \enum WAYPOINT. When the robot reaches this position \c the waypoint flag is set. */
      /*! \enum LANDMARK starting index for consecutively numbered landmarks. The storie's base landmark index must
          be added. */
      enum {FREE=0, WALL=1, PATH=2, WAYPOINT=3, ELEVATOR=4, START=5, LANDMARK=6};

      //! The actual image
      vector<int> data;

      int operator() (int x, int y) const {return data[x+width*y];}

      //! Returns, whether there is a double landmark at \c x,y
      /*! A double landmark is a landmark with another landmark as
          neighbor. True is only returned for the left/upper of the
          pair. */
      bool isDoubleLandmark (int x, int y) const;      

      //! Returns the one \c START
      /*! Orientation points to an adjacent \c PATH pixel. The purpose is to
          automatically determine a position for starting the robot. If no
          \c START or two \c START are found , \c x=y=orientation=-1 is returned.
      */
      void findStartPosition (int& x, int& y, int& orientation);      

      //! Returns whether this pixel is treated as a wall
      /*! A wall is a pixel with R, G, and B component all below 0x40. */
      bool isWallPixel (unsigned int pix);      
    };

  //! Stories of the building. 
  /*! The same \c Story object can be referenced in several stories. */
  vector<Story*> story;

  //! Base landmark id of different stories story nr \c storyNr
  /*! The value \c baseLandmarkId[storyNr] must be added to (p-Story::LANDMARK) for a landmark
      entry in \c story[storyNr] to get the final landmark id. This is necessary since the
      same \c Story data structure can be used for different -- identical -- stories.

      \c baseLandmarkId contains one entry more than \c story where the last entry consistently
      contains the total number of landmarks
  */
  vector<int> baseLandmarkId;

  //! List of different stories of the building
  /*! Each \c Story object is only once, so by freeing all objects from
    this list memory can be deallocated. */
  vector<Story*> originalStory;  

  //! How the robot moves from story to story.
  /*! Whenever the robot hits an \c ELEVATOR pixel, it is moved to the
      same \c x/\c y position on \c storyTransition[nextTransition]
      and \c nextTransition is incremented. This way the robot is
      guided through the stories, whereas in a single story it is
      guided by the red line.
  */
  vector<int> storyTransition;  

  //! Go to \c storyTransition[nextTransition] when hitting the next \c ELEVATOR
  /*! When \c nextTransition==storyTransition.size() this does not mean to immediately
      stop but still to go to the next elevator and stop there. So when hitting the
      elevator \c nextTransition is incremented to \c storyTransition.size()+1.
  */
  int nextTransition;  

  //! Robot position in pixel  (y pointing downward)
  int robotX, robotY;

  //! Robot position in stories
  int robotZ;

  //! Robot position before transition
  int robotOldZ;  

  //! orientation in 45 degree steps
  int robotOrientation;

  //! Whether the robot has hit (and used) a waypoint or elevator in the last step
  /*! This flag can be used to save a snapshot at regular intervals. */
  bool hasHitWaypoint;  


//  //! All landmark id's in the current story must be increased by \c baseLandmarkId
//  int baseLandmarkId;

  //! Size of one pixel
  double scale;  

  //! Height of one story
  /*! Used to compute 3D coordinates for observations going over several stories. */
  double storyHeight;    
  
  //! Steps since last elevator
  int nrOfStepsSinceLastElevator;  

  //! Whether the robot is still at the starting position
  bool isFirst;  

  //! Returns the bounding box (world coordinates) of story \c storyNr
  void boundingBox (int storyNr, double& loX, double& hiX, double& loY, double& hiY);
  

  //! Random number generator used for the simulation
  VmRandom random;

  //! The robots simulation characteristic
  class Robot 
  {
  public:
    //! Whether the robot is there
    bool isThere;
    
  
    //! Scale for the covariance passed to the algorithm
    /*! The covariance used to generate artificial noise is scaled by
      \c _assumeCovarianceScale* before it is passed to the estimation
      algorithm. The usual case is \c _assumeCovarianceScale*==1, so
      the estimation is based on the same covariance used for
      simulation.
      
      \c _assumeCovarianceScaleOdo is used for odometry, \c
      _assumeCovarianceScaleLandmark for landmark measurements.
    */
    double assumeCovarianceScaleOdo, assumeCovarianceScaleLandmark;
    
    //! Robot odometry characteristic
    /*! \c _odoBias is the odometry bias in angle/distance.  \c
      _odoSigma is the std deviation of the translational error per
      sqrt of distance moved. The corresponding std deviation of of
      the orientational error is \c _odoSigma/_rotationRadius. This
      would be the case, if the robot had wheels with a distance of \c
      2* _rotationRadius.  Consequently the distance move consists of
      the translational distance plus the orientational distance
      multiplied by \c _rotationRadius.
    */
    double odoBias, odoSigma, rotationRadius;
    

    //! Robot landmark sensor range characteristic
    /*! The robot can observe landmarks, if they are within an angle
      of \c +/-_fieldOfView/2 and at a distance of at most \c
      _maxDistance.  It certainly cannot observe landmarks through
      walls.
    */
    double fieldOfView, maxDistance;

    //! Robot landmark sensor noise characteristics
    /*! \c The uncertainty of the measurement of a landmark of
      distance \c d is assumed to be \c
      (_sigma[0][0]+_sigma[0][1]*d+_sigma[0][2]*d^2) for the
      angular and \c
      (_sigma[1][0]+_sigma[1][1]*d+_sigma[1][2]*d^2) for the
      distance error.q
      
      \c _observationProb determines the chance for observing
      a landmark that is in range and unoccluded by walls.
    */
    double landmarkSigma[2][3], observationProb;

    //! Bias in the landmark observation
    /*! \c _bias*angle is added to the distance, if angle is the angle
      relative to the center of the field of view. \c _landmarkBias is
      not directly defined by the config file but drawn from a
      gaussian distribution defined by \c _landmarkBiasSigma. This
      allows to make a statistic over several different bias values.
    */
    double landmarkBias;

    //! The bias is constant for a single robot but random over several robots
    /*! This is achieved by drawing \c _landmarkBias from a 0 mean
      gaussian distribution with stddev \c _landmarkBiasSigma upon
      initialisation.
    */
    double landmarkBiasSigma;


    //! Whether the robot see both stories in the same view when going through an elevator
    bool canSeeThroughElevator;    

    //! 3D translation and rotation noise added to the 2D robot pose to get a 3D observer position
    double random3DTranslation, random3DRotation;    


    //! Whether sensor noise is activated
    bool isNoiseActive;
    
    //! Constructs an invalid robot (\c _isThere==false)
    Robot();

    //! Set \c _isThere to false
    void clear();
        
    //! Returns the noise parameter at a distance of \c dist.
    void sigma (double dist, double& angSig, double& distSig);

    //! Starts a new run and draws a \c landmarkBias with std. dev. \c landmarkBiasSigma
    void init (VmRandom& random);    

    //! Disturbs the origional \c pos landmark measurement by gaussian noise and bias
    /*! The noise and bias depend on the parameters defined. \c pos is updated
        with the distorted measurement and \c posCov is the covariance passed
        to the SLAM algorithm.
    */
    void disturbLandmarkMeasurement (VmRandom& random, VmVector2& pos, VmMatrix2x2& posCov);

    //! Same for a 3D measurement
    /*! Again the noise consists of a distance and angle component. No landmarkBias is supported.
    */
    void disturbLandmarkMeasurement3D (VmRandom& random, VmVector3& pos, VmMatrix3x3& posCov);

    //! Disturbs the original \c pos odometry measurement by gaussian noise and bias
    /*! The noise and bias depend on the parameters defined. \c pos is updated
        with the distorted measurement and \c posCov is the covariance passed
        to the SLAM algorithm.
    */
    void disturbOdometryMeasurement (VmRandom& random, VmVector3& pos, VmMatrix3x3& posCov);    

    //! Applies \c random3DRotation and \c random3DTranslation to \c pose
    void disturbPose (VmRandom& random, VmMatrix4x4& pose) const;    

    //! Computes a reasonable heuristic covariance for a relative movement of \c d
    /*! \c odoSigma specifies the std deviation of the error of a
        single wheel when moving a distance of 1. \c radius specifies
        the radius of the robot, which should be set as half the
        maximum distance between two wheels if the wheels are not
        placed symmetrically. This parameter governs the scaling
        between orientation and translation errors, since the
        odometric orientation error is generated by the error in the
        wheel measurements.

        The result is returned in \c dCov. Take a look at the source code
        to see a derivation of the formula.
    */
    void stdCovariance (const VmVector3& d, VmMatrix3x3& dCov,
                        double radius, double odoSigma);
  protected:
    //! Rotates \c pose by angle in the \c i-j plane.
    static void rotate (VmMatrix4x4& pose, int i, int j, double angle);    
  };


  //! The robot's simulation parameter
  Robot robot;  

  //! Returns the true robot pose
  void trueRobotPose (VmVector3& p);
  
  //! Returns the true robot pose (without RANDOMMOTION3D camera jerk) in 3D
  /*! It returns a homogenous 4*4 matrix. */
  void trueRobotPose3D (VmMatrix4x4& pose);  

  //! True 3D position of the landmark at \c x, y in story \c story
  /*! The routine both incorporates \c scale and \c storyHeight but
      also a deterministic pseudo-random number derived from the index
      of the landmark there that determines the z coordinate between
      stories.
  */
  void landmarkPosition3D (VmVector3& p, int x, int y, int storyNr);  

  //! Move one step
  /*! The robot follows the red (\c PATH) line. It first tries to go
      forward, if there is no \c PATH pixel in that direction it
      attemps a (-45, +45, -90, +90deg) rotation. This way the red
      path can cross itself orthogonally without the robot loosing
      track. When the robot hits an \c ELEVATOR, it is moved to story
      number \c storyTransition[nextTransition] and \c nextTransition
      is incremented. This way the robot is
      guided through the stories, whereas in a single story it is
      guided by the red line.

      In \c odometry the \c (deltaX, deltaY, deltaTheta) movement is
      returned. It is perturbed by noise and bias (to be implemented
      yet). In \c odometryCov the corresponding covariance is returned.
  */
  void step (VmVector3& odometry, VmMatrix3x3& odometryCov);

  //! Observation of a landmark with known identity
  class Observation 
    {
    public:
      //! Landmark id
      int id;      
      //! Position relative to the robot
      VmVector2 pos;
      //! Position measurement covariance
      VmMatrix2x2 posCov;      

      Observation ()
        :id(-1)
        {
          vmZero (pos);
          vmZero (posCov);
        }
      

      Observation (int id, const VmVector2& pos, const VmMatrix2x2& posCov)
        :id(id)
        {
          vmCopy (pos, this->pos);
          vmCopy (posCov, this->posCov);          
        }
    };

  typedef Observation Observation2D;  

  class Observation3D 
    {
    public:
      //! Landmark id
      int id;      
      //! Position relative to the robot
      VmVector3 pos;
      //! Position measurement covariance
      VmMatrix3x3 posCov;      

      Observation3D ()
        :id(-1)
        {
          vmZero (pos);
          vmZero (posCov);
        }
      

      Observation3D (int id, const VmVector3& pos, const VmMatrix3x3& posCov)
        :id(id)
        {
          vmCopy (pos, this->pos);
          vmCopy (posCov, this->posCov);          
        }
    };
    
      

  //! Returns measurements for all landmarks visible from the current robot pose
  /*! A landmark can be observed if it is on the same story within \c maxRange
      distance and in an angle of \c +/-halfFieldOfView with respect to the
      robot's front. Measurements are perturbed by noise and bias.
  */
  void observe (vector<Observation>& observation);

  //! Returns 3D measurements for all landmarks visible from the current robot pose
  /*! Each individual measurement is disturbed by specified noise. Additionally the
      pose from which the observations are taken is corrupted by noise. As a special
      rule this does not happen for the first pose.
   */
  void observe3D (vector<Observation3D>& observation);  

  //! Returns the range of landmarks \c [from..to-1] in the current story
  void landmarksInStory (int& from, int& to) const;  

  //! Returns the area (true coordinates) covered by the simulated building
  void area (double& loX, double& hiX, double& loY, double& hiY) const;

  //! Returns the story number of which \c landmarkId belongs
  int storyOfLandmark (int landmarkId) const;
  
    

  //! Whether the simulation is finished
  /*! It is finished when the robot reaches an the last \c ELEVATOR (i.e.
    \c nextTransition==storyTransition.size().
  */
  bool isFinished () const
    {return nextTransition>(int) storyTransition.size();}  

  //! Returns the pixel on which the robot is currently
  int pixel () const;  

  //! Returns pixel x, y in the current story
  int pixel (int x, int y) const;  

  //! Returns the robot's neighbor pixel in direction \c
  //! relOrientation relative to it's own orientation.
  int pixelRel (int relOrientation) const;  

  //! Returns whether \c pixelRel (relOrientation) is \c PATH or \c ELEVATOR
  bool canMove (int relOrientation) const;

  //! Moves the robot directly one step ahead.
  void move();  

  //! Sine and Cosine of the robot's orientation
  void sinCos (double& s, double& c);

  //! Computes \c baseLandmarkId and \c myNrOfLandmarks
  void computeBaseLandmarkIds ();  

  int computeBaseLandmarkId (int storyNr) const;  

  //! Total number of landmarks in all stories
  int nrOfLandmarks () const 
  {
    if (!baseLandmarkId.empty()) return baseLandmarkId.back();
    else return 0;
  }
  

  //! Returns whether there is no \c WALL in the line between the robot and \c x,y
  bool isLineOfSightFree (int x, int y);  

  //! Auxiliary class storing a position inside a story (pixel)
  class Position 
    {
    public:
      int x, y;
      Position (int x, int y):x(x), y(y){};
      Position ():x(0), y(0){};
      double distance (const Position& p2) const
        {return sqrt((double)((x-p2.x)*(x-p2.x) + (y-p2.y)*(y-p2.y)));}      
    };  
      
 protected:
  //! Internal subroutine for \c observe
  /*! Adds observations from \c fromX, \c fromY, \c fromStory to \c observation.
   */
  void internalAddObservations (vector<Observation>& observation, int fromX, int fromY, int fromStory);  

  //! Internal subroutine for \c observe3D
  /*! Adds observations from \c fromX, \c fromY, \c fromStory to \c observation.
      Use \c pose as observer2World transform.
   */
  void internalAddObservations3D (vector<Observation3D>& observation, int fromX, int fromY, int fromStory, const VmMatrix4x4& pose);  

  //! Computes the range of cells that may potentially be seend from \c x, y
  /*! Clipped to the area of \c story[fromStory]. */
  void visibleRectangle (int& loX, int& hiX, int& loY, int& hiY, int x, int y, int fromStory);  

};
  
#endif
