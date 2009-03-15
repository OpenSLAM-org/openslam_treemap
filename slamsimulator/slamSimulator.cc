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
/*!\file slamSimulator.cc 
   \brief Implementation of \c SlamSimulator
   \author Udo Frese

   Contains the implementation of class \c
   SlamSimulator that can simulate robot motion and landmark perception
   in a multi-story building specified by a bitmap.
*/

#include "slamSimulator.h"
#include <stdexcept>
#include <math.h>

//! Used to convert directions [0..7] into vectors \c (dX,dY)
int dX[8] = {1, 1, 0, -1, -1, -1, 0, 1};
//! Used to convert directions [0..7] into vectors \c (dX,dY)
int dY[8] = {0,-1,-1, -1,  0,  1, 1, 1};

SlamSimulator::SlamSimulator()
  :story(), originalStory(), storyTransition(), nextTransition(0), 
   robotX(0), robotY(0), robotZ(0), robotOrientation(0), isFirst (true),
   random()
{}

  
SlamSimulator::SlamSimulator(const char* specification)
  :story(), originalStory(), storyTransition(), nextTransition(0), 
   robotX(0), robotY(0), robotZ(0), robotOrientation(0), isFirst (true),
   random()
{  
  load (specification);  
}


void SlamSimulator::load (const char* configFilename)
{
  FILE* f = fopen (configFilename, "r");
  if (f==NULL) throw runtime_error ("Could not open "+string(configFilename));
  clear();  
  try {
    robot.assumeCovarianceScaleOdo = robot.assumeCovarianceScaleLandmark = 1;
    robot.isNoiseActive = false;    
    scale = 1;    
    while (!feof(f)) {
      char buffer[1010];
      char comment[1010];
      char fname[1000];
      int storyNr, hv;
      
      buffer[0]='\0';
      fgets (buffer, 999, f);
      if (buffer[0]=='\0') break; //! EOF
      int n = strlen(buffer);
      if (n>0 && buffer[n-1]=='\n') n--;
      strcpy (buffer+n, " #");
      
      if (sscanf (buffer, " %[#]", comment)==1) {
        // Empty or comment line
      }
      else if (sscanf (buffer, "STORY %d %s %[#]",
                       &storyNr, fname, comment)==3) {
        // Load a story from fname
        if (storyNr>=(int) story.size()) story.resize (storyNr+1, NULL);
        int i;
        for (i=strlen(configFilename);i>=0 && configFilename[i]!='/'; i--);
        i++;
        for (int j=strlen(fname); j>=0; j--) fname[i+j]=fname[j];
        for (int j=0; j<i; j++) fname[j] = configFilename[j];        
        i=0;
        while (i<(int) originalStory.size() && originalStory[i]->ppmFilename!=fname) i++;
        if (i==(int) originalStory.size()) originalStory.push_back (new Story(fname));
        story[storyNr] = originalStory[i];          
      }      
      else if (sscanf (buffer, "ELEVATETO %d %[#]",
                       &storyNr, comment)==2) {
        // add an \c storyTransition entry
        if (storyNr<0 || (int) story.size()<=storyNr) throw runtime_error ("Illegal story number.");
        storyTransition.push_back (storyNr);        
      }      
      else if (sscanf (buffer, "ASSUMECOVARIANCESCALE %lf %lf %[#]", 
                       &robot.assumeCovarianceScaleOdo, 
                       &robot.assumeCovarianceScaleLandmark, 
                       comment)==3) {
        // Define covariance scale
      }        
      else if (sscanf (buffer, "ODOMETRY %lf %lf %lf %[#]", 
                       &robot.odoBias, &robot.odoSigma, 
                       &robot.rotationRadius, comment)==4) {
        robot.odoBias = vmRad(robot.odoBias);
        // Define odometry parameter
      }
      else if (sscanf (buffer, "SENSOR %lf %lf %lf %[#]", 
                       &robot.fieldOfView, &robot.maxDistance, &robot.observationProb, comment)==4) {
        // Define sensor range characteristics
        robot.fieldOfView = vmRad(robot.fieldOfView);
      }
      else if (sscanf (buffer, "SENSORNOISE %lf %lf %lf %lf %lf %lf %lf %[#]", 
                       &robot.landmarkSigma[0][0], &robot.landmarkSigma[0][1], &robot.landmarkSigma[0][2], 
                       &robot.landmarkSigma[1][0], &robot.landmarkSigma[1][1], &robot.landmarkSigma[1][2], 
                       &robot.landmarkBiasSigma, 
                       comment)==8) {
        // Define sensor noise characteristics
        robot.landmarkSigma[0][0] = vmRad(robot.landmarkSigma[0][0]);
        robot.landmarkSigma[0][1] = vmRad(robot.landmarkSigma[0][1]);
        robot.landmarkSigma[0][2] = vmRad(robot.landmarkSigma[0][2]);            
        robot.landmarkBiasSigma *= 180/M_PI; // convert from 1/deg to 1/rad
      }
      else if (sscanf (buffer, "CANSEETHROUGHELEVATOR %[#]", comment)==1) {
        robot.canSeeThroughElevator = true;        
      }      
      else if (sscanf (buffer, "RANDOMMOTION3D %lf %lf %[#]", &robot.random3DRotation, &robot.random3DTranslation, comment)==3) {
        // Define additional jerking motion of the camera to define a 3D trajectory
        robot.random3DRotation = vmRad(robot.random3DRotation);
      }      
      else if (sscanf (buffer, "STORYHEIGHT %lf %[#]", &storyHeight, comment)==2) {
        //Set storyHeight
      }      
      else if (sscanf (buffer, "ACTIVATENOISE %d %[#]", &hv, comment)==2) {
        // Switch the sensor noise on
        random.create (hv);        
        robot.isNoiseActive = true;        
      }
      else if (sscanf (buffer, "SCALE %lf %[#]", &scale, comment)==2) {
        // Set scale
      }      
      else throw runtime_error ("Could not parse "+string(buffer));
    }
    // Set robot to initial position
    nextTransition = 0;    
    robotZ = 0;    
    robotOldZ = 0;    
    nrOfStepsSinceLastElevator = 0;
    if (story.empty()) throw runtime_error ("No story defined");    
    story[0]->findStartPosition (robotX, robotY, robotOrientation);
    if (robotX<0) throw runtime_error ("Could not find starting position.");    
    computeBaseLandmarkIds ();
    robot.init (random);
    hasHitWaypoint = false;    
    isFirst = true;    
  }  
  catch (...) {
    clear();
    if (f!=NULL) fclose(f);
    throw;
  }
  fclose(f);  
}


void SlamSimulator::clear()
{  
  for (int i=0; i<(int) originalStory.size(); i++)
    delete originalStory[i];
  originalStory.clear();
  story.clear();
  storyTransition.clear();  
  nextTransition=0;
  robotX=0;
  robotY=0;
  robotZ=0;
  robotOldZ = 0;  
  robotOrientation=0;
  robot.clear();
}


SlamSimulator::~SlamSimulator()
{
  clear();  
}


int SlamSimulator::pixel () const
{
  return (*story[robotZ]) (robotX, robotY);
}

int SlamSimulator::pixel (int x, int y) const
{
  return (*story[robotZ]) (x, y);
}


int SlamSimulator::pixelRel (int relOrientation) const
{
  int ori = (robotOrientation+relOrientation) & 7;
  return (*story[robotZ]) (robotX+dX[ori], robotY+dY[ori]);  
}


bool SlamSimulator::canMove (int relOrientation) const
{
  int pr = pixelRel (relOrientation);
  return pr==Story::PATH || pr==Story::ELEVATOR || pr==Story::START || pr==Story::WAYPOINT;
}


void SlamSimulator::move()
{
  robotX += dX[robotOrientation&7];
  robotY += dY[robotOrientation&7];
  isFirst = false;  
}


void SlamSimulator::boundingBox (int storyNr, double& loX, double& hiX, double& loY, double& hiY)
{
  if (0<=storyNr && storyNr<(int) story.size()) {
    loX = hiX = loY = hiY = 0;
    return;
  }
  loX = 0;
  hiX =  scale*story[storyNr]->width;  
  loY = -scale*story[storyNr]->height;
  hiY = 0;  
}


void SlamSimulator::step (VmVector3& odometry, VmMatrix3x3& odometryCov)
{
  int oldOrientation = robotOrientation;
  int oldX = robotX;
  int oldY = robotY;  
  double c, s;
  sinCos (s, c);  

  if (canMove (0)) {}  
  else if (canMove(1)) robotOrientation +=1;
  else if (canMove(-1)) robotOrientation -=1;
  else if (canMove(2)) robotOrientation +=2;
  else if (canMove(-2)) robotOrientation -=2;
  else {
      char txt[1000];
      sprintf(txt, "Path ended suddenly (%d/%d/%d)", robotX, robotY, robotZ);      
      throw runtime_error (string(txt));
  }
  move();    
  nrOfStepsSinceLastElevator++;
  assert (nrOfStepsSinceLastElevator<=story[robotZ]->nrOfPath+1);  
  robotOldZ = robotZ;    
  if (pixel()==Story::ELEVATOR || pixel()==Story::START) {
    if (nextTransition<(int) storyTransition.size())
      robotZ = storyTransition[nextTransition];
    if (!pixel()==Story::ELEVATOR || pixel()==Story::START) {
      char txt[1000];
      sprintf(txt, "ELEVATOR's not on top of each other (%d/%d/%d-->%d)",
              robotX, robotY, robotOldZ, robotZ);      
      throw runtime_error (string(txt));
    }    
    char* txt[] = {"E", "NE", "N", "NW", "W", "SW", "S", "SE"};                
    printf ("ELEVATOR transition %d/%d (%d, %d, %s) %d-->%d\n",
            nextTransition, storyTransition.size(),
            robotX, robotY, txt[robotOrientation&0x7], robotOldZ, robotZ);    
    nrOfStepsSinceLastElevator = 0;
    nextTransition++;
    hasHitWaypoint = true;    
  }  
  else if (pixel()==Story::WAYPOINT) hasHitWaypoint = true;
  else hasHitWaypoint = false;  

  // we use - for y, because pixel coordinates are negative of world coordinates
  odometry[0] =  scale* ( (robotX-oldX)*c + -(robotY-oldY)*s);
  odometry[1] =  scale* (-(robotX-oldX)*s + -(robotY-oldY)*c);  
  odometry[2] = (robotOrientation-oldOrientation)*M_PI/4;

  robot.disturbOdometryMeasurement (random, odometry, odometryCov);  
}


void SlamSimulator::sinCos (double& s, double& c)
{
  c = cos(robotOrientation*M_PI/4);
  s = sin(robotOrientation*M_PI/4);  
}



void SlamSimulator::visibleRectangle (int& loX, int& hiX, int& loY, int& hiY, int x, int y, int fromStory)
{
  if (fromStory<0 || (int) story.size()<=fromStory) {
    loX= loY = 0;
    hiX= hiY = -1;
  }  
  Story* st = story[fromStory];  
  loX = (int) floor (x-robot.maxDistance/scale);   
  if (loX<0) loX = 0;  
  hiX = (int) ceil  (x+robot.maxDistance/scale);
  if (hiX>=st->width) hiX = st->width-1;  
  loY = (int) floor (y-robot.maxDistance/scale);
  if (loY<0) loY = 0;  
  hiY = (int) ceil  (y+robot.maxDistance/scale);  
  if (hiY>=st->height) hiY = st->height-1;
}



void SlamSimulator::internalAddObservations2D (vector<Observation2D>& observation, int fromX, int fromY, int fromStory)
{
  if (fromStory<0 || (int) story.size()<=fromStory) return;  
  Story* st = story[fromStory];  
  int loX, loY, hiX, hiY;  
  visibleRectangle (loX, hiX, loY, hiY, fromX, fromY, fromStory);  

  for (int x=loX;x<=hiX;x++)
    for (int y=loY;y<=hiY;y++) {
      int p = (*st) (x, y);
      if (p>=Story::LANDMARK && isLineOfSightFree (x, y)) { // found a visible landmark
        VmVector2 pos;
        double c,s;
        sinCos (s, c);
        //! - because world coordinates are - pixel y
        pos[0] = scale*( (x-robotX)*c + -(y-robotY)*s);
        pos[1] = scale*(-(x-robotX)*s + -(y-robotY)*c);
        
        double relativeAngle = atan2 (pos[1], pos[0]);
        double dist = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
        if (dist>robot.maxDistance || fabs(relativeAngle)>robot.fieldOfView/2) continue;        
        VmMatrix2x2 posCov;
        robot.disturbLandmarkMeasurement (random, pos, posCov);
        int id = baseLandmarkId[fromStory]+p-Story::LANDMARK;
        assert (0<=id && id<nrOfLandmarks());        
        observation.push_back (Observation2D (id, pos, posCov));
      }
    }
}


void SlamSimulator::observe2D (vector<Observation2D>& observation)
{
  observation.clear();
  internalAddObservations2D (observation, robotX, robotY, robotZ);
  if (robot.canSeeThroughElevator && robotOldZ!=robotZ) 
    internalAddObservations2D (observation, robotX, robotY, robotOldZ);  
}


void SlamSimulator::internalAddObservations3D (vector<Observation3D>& observation, int fromX, int fromY, int fromStory, const VmMatrix4x4& pose)
{
  if (fromStory<0 || (int) story.size()<=fromStory) return;  
  Story* st = story[fromStory];  
  int loX, loY, hiX, hiY;  
  visibleRectangle (loX, hiX, loY, hiY, fromX, fromY, fromStory);  

  for (int x=loX;x<=hiX;x++)
    for (int y=loY;y<=hiY;y++) {
      int p = (*st) (x, y);
      if (p>=Story::LANDMARK && isLineOfSightFree (x, y)) { // found a visible landmark
        VmVector3 posWorld, posRobot;
        landmarkPosition3D (posWorld, x, y, fromStory);
        vmApplyInverseTransformToPoint (posRobot, pose, posWorld);    
        VmMatrix3x3 posCov;
        robot.disturbLandmarkMeasurement3D (random, posRobot, posCov);        

        int id = baseLandmarkId[fromStory]+p-Story::LANDMARK;
        assert (0<=id && id<nrOfLandmarks());        
        observation.push_back (Observation3D (id, posRobot, posCov));
      }
    }
}


void SlamSimulator::observe3D (vector<Observation3D>& observation)
{
  observation.clear();
  VmMatrix4x4 pose;
  trueRobotPose3D (pose);
  if (!isFirst) robot.disturbPose (random, pose);

  internalAddObservations3D (observation, robotX, robotY, robotZ, pose);
  if (robot.canSeeThroughElevator && robotOldZ!=robotZ) {
    internalAddObservations3D (observation, robotX, robotY, robotOldZ, pose);  
  }  
}


bool SlamSimulator::isLineOfSightFree (int x, int y)
{  
  int dX = robotX-x;
  int dY = robotY-y;
  // maintain 'equation==0' as close to as possible
  int equationBound = max(abs(dX), abs(dY));  
  while ((x!=robotX) || (y!=robotY)) {
    int equation = (x-robotX)*dY - (y-robotY)*dX;  // not very efficient but who cares
    assert (abs(equation) <= equationBound);
    if (pixel(x, y)==Story::WALL) return false;
    if (equation>0) {
      if (dX>=0 && dY>0) y++;
      else if (dX>0 && dY<=0) x++;
      else if (dX<0 && dY>=0) x--;
      else if (dX<=0 && dY<0) y--;
    }
    else {
      if (dX>0 && dY>=0) x++;
      else if (dX>=0 && dY<0) y--;
      else if (dX<=0 && dY>0) y++;
      else if (dX<0 && dY<=0) x--;
    }
  }
  if (pixel(x, y)==Story::WALL) return false;
  return true;
}


void SlamSimulator::computeBaseLandmarkIds ()
{
  baseLandmarkId.clear();
  int ctr = 0;
  for (int i=0;i<(int) story.size();i++) {
    baseLandmarkId.push_back (ctr);
    ctr += story[i]->nrOfLandmarks;
  }
  baseLandmarkId.push_back (ctr);  
}


void SlamSimulator::landmarksInStory (int& from, int& to) const
{
  from = baseLandmarkId[robotZ];
  to   = baseLandmarkId[robotZ+1];  
}


//************* SlamSimulator::Story


SlamSimulator::Story::Story ()
  :ppmFilename(), nrOfLandmarks (0), width(0), height(0), data()
{}


SlamSimulator::Story::Story (char* ppmFilename)
  :ppmFilename(), nrOfLandmarks (0), width(0), height(0), data()
{
  load (ppmFilename);
}


void SlamSimulator::Story::findStartPosition (int& x, int& y, int& orientation)
{
  int ctr = 0;  
  for (int yy=0; yy<height; yy++)
    for (int xx=0; xx<width; xx++)
      if ((*this)(xx, yy)==START) {
        ctr++;        
        x = xx;
        y = yy;
        for (int i=0; i<8; i++)
          if ((*this)(xx+dX[i], yy+dY[i])==PATH) {
            orientation = i;
            break;
          }
      }
  if (ctr!=1) x = y = orientation = -1;
}


bool SlamSimulator::Story::isWallPixel (unsigned int pix)
{
  return (pix & 0xff) <= 0x40 && (pix & 0xff00) <= 0x4000 && (pix & 0xff0000) <= 0x400000;  
}



void SlamSimulator::Story::load (char* ppmFilename)
{
  data.clear();  
  this->ppmFilename = ppmFilename;
 
  FILE* f;
  int dummy, buflen;
  char buf[2000];

  nrOfLandmarks = nrOfPath = 0;  
  buflen=2000;
  f = fopen (ppmFilename, "r");
  if (f==NULL) throw runtime_error ("Could not open .ppm file "+string(ppmFilename));

  fgets( buf, buflen, f );
  if( strcmp( buf, "P6\n" )!=0 ) throw runtime_error ("Illegal ppm format");

  do { 
    fgets( buf, buflen, f );
  } while( buf[0]=='#' );
  if (sscanf( buf, "%d %d", &width, &height )<2)
    throw runtime_error ("Could not find image dimensions");
  data.reserve (width*height);  
 
  do { fgets( buf, buflen, f ); } while( buf[0]=='#' );
  sscanf( buf,"%d", &dummy );
    
  unsigned char* line = new unsigned char[3*width];
  for (int y=0; y<height; y++) {
    fread (line, 3, width, f);
    for (int x=0; x<width; x++) {
      unsigned char* p = line+3*x;
      int color = (p[0]<<16) + (p[1]<<8) + (p[2]<<0);
      if (isWallPixel(color)) data.push_back (WALL);
      else if (color==0xff0000) { // black or dark grey
        data.push_back (PATH);
        nrOfPath++;
      }      
      else if (color==0x00ff00) data.push_back (ELEVATOR); // green
      else if (color==0xffff00) data.push_back (WAYPOINT); // yellow
      else if (color==0x00ffff) data.push_back (START);    // cyan
      else if (color==0x0000ff) { // blue
        data.push_back (LANDMARK+nrOfLandmarks);
        nrOfLandmarks++;
      }
      else data.push_back (FREE);      
    }
  }
  fclose(f);
  delete[] line;
}


bool SlamSimulator::Story::isDoubleLandmark (int x, int y) const
{
  if ((*this)(x,y)<LANDMARK) return false;
  if (x+1<width  && (*this)(x+1,y)>=LANDMARK) return true;
  if (y+1<height && (*this)(x,y+1)>=LANDMARK) return true;
  return false;
}


// *************** SlamSimulator::Robot

SlamSimulator::Robot::Robot()
        :isThere(false), 
         assumeCovarianceScaleOdo(1), assumeCovarianceScaleLandmark(1),
         odoSigma(0), rotationRadius(0), fieldOfView(0), maxDistance(0),
         canSeeThroughElevator (false), 
         random3DTranslation(0), random3DRotation(0), isNoiseActive(false)
{}


void SlamSimulator::Robot::clear()
{
    isThere = false;
}


void SlamSimulator::Robot::init (VmRandom& random)
{
  if (isNoiseActive) landmarkBias = landmarkBiasSigma*random.gauss();
  else landmarkBias = 0;  
}


void SlamSimulator::Robot::sigma (double dist, double& angSig, double& distSig)
{
    angSig  = landmarkSigma[0][0] + dist*landmarkSigma[0][1] + dist*dist*landmarkSigma[0][2];
    distSig = landmarkSigma[1][0] + dist*landmarkSigma[1][1] + dist*dist*landmarkSigma[1][2];
}


void SlamSimulator::Robot::disturbLandmarkMeasurement (VmRandom& random, VmVector2& pos, VmMatrix2x2& posCov)
{ 
  double dist  = sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
  double angle = atan2(pos[1], pos[0]);
  
  double angSigma, distSigma;
  sigma (dist, angSigma, distSigma);
    
  VmMatrix2x2 jac =   // Jacobian of polar->cartesian mapping (dist, angle) input
    {{cos(angle), -dist*sin(angle)},
     {sin(angle), dist*cos(angle)}};
  VmMatrix2x2 cov =
    {{distSigma*distSigma, 0},
     {0, angSigma*angSigma}};
  VmMatrix2x2 dCov;
  vmJCKt (dCov, jac, cov, jac);
  vmScale (posCov, dCov, assumeCovarianceScaleLandmark);
  if (isNoiseActive) {
    VmVector2 dErr;
    random.gauss (dErr, dCov);
    dist += landmarkBias*angle;
    pos[0] = dist*cos(angle) + dErr[0];
    pos[1] = dist*sin(angle) + dErr[1];
  }
  else {
    pos[0] = dist*cos(angle);
    pos[1] = dist*sin(angle);
  }
}


void SlamSimulator::Robot::disturbLandmarkMeasurement3D (VmRandom& random, VmVector3& pos, VmMatrix3x3& posCov)
{
  double dist  = vmLength(pos);  
  double angSigma, distSigma;
  sigma (dist, angSigma, distSigma);

  double effLateralAngleSigma = dist*angSigma; // lateral effect [m] of the angular error
  VmMatrix3x3 dCov;
  vmZero (dCov);
  dCov[0][0] = dCov[1][1] = dCov[2][2] = effLateralAngleSigma*effLateralAngleSigma;
  double cov = distSigma*distSigma - effLateralAngleSigma*effLateralAngleSigma; 
  // may be negative since it cancels out that \c effLateralAngleSigma is applied in all
  // not only in lateral directions. The overall matrix is always SPD.
  vmJCKtAdd (dCov, pos, cov/(dist*dist), pos);  

  assert (landmarkBias==0);  

  vmScale (posCov, dCov, assumeCovarianceScaleLandmark);
  if (isNoiseActive) {
    VmVector3 dErr;
    random.gauss (dErr, dCov);
    vmAdd (pos, dErr);    
  }
}


void SlamSimulator::Robot::stdCovariance (const VmVector3& d, VmMatrix3x3& dCov,
                                         double radius, double odoSigma)
{
    double dO = d[2]*radius;
    double dX = d[0];
    double dY = d[1];
    double dist = sqrt(dX*dX + dY*dY + dO*dO);
// dist is approximately the average distance a wheel has travelled
    double tSigma2 = odoSigma*odoSigma*dist;
// which results in a translational error of sqrt(tSigma2)
    double oSigma2 = tSigma2 / (radius*radius);
// and an orientational error of sqrt(oSigma2)

// The movement is parametrized as x = lambda*dX/dist, y = lambda*dY/dist.
// lambda=[0..dist].
// We assume, that at every moment of the trajectory, the wheels have make
// a translational error with stddev 1 (later scaled to 'odoSigma').
// This results directly in a translational covariance contribution of
// {{1,0,0},{0,1,0},{0,0,0}}. Further it generates an orientation error
// of 1/radius which in turn will generate a correlated orientation/position 
// error as the robot continues moving.
// This error is: 
//      x(lambda):={-dY*(dist-lambda), dX*(dist-lambda), 1 }/radius. 
// The the contribution of the error made in a specific moment
// lambda is to the covariance is:
//      C(lambda):=x(lambda)x(lambda)^T + {{1,0,0},{0,1,0},{0,0,0}}
// The final result is derived by integrating C(lambda), lambda=[0..dist]
// and scaling by odoSigma:
    vmSet (dCov, 
           oSigma2* dY*dY/3+tSigma2, oSigma2*-dX*dY/3        , oSigma2* -dY/2,
           oSigma2*-dX*dY/3        , oSigma2* dX*dX/3+tSigma2, oSigma2*  dX/2,
           oSigma2*   -dY/2        , oSigma2*    dX/2        , oSigma2);
}


void SlamSimulator::Robot::disturbOdometryMeasurement (VmRandom& random, VmVector3& pos, VmMatrix3x3& posCov)
{ 
  double dO = pos[2]*rotationRadius;
  double dX = pos[0];
  double dY = pos[1];
  double dist = sqrt(dX*dX + dY*dY + dO*dO);
  double refLength = 1; 
  double tSigma2 = odoSigma*odoSigma*dist;
  double oSigma2 = tSigma2 / (rotationRadius*rotationRadius);
  double oBias = odoBias*dist;
  VmMatrix3x3 dCov;  
  vmSet (dCov, 
         oSigma2* dY*dY/3+tSigma2, oSigma2*-dX*dY/3        , oSigma2* -dY/2,
         oSigma2*-dX*dY/3        , oSigma2* dX*dX/3+tSigma2, oSigma2*  dX/2,
         oSigma2*   -dY/2        , oSigma2*    dX/2        , oSigma2);
  vmScale (posCov, dCov, assumeCovarianceScaleOdo);
  posCov[2][2] += dist*refLength*oBias*oBias;
  if (isNoiseActive) {
    VmVector3 noise;
    random.gauss (noise, dCov);
    noise[2] += oBias;    
    vmAdd (pos, noise);    
  }
}


void SlamSimulator::Robot::rotate (VmMatrix4x4& pose, int i, int j, double angle)
{
  double c = cos(angle), s = sin(angle);
  for (int k=0; k<3; k++) {
    double pKI = pose[k][i], pKJ=pose[k][j];    
    pose[k][i] = c*pKI - s*pKJ;
    pose[k][j] = s*pKI + c*pKJ;
  }  
}


void SlamSimulator::Robot::disturbPose (VmRandom& random, VmMatrix4x4& pose) const
{
  if (isNoiseActive) {    
    pose[0][3] += random3DTranslation*random.gauss();
    pose[1][3] += random3DTranslation*random.gauss();
    pose[2][3] += random3DTranslation*random.gauss();  
    rotate (pose, 1, 2, random3DRotation*random.gauss());
    rotate (pose, 2, 0, random3DRotation*random.gauss());
    rotate (pose, 0, 1, random3DRotation*random.gauss());
  }  
}


void SlamSimulator::area (double& loX, double& hiX, double& loY, double& hiY) const
{
  int wMax = 0;
  int hMax = 0;
  for (int i=0; i<(int) story.size(); i++) {
    if (story[i]->width  > wMax) wMax = story[i]->width;
    if (story[i]->height > hMax) hMax = story[i]->height;
  }
  loX = 0;
  hiX = wMax*scale;
  loY = -hMax*scale;
  hiY = 0;
}


void SlamSimulator::trueRobotPose (VmVector3& p)
{
  p[0] = scale*robotX;
  p[1] = -scale*robotY;
  p[2] = M_PI/4*robotOrientation;
}

void SlamSimulator::trueRobotPose3D (VmMatrix4x4& pose)
{
  double s, c;  
  sinCos (s, c);  
  pose[0][0] = s;  pose[0][1] = 0;  pose[0][2] = c; pose[0][3] = scale*robotX;
  pose[1][0] =-c;  pose[1][1] = 0;  pose[1][2] = s; pose[1][3] = -scale*robotY;
  pose[2][0] = 0;  pose[2][1] =-1;  pose[2][2] = 0; pose[2][3] = (robotZ+0.5)*storyHeight;
  pose[3][0] = 0;  pose[3][1] = 0;  pose[3][2] = 0; pose[3][3] = 1;    
}


void SlamSimulator::landmarkPosition3D (VmVector3& p, int x, int y, int storyNr)
{
  int lmIndex = (*story[storyNr]) (x, y);
  assert (lmIndex>=Story::LANDMARK);

  p[0] = scale*x;
  p[1] = -scale*y;
  p[2] = (VmRandom::deterministicRandom (0, x + y * 317 + 17*storyNr) + storyNr)*storyHeight;  
}


int SlamSimulator::storyOfLandmark (int landmarkId) const
{
  int storyNr;  
  for (storyNr=0; storyNr<(int) baseLandmarkId.size() && landmarkId>=baseLandmarkId[storyNr]; storyNr++);  
  return storyNr-1;  
}

