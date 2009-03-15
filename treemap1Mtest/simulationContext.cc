#include "simulationContext.h"
#include <treemapgui/tmgSlam2DView.h>

SimulationContext::SimulationContext ()
  :storeHistory (false)
{}


void SimulationContext::initTmgSimulator (const char* file)
{
  sim.load (file);
  int n = sim.nrOfLandmarks();  
  VmVector3 initialPose;
  sim.trueRobotPose (initialPose);
  VmMatrix3x3 initialPoseCov = 
    {{0.0001, 0, 0},
     {0, 0.0001, 0},
     {0, 0, 0.0001}};
  treemap2D.clear();
  treemap2D.push_back (TmSlamDriver2DL()); // start with empty treemap
  // reserve n/200 robot poses, 
  // 4 optimization step per SLAM step, maximally 3 unsuccessful optimization steps,
  // keep 3 nonlinear leaves, sparsification every 4m
  treemap().create (n, max(1000,n/10), initialPose, initialPoseCov, 4, 3, 3, 4);  
//  treemap().create (n, max(1000,n/10), initialPose, initialPoseCov, 4, 3, 40, 4);  
}


void SimulationContext::saveEstimate (const char* filename) const
{
  const TmSlamDriver2DL& tm = treemap();  
  XycVector<Landmark> landmarks;
  landmarks.reserve (tm.landmark.size()+oldPose.size());  
  for (int i=0; i<tm.landmark.size(); i++) {    
    const TmSlamDriver2DL::Landmark& lm = tm.landmark[i];
    if (lm.featureId<0) continue;
    VmVector3 p;    
    p[0] = tm.feature[lm.featureId  ].est;
    p[1] = tm.feature[lm.featureId+1].est;
    p[2] = lm.level * sim.storyHeight;
    landmarks.push_back (Landmark (p, 0, lm.level));
  }
  for (int i=0; i<oldPose.size(); i++) landmarks.push_back (oldPose[i]);
  TmgBasicSimulationContext::saveEstimate (filename, landmarks);  
}


void SimulationContext::oneStepWithSimulator (VmVector3& odo, VmMatrix3x3& odoCov, XycVector<TmSlamDriver2DL::Observation>& obs)
{
  sim.step (odo, odoCov);
  vector<TmgSimulator::Observation> obs2; 
  sim.observe (obs2);
  obs.clear();  
  for (int i=0; i<(int) obs2.size(); i++) 
    obs.push_back (TmSlamDriver2DL::Observation (obs2[i].id, obs2[i].pos, obs2[i].posCov));
}


void SimulationContext::simStep (StatisticEntry& stat)
{
  VmVector3 odo;
  VmMatrix3x3 odoCov;    
  sim.step (odo, odoCov);  
  vector<TmgSimulator::Observation> obs; 
  sim.observe (obs);
  XycVector<TmSlamDriver2DL::Observation> obs2;
  double len, maxLen = 1;  
  for (int i=0; i<(int) obs.size(); i++) { 
    len = sqrt(obs[i].pos[0]*obs[i].pos[0] + obs[i].pos[1]*obs[i].pos[1]);
    if (len>maxLen) maxLen = len;    
    obs2.push_back (TmSlamDriver2DL::Observation (obs[i].id, obs[i].pos, obs[i].posCov));
  }  

  
  double t4, t3, t2, t1, t0 = treemap().time();      
  treemap().step (odo, odoCov);
  double t0a = treemap().time();  
  treemap().setLevel (sim.robotZ);  
  double t0b = treemap().time();  
  double chi2Before = treemap().chi2 (obs2);  
  double t0c = treemap().time();  
  treemap().observe (obs2);
  t1 = treemap().time();      
  treemap().optimizeFullRuns ();
  t2 = treemap().time();
  treemap().updateGaussians ();      
  t3 = treemap().time();
/*  if (t2-t0>t3-t2 && t3-t0>0.02) {
    printf("%9d Strange timings: %6.2fms %6.2fms %6.2fms\n", treemap().slamStatistics().p, 1000*(t1-t0), 1000*(t2-t1), 1000*(t3-t2));
    printf("%6.2fms %6.2fms %6.2fms %6.2fms\n", 1000*(t0a-t0), 1000*(t0b-t0a), 1000*(t0c-t0b), 1000*(t1-t0c));    
    }  */
  if (sim.hasHitWaypoint) treemap().updateAllEstimates();
  else treemap().onlyUpdateLevel (sim.robotZ);
  treemap().computeLinearEstimate ();
  if (sim.hasHitWaypoint) { // To focus the cache on the next level
    treemap().onlyUpdateLevel (sim.robotZ);
    treemap().computeLinearEstimate ();        
  }  
  t4 = treemap().time();

  stat.fromTreemap (treemap());
  stat.timeBookkeeping = t2-t0;
  if (sim.hasHitWaypoint) {
    stat.timeFullEstimation = t4-t2;
    stat.timeEstimation     = 0;
  }  
  else {
    stat.timeEstimation  = t4-t2;
    stat.timeFullEstimation = 0;
  }  
  stat.timeTotal = stat.timeBookkeeping + stat.timeEstimation;  
  VmVector3 p;  
  double theta;  
  treemap().robotEstimate (p[0], p[1], theta);  
  p[2] = sim.robotZ*sim.storyHeight;  
  oldPose.push_back (Landmark(p, 1, sim.robotZ));
  double chi2 = treemap().chi2 (obs2);
  if (chi2>6*obs2.size()) {    
    printf ("Warning: Chi^2 error of %f (%f before) with %d observations (p=%d)\n", 
            chi2, chi2Before, obs2.size(), treemap().slamStatistics().p);
  }
}



void SimulationContext::savePng (int level, char* filename)
{
  double loX, hiX, loY, hiY;  
  sim.boundingBox (level, loX, hiX, loY, hiY);  
  double w = hiX-loX;
  double h = hiY-loY;
  loX -= 0.1*w;
  hiX += 0.1*w;
  loY -= 0.1*h;
  hiY += 0.1*h;  
  TmgSlam2DView::savePng (treemap(), level, loX, hiX, loY, hiY, 10, filename);
}



void SimulationContext::runBatchExperiment (char* filename)
{
  char fn[1000];  
  if (!TmTreemap::isCompiledWithOptimization()) printf ("WARNING: Not compiled with optimization. \n");  
  replaceSuffix (fn, filename, ".dat");  
  FILE* logFile = fopen (fn, "w");  
  printStatComment (logFile);  
  int ctr = 0;  
  try {
    initTmgSimulator (filename);
    StatisticEntry minStat, maxStat;
    while (!sim.isFinished()) {
      StatisticEntry stat;      
      simStep (stat);
      if (sim.hasHitWaypoint || stat.p%200==0) {
        printStat (logFile, minStat, maxStat);
        printStat (stdout, minStat, maxStat);        
        fflush (stdout);        
        minStat = maxStat = stat;
      }      
      else {
        minStat.min (stat);
        maxStat.max (stat);
      }      
      
      if (sim.hasHitWaypoint) { // Save a snapshot
//        treemap().printFeatureFragmentation();        
        printf ("\n\n");        
        printf ("Saving snapshot %d (Mem: %7.3fMB) of story %d\n", ctr, treemap().memory()/1.0E6, sim.robotOldZ);
        char fname[1000], txt[1000];
        sprintf (fname, ".%05d.map.png", ctr);
        replaceSuffix (txt, filename, fname);
        savePng (sim.robotOldZ, txt);        
        printf(".\n");        
        sprintf (fname, ".%05d.map.dat", ctr);
        replaceSuffix (txt, filename, fname);
        saveEstimate (txt);
        oldPose.clear();        
        printf(".\n");        
        printStatComment (stdout);        
        ctr++;
        fflush (stdout);
        sleep(1);        
      }
    }
    fclose (logFile);    
  }
  catch (runtime_error err) {
    fprintf (stderr, "%s", err.what());
    return;    
  }
}

void SimulationContext::clear()
{
  storeHistory = false;
  treemap2D.clear();
  TmgBasicSimulationContext::clear();  
}

