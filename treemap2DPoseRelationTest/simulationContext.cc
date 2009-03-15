#include "simulationContext.h"
#include <treemapgui/tmgSlam2DView.h>

SimulationContext::SimulationContext ()
  :treemap2D(), storeHistory (false), log(), replayCtr(0)
{}


//! Compare larger link index predicate used for sorting with STL
/*!
  Use 
  \c sort (v.begin()+i, v.begin()+j+1, LessOnLinks());
 */
class LessOnLinks
{
 public:
  LessOnLinks () {}
  bool operator() (const TmSlamDriver2DP::Link& a, const TmSlamDriver2DP::Link& b)
    {
      return a.largerPose() < b.largerPose();
    }  
};


void SimulationContext::loadLogFile (const char* filename)
{  
  clear();
  // Load the log from \c filename

  FILE* f = fopen(filename, "r");
  if (f==NULL) throw runtime_error ("Could not open "+string(filename)+" for reading");
  while (true) {
    char buffer[1024];
    char* result;
    result = fgets (buffer, 1023, f);
    if (result==NULL) break;
    buffer[1023]='\0';
    int lastNonSpace=-1;
    bool hasComment = false;
    for (int i=0; buffer[i]!='\0'; i++) {
      // read line
      char ch = buffer[i];
      if (ch=='#') {
	buffer[lastNonSpace+1]='\0';
	hasComment = true;
	break;
      }
      else if (ch=='\n') {
	buffer[lastNonSpace+1]='\0';
	break;
      }
      else if (ch!=' ') lastNonSpace=i;
    }
    buffer[lastNonSpace+1]='\0';

    TmSlamDriver2DP::Link link;    
    int id;
    char dummy[100];    
    // parse link
    if (sscanf (buffer, "RELATION %d %d %d%[ ,] %lf %lf %lf%[ ,] %lf %lf %lf %lf %lf %lf",
                &link.poseA, &link.poseB, &id, dummy,
                &link.d[0], &link.d[1], &link.d[2], dummy,
                &link.dCov[0][0], &link.dCov[1][1], &link.dCov[2][2],
                &link.dCov[0][1], &link.dCov[0][2], &link.dCov[1][2])==14) {
      link.dCov[1][0] = link.dCov[0][1];
      link.dCov[2][0] = link.dCov[0][2];
      link.dCov[2][1] = link.dCov[1][2];      
      double factor = M_PI/180;      
      link.d [2] *= factor;
      link.dCov [0][2] *= factor;
      link.dCov [1][2] *= factor;
      link.dCov [2][2] *= factor;      
      link.dCov [2][0] *= factor;
      link.dCov [2][1] *= factor;
      link.dCov [2][2] *= factor; // that's actually factor^2
      log.push_back (link);      
    }
    else if (strlen(buffer)>0) throw runtime_error ("Could not parse: "+string(buffer));    
  }
  fclose(f);

  // Add a link of pose 0 to -1 if it is not existing yet
  bool hasGround=false;
  VmMatrix3x3 meanCov;
  vmZero (meanCov);  
  int i;  
  for (i=0; i<log.size(); i++) {
    vmAdd (meanCov, meanCov, log[i].dCov);    
    if (log[i].poseB==-1) hasGround=true;
  }  
  if (log.size()>0) vmScale (meanCov, 1.0/log.size());  
  if (!hasGround) {
    VmVector3 d;
    vmZero (d);        
    log.push_back (TmSlamDriver2DP::Link(0, -1, d, meanCov));
  }  

  // The log may be from a batch mode scenario so we have to sort it according to some
  // artificial chronology
  sort (log.begin(), log.end(), LessOnLinks());  

  // Initialize treemap
  treemap2D.clear();
  treemap2D.push_back (TmSlamDriver2DP()); // start with empty treemap
  treemap().create ();  
}

void SimulationContext::saveEstimate (const char* filename) const
{
  // not implemented
  assert (false);  
}


void SimulationContext::oneStepOnLog (XycVector<TmSlamDriver2DP::Link>& link)
{
  link.clear();
  if (isFinished()) return;  
  int largerPose = log[replayCtr].largerPose();
  while (replayCtr<log.size() && log[replayCtr].largerPose()<=largerPose) {    
    link.push_back (log[replayCtr]);
    replayCtr++;
  }  
}


void SimulationContext::logStep ()
{
  XycVector<TmSlamDriver2DP::Link> links;
  oneStepOnLog (links);
  
  for (int i=0; i<links.size(); i++) treemap().addLink (links[i]);

  treemap().optimizeFullRuns ();
  treemap().updateGaussians ();      
  treemap().computeLinearEstimate ();
}



void SimulationContext::savePng (char* filename)
{
  // not implemented
  // TmgSlam2DView::savePng (treemap(), level, loX, hiX, loY, hiY, 10, filename);
  assert (false);  
}



void SimulationContext::runBatchExperiment (char* filename)
{
/*  char fn[1000];  
  if (!TmTreemap::isCompiledWithOptimization()) printf ("WARNING: Not compiled with optimization. \n");  
  replaceSuffix (fn, filename, ".dat");  
  FILE* logFile = fopen (fn, "w");  
//  printStatComment (logFile);  
  int ctr = 0;  
  try {
    loadLogFile (filename);
    StatisticEntry minStat, maxStat;
    while (!sim.isFinished()) {
//      StatisticEntry stat;      
      simStep ();
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
*/
}

void SimulationContext::clear()
{
  storeHistory = false;
  treemap2D.clear();
  log.clear();
  replayCtr = 0;  
}

