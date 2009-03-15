#include "tmgBasicSimulationContext.h"

TmgBasicSimulationContext::TmgBasicSimulationContext ()
{
}

void TmgBasicSimulationContext::clear()
{
  sim.clear();
  oldPose.clear();  
}


void TmgBasicSimulationContext::saveEstimateBinary (const char* filename, const XycVector<Landmark>& landmarks) const
{
  FILE* f = fopen (filename, "w");
  if (f==NULL) throw runtime_error ("Could not open "+string(filename)+" for writing");

  // Find bounding box
  VmVector3 min, max;  
  bool first=true;  
  for (int i=0; i<landmarks.size(); i++) {
    VmVector3 p;
    vmCopy (landmarks[i].v, p);    
    if (first || p[0]<min[0]) min[0] = p[0];
    if (first || p[0]>max[0]) max[0] = p[0];
    if (first || p[1]<min[1]) min[1] = p[1];
    if (first || p[1]>max[1]) max[1] = p[1];
    if (first || p[2]<min[2]) min[2] = p[2];
    if (first || p[2]>max[2]) max[2] = p[2];
    first = false;    
  } 
  VmVector3 dim;
  dim[0] = max[0]-min[0];
  dim[1] = max[1]-min[1];
  dim[2] = max[2]-min[2];
  double maxDim;
  if (dim[0]>dim[1]) maxDim = dim[0];
  else maxDim = dim[1];
  if (dim[2]>maxDim) maxDim = dim[2];  
  for (int i=0; i<3; i++) {
    min[i] -= 0.1*dim[i]+0.01*maxDim;
    max[i] += 0.1*dim[i]+0.01*maxDim;
  }  

  // Store box
  for (int i=0; i<3; i++) {
    fwrite (&min[i], sizeof(double), 1, f);
    fwrite (&max[i], sizeof(double), 1, f);
  }  

  // Store landmarks
  for (int i=0; i<landmarks.size(); i++) {
    VmVector3 p;
    vmCopy (landmarks[i].v, p);
    unsigned short xx = (unsigned short) round((p[0]-min[0])*65536/(max[0]-min[0]));
    unsigned short yy = (unsigned short) round((p[1]-min[1])*65536/(max[1]-min[1]));
    unsigned short zz = (unsigned short) round((p[2]-min[2])*65536/(max[2]-min[2]));
    unsigned char flag = landmarks[i].flag;
    unsigned char story = landmarks[i].story;
    fwrite (&xx, sizeof(xx), 1, f);
    fwrite (&yy, sizeof(yy), 1, f);
    fwrite (&zz, sizeof(zz), 1, f);
    fwrite (&flag, sizeof(flag), 1, f);    
    fwrite (&story, sizeof(story), 1, f);    
  }

  fclose(f);    
}


void TmgBasicSimulationContext::saveEstimate (const char* filename, const XycVector<Landmark>& landmarks) const
{
  FILE* f = fopen (filename, "w");
  if (f==NULL) throw runtime_error ("Could not open "+string(filename)+" for writing");

  // Store landmarks
  for (int i=0; i<landmarks.size(); i++) {
    VmVector3 p;
    vmCopy (landmarks[i].v, p);
    fprintf (f, "%10.8f %10.8f %10.8f 1\n", p[0], p[1], p[2]);
  }

  fclose(f);    
}


void TmgBasicSimulationContext::StatisticEntry::fromTreemap (TmTreemap& tm)
{
  TmTreemap::SlamStatistic stat = tm.slamStatistics();  
  n              = stat.n;
  m              = stat.m;
  p              = stat.p;
  pMarginalized  = stat.pMarginalized;
  pSparsified    = stat.pSparsified;
  if (tm.root!=NULL) worstCaseUpdateCost = tm.root->worstCaseUpdateCost;  
  else worstCaseUpdateCost = 0;  
  TmTreemap::TreemapStatistics stat2;
  tm.computeStatistics (stat2, false);  
  snapshotCtr       = 0;  
  report            = tm.getAndClearReport ();  
  nrOfNodes         = stat2.nrOfNodes;
  nrOfToBeOptimized = stat2.nrOfNodesToBeOptimized;  
  mem               = 0; // to expensive  
}


void TmgBasicSimulationContext::StatisticEntry::min (const StatisticEntry& s2)
{
  if (s2.n                   < n                  ) n                   = s2.n;
  if (s2.m                   < m                  ) m                   = s2.m;
  if (s2.p                   < p                  ) p                   = s2.p;
  if (s2.pMarginalized       < pMarginalized      ) pMarginalized       = s2.pMarginalized;
  if (s2.pSparsified         < pSparsified        ) pSparsified         = s2.pSparsified;
  if (s2.worstCaseUpdateCost < worstCaseUpdateCost) worstCaseUpdateCost = s2.worstCaseUpdateCost;
  if (s2.timeBookkeeping     < timeBookkeeping    ) timeBookkeeping     = s2.timeBookkeeping;
  if (s2.timeEstimation      < timeEstimation     ) timeEstimation      = s2.timeEstimation;
  if (s2.timeFullEstimation  < timeFullEstimation ) timeFullEstimation  = s2.timeFullEstimation;     
  if (s2.timeTotal           < timeTotal          ) timeTotal           = s2.timeTotal;  
  if (s2.snapshotCtr         < snapshotCtr        ) snapshotCtr         = s2.snapshotCtr;  
  if (s2.nrOfNodes           < nrOfNodes          ) nrOfNodes           = s2.nrOfNodes;
  if (s2.nrOfToBeOptimized   < nrOfToBeOptimized  ) nrOfToBeOptimized   = s2.nrOfToBeOptimized;
  if (s2.nrOfGaussianUpdates < nrOfGaussianUpdates) nrOfGaussianUpdates = s2.nrOfGaussianUpdates;  
  if (s2.mem                 < mem                ) mem                 = s2.mem;  
}


void TmgBasicSimulationContext::StatisticEntry::max (const StatisticEntry& s2)
{
  if (s2.n                   > n                  ) n                   = s2.n;
  if (s2.m                   > m                  ) m                   = s2.m;
  if (s2.p                   > p                  ) p                   = s2.p;
  if (s2.pMarginalized       > pMarginalized      ) pMarginalized       = s2.pMarginalized;
  if (s2.pSparsified         > pSparsified        ) pSparsified         = s2.pSparsified;
  if (s2.worstCaseUpdateCost > worstCaseUpdateCost) worstCaseUpdateCost = s2.worstCaseUpdateCost;
  if (s2.timeBookkeeping     > timeBookkeeping    ) timeBookkeeping     = s2.timeBookkeeping;
  if (s2.timeEstimation      > timeEstimation     ) timeEstimation      = s2.timeEstimation;
  if (s2.timeFullEstimation  > timeFullEstimation ) timeFullEstimation  = s2.timeFullEstimation;     
  if (s2.timeTotal           > timeTotal          ) timeTotal           = s2.timeTotal;  
  if (s2.snapshotCtr         > snapshotCtr        ) snapshotCtr         = s2.snapshotCtr;  
  if (s2.nrOfNodes           > nrOfNodes          ) nrOfNodes           = s2.nrOfNodes;
  if (s2.nrOfToBeOptimized   > nrOfToBeOptimized  ) nrOfToBeOptimized   = s2.nrOfToBeOptimized;  
  if (s2.nrOfGaussianUpdates > nrOfGaussianUpdates) nrOfGaussianUpdates = s2.nrOfGaussianUpdates;  
  if (s2.mem                 > mem                ) mem                 = s2.mem;  
}


void TmgBasicSimulationContext::replaceSuffix (char* result, char* filename, char* suffix)
{
  strcpy (result, filename);
  int i;  
  for (i=strlen(result); i>=0 && result[i]!='.';i--);  
  if (i<0) i=strlen(result);  
  strcpy (result+i, suffix);  
}


void TmgBasicSimulationContext::printStatComment (FILE* f)
{
  fprintf (f, "#$1ct    $2n      $3m      $4p   $5pMag $6pSpars $7-wcCost $8+wcCost  $9-tBook  $10+tBook "\
           "$11-tEst  $12+tEst $13+tFEst     $14-t     $15+t $16nod. $17ntO$18nGupd  $19mem      \n");  
}


void TmgBasicSimulationContext::printStat (FILE* f, const StatisticEntry& minStat, const StatisticEntry& maxStat)
{
  //           $1  $2  $3  $4  $5  $6   $7    $8    $9    $10   $11   $12   $13   $14  $15   $16 $17 $18 $19
  fprintf (f, "%4d %8d %8d %8d %8d %8d %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %6d %4d %4d %12d\n",
           maxStat.snapshotCtr, /* $1 */
           maxStat.n /*$2*/, maxStat.m /*$3*/, maxStat.p /*$4*/, maxStat.pMarginalized /*$5*/, maxStat.pSparsified /*$6*/,
           minStat.worstCaseUpdateCost /*$7*/, maxStat.worstCaseUpdateCost /*$8*/,
           minStat.timeBookkeeping     /*$9*/, maxStat.timeBookkeeping     /*$10*/,
           minStat.timeEstimation     /*$11*/, maxStat.timeEstimation      /*$12*/,
           maxStat.timeFullEstimation /*$13*/, 
           minStat.timeTotal /*$14*/, maxStat.timeTotal /*$15*/,
           maxStat.nrOfNodes /*$16*/, maxStat.nrOfToBeOptimized /*$17*/, maxStat.nrOfGaussianUpdates /*$18*/,
           maxStat.mem /*$19*/
           );  
}
