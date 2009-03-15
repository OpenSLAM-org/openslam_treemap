#ifndef TMGBASICSIMULATIONCONTEXT_H
#define TMGBASICSIMULATIONCONTEXT_H
//!\author Udo Frese

/*!\file tmgBasicSimulationContext.h Contains some usefull routine for
   writing batch mode simulation programs based on \c TmgSimulator. Especially class \c TmgBasicSimulationContext
*/

#include "treemapgui/tmgSimulator.h"
#include "vectormath/vectormath.h"
#include <treemap/tmSlamDriver2DL.h>
#include <xycontainer/xycVector.h>

class TmgBasicSimulationContext
{
 public:
  //! Default constructor
  TmgBasicSimulationContext ();  

  class AlgorithmAbortedException : public runtime_error
    {
    public:
      AlgorithmAbortedException ()
        :runtime_error ("Algorithm aborted")
        {}
    };  

  //! The simulator used
  TmgSimulator sim;  

  class Landmark
  {
  public:
    VmVector3 v;
    int flag;
    int story;    
    
    Landmark () :story(0){v[0]=v[1]=v[2]=0;}
    Landmark (const VmVector3& v2) :flag(0), story(0) {v[0]=v2[0]; v[1]=v2[1]; v[2]=v2[2];}
    Landmark (const VmVector3& v2,int flag, int story) :flag(flag), story(story) {v[0]=v2[0]; v[1]=v2[1]; v[2]=v2[2];}
    VmReal& operator [] (int i) {return v[i];}
    VmReal  operator [] (int i) const {return v[i];}
  };


  //! Old poses of the robot as landmarks with flag=1
  XycVector<Landmark> oldPose;  

  class StatisticEntry : public TmTreemap::SlamStatistic
    {
    public:
      //! as reported by \c TmTreemap::root->worstCaseUpdateCost()
      double worstCaseUpdateCost;      
      //! Time spent in optimizing the tree without any linear algebra
      double timeBookkeeping;
      //! Time spent in updating Gaussians (upward pass)
      double timeEstimation;
      //! Time spent in updating Gaussians and computing the estimate (up and down)
      /*! If not a complete estimate is computed this time is left 0. */
      double timeFullEstimation;      
      //! Total computationTime
      double timeTotal; 

      //! Nr the line in the .dat file and nr. of the image / map data stored
      int snapshotCtr;

      //! Overall number of nodes in the tree
      int nrOfNodes;
      //! Number of nodes in the optimization queue
      int nrOfToBeOptimized;
      //! Nr of Gaussians updated in the upward pass
      int nrOfGaussianUpdates;      
      //! Memory consumption in bytes
      int mem;      

      string report;
      
      StatisticEntry ()
        :SlamStatistic(), worstCaseUpdateCost(0),
        timeBookkeeping(0), timeEstimation(0), timeFullEstimation(0), timeTotal(0),
        snapshotCtr (0), nrOfNodes (0), nrOfToBeOptimized (0), nrOfGaussianUpdates (0),
        mem (0)
        {}      

      void fromTreemap (TmTreemap& tm);      

      void min (const StatisticEntry& s2);
      void max (const StatisticEntry& s2);      
    };

  //! Saves \c landmarks and \c poses into \c filename in the special format we use
  void saveEstimateBinary (const char* filename, const XycVector<Landmark>& landmarks) const;  

  //! Saves \c landmarks and \c poses into \c filename in a text form that can be plotted by gnuplot
  void saveEstimate (const char* filename, const XycVector<Landmark>& landmarks) const;  

  //! Returns to empty state
  void clear();  

  //! Replaces the suffix of \c filename by \c suffix
  static void replaceSuffix (char* result, char* filename, char* suffix);    

  //! Prints one line into the performance log FILE \c f
  /*! \c minStat is the minimum time/... encountered. \c maxStat the corresponding maximum */
  static void printStat (FILE* f, const StatisticEntry& minStat, const StatisticEntry& maxStat);  

  //! Prints a comment line labeling the different columns printed by \c printStat
  static void printStatComment (FILE* f);  
};


#endif
