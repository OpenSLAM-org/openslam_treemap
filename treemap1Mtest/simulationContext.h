#ifndef SIMULATION_CONTEXT_H
#define SIMULATION_CONTEXT_H

#include "treemapgui/tmgSimulator.h"
#include "treemapgui/tmgBasicSimulationContext.h"
#include "vectormath/vectormath.h"
#include <treemap/tmSlamDriver2DL.h>
#include <xycontainer/xycVector.h>


class SimulationContext : public TmgBasicSimulationContext
{
 public:
  //! Default constructor
  SimulationContext ();  

  //! The treemap on which the algorithm works
  /*! Computation is performed on the treemap2D.back() element.
    If \c storeHistory is true, we store older results to allow to view them. */
  XycVector<TmSlamDriver2DL> treemap2D;

  //! Whether to store old results in \c treemap2D
  bool storeHistory;  

  TmSlamDriver2DL& treemap() {return treemap2D.back();}
  const TmSlamDriver2DL& treemap() const {return treemap2D.back();}


  bool isRunning () const {return !treemap2D.empty();} 
  
  
  //! Runs the experiment configured in \c filename batch mode
  void runBatchExperiment (char* filename);  
  
  //! Saves the treemap estimate for all landmarks and \c robotTrajectory into filename
  /*! File format: 
   */
  void saveEstimate (const char* filename) const;  

  void initTmgSimulator (const char* file);

  //! Performs one movement & observation, updates treemap
  void simStep (StatisticEntry& stat);  

  //! Performs one movement & observation
  void oneStepWithSimulator (VmVector3& odo, VmMatrix3x3& odoCov, XycVector<TmSlamDriver2DL::Observation>& obs);  

  //! Saves an .png image of the current estimate of \c level
  void savePng (int level, char* filename);  

  //! Stops the simulation
  void clear();  
};


#endif
