#ifndef SIMULATION_CONTEXT_H
#define SIMULATION_CONTEXT_H

#include "vectormath/vectormath.h"
#include <treemap/tmSlamDriver2DP.h>
#include <xycontainer/xycVector.h>


class SimulationContext 
{
 public:
  //! Default constructor
  SimulationContext ();  

  //! The treemap on which the algorithm works
  /*! Computation is performed on the treemap2D.back() element.
    If \c storeHistory is true, we store older results to allow to view them. */
  XycVector<TmSlamDriver2DP> treemap2D;

  //! Whether to store old results in \c treemap2D
  bool storeHistory;  

  TmSlamDriver2DP& treemap() {return treemap2D.back();}
  const TmSlamDriver2DP& treemap() const {return treemap2D.back();}

  //! The log file as a list of links
  /*! The links are sorted according to the maximum of the poses involved. */
  XycVector<TmSlamDriver2DP::Link> log;  

  //! The next link to be replayed from the log is \c log[replayCtr]
  int replayCtr;  

  class AlgorithmAbortedException : public runtime_error
    {
    public:
      AlgorithmAbortedException ()
        :runtime_error ("Algorithm aborted")
        {}
    };  

  bool isRunning () const {return !treemap2D.empty();} 
  
  
  //! Runs the experiment configured in \c filename batch mode
  void runBatchExperiment (char* filename);  
  
  //! Saves the treemap estimate for all poses into filename
  /*! File format: 
   */
  void saveEstimate (const char* filename) const;  

  bool isFinished () const {return replayCtr>=log.size();}
  

  //! Loads a log file and sort the links according to the larger pose id
  void loadLogFile (const char* filename);  

  //! Performs one movement & observation, updates treemap
  void logStep ();  

  //! Performs one movement and returns all the new links
  void oneStepOnLog (XycVector<TmSlamDriver2DP::Link>& link);  

  //! Saves an .png image of the current estimate of \c level
  void savePng (char* filename);  

  //! Stops the simulation
  void clear();  
};


#endif
