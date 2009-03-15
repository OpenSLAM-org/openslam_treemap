#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "simulationContext.h"
#include <treemapgui/tmgTreemapView.h>
#include <treemapgui/tmgMapView.h>
#include <treemapgui/tmgSlam2DView.h>
#include <treemapgui/tmgSimulator.h>
#include <vectormath/vectormath.h>

#include <qmainwindow.h>
#include <qaction.h>
#include <qsplitter.h>
#include <stdio.h>

class MainWindow : public QMainWindow 
{
  Q_OBJECT
    public:

  MainWindow ();  

  SimulationContext sim;  

  TmSlamDriver2DL& treemap() {return sim.treemap();}
  

  //! Widget deviding tree view from map view
  QSplitter* vSplitter;  

  //! The widget showing \c treemap
  TmgTreemapView* treemapView;

  //! The widet showing \c treemap2D as a map
  TmgSlam2DView* mapView;
  
  //! Caption of the main window (without the filename shown)
  QString mainCaption;  

  //! Which index of \c treemap2D is currently shown
  int treemapShown;  

  bool isStepPressed;  

  bool autoRun;  

  bool abortAlgorithm;  


  void nextStep ();  


  public slots:

  //! Start assorted algorithm test codes
  void testQRSpeed ();  
  void testVmRandom ();  
  void testSimulator ();  
  void runSLAM (const char* filename);
  void recompute();  
  void showNodeFeatures ();
  void showRecursiveFeatures ();  
  void updateViews(int level=-1);  
  void updateViewsAndMap (int level=-1);
  void ensureRobotVisible ();  

  void printStatistics();

  void optimalKLStep();  

  void moveSubtreeDialog();  
  void optimizationStep ();  
  

  //! Called when ctrl Enter is pressed
  void goStep ();  

  void autoRunToggled (bool on);  
  void autoRunOff ();  

  void storeHistoryToggled (bool on);  

  void changeScale ();  

  void changeLevel ();  

  //! Called when the algorithm should be aborted
  void abortTheAlgorithm ();  

  //! Called when ctrl-Backspace is pressed
  void backStep ();  
  
  //! Refresh the GUI and wait for the user pressing ctrl-A
  void step ();  

  //! Prints the whole treemap estimate
  void printEst ();  

  void assertEst ();  

  void fileOpenBui ();  
  void fileSavePngMap ();  
  void fileSavePngTree ();  
  void fileQuit(); 

  void fullUpdate (bool doPrintEst=true);  

 protected:

  //! Sets \c mapView->setArea according to the area of the building in \c sim
  void setArea ();

  void moveSubtree (int subtree, int above);

  //! Asserts that the optimal step for node lca is \c subtree->\c above
  void checkOptimalStep (int lca, int subtree, int above);
};


#endif
