#include "mainWindow.h"
#include <qpopupmenu.h>
#include <qmenubar.h>
#include <qapplication.h>
#include <qinputdialog.h>
#include <qstatusbar.h>
#include <qfiledialog.h>
#include <qstring.h>

#include "treemap/tmSlamDriver2DP.h"
#include <vectormath/vmRandom.h>


MainWindow::MainWindow ()
  :QMainWindow (), sim(), mainCaption ("DFKI Treemap 2D with Landmarks"),
   treemapView(NULL), treemapShown(0),
   isStepPressed (false), autoRun (false), abortAlgorithm(false)
{
  setCaption (mainCaption);

  QAction*  action;  
  QPopupMenu* filePopup = new QPopupMenu (this, "MainWindow::FilePopup");

  action = new QAction ("Opens a building (.bui) file and start the SLAM algorithm on that building.", "Load building", 0, this);
  connect (action, SIGNAL(activated()), SLOT(fileOpenBui()));
  action->addTo (filePopup);  

  action = new QAction ("Save map as .png", "Save map", 0, this);
  connect (action, SIGNAL(activated()), SLOT(fileSavePngMap()));
  action->addTo (filePopup);  
  
  action = new QAction ("Save tree as .png", "Save tree", 0, this);
  connect (action, SIGNAL(activated()), SLOT(fileSavePngTree()));
  action->addTo (filePopup);  

  action = new QAction ("Quit", "&Quit", Qt::CTRL+Qt::Key_Q, this);
  connect (action, SIGNAL(activated()), SLOT(fileQuit()));
  action->addTo (filePopup);  
  menuBar()->insertItem ("&File", filePopup);  



  QPopupMenu* algorithmPopup = new QPopupMenu (this, "MainWindow::AlgorithmPopup");
  action = new QAction ("Run a single update step of the SLAM algorithm", "&Step", Qt::CTRL + Qt::Key_Return, this);
  connect (action, SIGNAL(activated()), SLOT(goStep()));
  action->addTo (algorithmPopup);

  action = new QAction ("Run the SLAM algorithm continuously without user intervention (ESC to break)", "&Autorun", 0, this);
  action->setToggleAction (true);  
  connect (action, SIGNAL(toggled(bool)), SLOT(autoRunToggled(bool)));
  action->addTo (algorithmPopup);

  action = new QAction ("Stops autorun", "Autorun off", Qt::Key_Escape, this);
  connect (action, SIGNAL(activated()), SLOT(autoRunOff()));
  action->addTo (algorithmPopup);

  action = new QAction ("Store old treemaps.", "&Store History", 0, this);
  action->setToggleAction (true);
  connect (action, SIGNAL(toggled(bool)), SLOT(storeHistoryToggled(bool)));
  action->addTo (algorithmPopup); 

  action = new QAction ("Undo the last step (only possible if history is activated)", "&Back", Qt::CTRL + Qt::Key_BackSpace, this);
  connect (action, SIGNAL(activated()), SLOT(backStep()));
  action->addTo (algorithmPopup);  

  action = new QAction ("Abort the running algorithm", "&Abort", 0, this);
  connect (action, SIGNAL(activated()), SLOT(abortTheAlgorithm()));  
  action->addTo (algorithmPopup);  

  algorithmPopup->insertSeparator();


  action = new QAction ("Calibrate the QR performance", "&Calibrate", 0, this);
  connect (action, SIGNAL(activated()), SLOT(testQRSpeed()));
  action->addTo (algorithmPopup);  

  action = new QAction ("Test random number generator", "&Test Randomnumbers", 0, this);
  connect (action, SIGNAL(activated()), SLOT(testVmRandom()));
  action->addTo (algorithmPopup);  


  action = new QAction ("Find the optimal KL-step of a node", "optimal KL-step", 0, this);
  connect (action, SIGNAL(activated()), SLOT(optimalKLStep()));
  action->addTo (algorithmPopup);

  action = new QAction ("Moves subtree inside the tree", "&Move subtree", 0, this);
  connect (action, SIGNAL(activated()), SLOT(moveSubtreeDialog()));
  action->addTo (algorithmPopup);

  action = new QAction ("Optimization step", "&Optimize", Qt::CTRL + Qt::Key_O, this);
  connect (action, SIGNAL(activated()), SLOT(optimizationStep()));  
  action->addTo (algorithmPopup);

  action = new QAction ("Recompute features / Gaussians from leaves", "&Recompute", 0, this);
  connect (action, SIGNAL(activated()), SLOT(recompute()));  
  action->addTo (algorithmPopup);  

  menuBar()->insertItem ("&Algorithm", algorithmPopup);  



  QPopupMenu* viewPopup = new QPopupMenu (this, "MainWindow::ViewPopup");
  action = new QAction ("Show the features involved IN a node", "&Show node features", Qt::CTRL+Qt::Key_F, this);
  connect (action, SIGNAL(activated()), SLOT(showNodeFeatures()));
  action->addTo (viewPopup);  

  action = new QAction ("Show the features involved BELOW a node", "Show recursive features", Qt::CTRL+Qt::Key_R, this);
  connect (action, SIGNAL(activated()), SLOT(showRecursiveFeatures()));
  action->addTo (viewPopup);  

  action = new QAction ("Show several statistics", "statistics", 0, this);
  connect (action, SIGNAL(activated()), SLOT(printStatistics()));
  action->addTo (viewPopup);

  action = new QAction ("Change scale", "scale", 0, this);
  connect (action, SIGNAL(activated()), SLOT(changeScale()));
  action->addTo (viewPopup);

  action = new QAction ("Draw nr of levels", "Nr of levels", 0, this);
  connect (action, SIGNAL(activated()), SLOT(changeLevel()));
  action->addTo (viewPopup);

  action = new QAction ("View names of poses/landmarks", "View names", 0, this);
  action->setToggleAction (true);  
  action->addTo (viewPopup);  
  

  menuBar()->insertItem ("&View", viewPopup);

  vSplitter = new QSplitter (Qt::Vertical, this);
  setCentralWidget (vSplitter);
  
  treemapView = new TmgTreemapView (vSplitter, "MainWindow::treemapView");
  treemapView->setDrawNrOfLevels (3);  
  mapView = new TmgSlam2DView (vSplitter, "MainWindow::mapView");
  mapView->connect (action, SIGNAL(toggled(bool)), SLOT(setDoShowNames(bool)));

  mapView->setScale (4);  
  mapView->setDoShowNames (false);  

  QValueList<int> sizes;
  sizes.append(2*height()/3);
  sizes.append(1*height()/3);
  vSplitter->setSizes (sizes);

  statusBar();  
}


void MainWindow::fileOpenBui ()
{
  if(sim.isRunning()) {
    statusBar()->message ("Please abort the running algorithm first");    
    return;  
  }  
  QString s( QFileDialog::getOpenFileName( QString::null, "Building *.bui", this ) );
  if ( s.isEmpty() ) return;
  runSLAM (s);  
}


void MainWindow::fileSavePngMap ()
{
    QString s( QFileDialog::getSaveFileName( QString::null, "Image *.png", this ) );
    if ( s.isEmpty() ) return;
    ((TmgMapView*)mapView)->savePng (s);    
}


void MainWindow::fileSavePngTree ()
{
    QString s( QFileDialog::getSaveFileName( QString::null, "Image *.png", this ) );
    if ( s.isEmpty() ) return;
    treemapView->savePng (s);    
}


void MainWindow::fileQuit()
{
  //int really = QMessageBox::information (this, "Quit", "Do you really want to quit?", "&Yes", "&No");
//    if (really==0) 
  abortAlgorithm = true;  
  qApp->quit();
}


void MainWindow::checkOptimalStep (int lca, int subtree, int above)
{
  TmNode* n = treemap().getNode(lca);
  if (n==NULL) throw runtime_error ("Could not find node.");  
  TmTreemap::Move move;  
  treemap().optimalKLStep (n, INT_MAX, move);
  char txt[1000];
  if (!move.isEmpty()) sprintf (txt, "%d: %d-->%d cost: %9.4fms-->%9.4fms", 
                                n->index, move.subtree->index, move.above->index, 
                                n->worstCaseUpdateCost*1000, move.cost*1000);
  else sprintf (txt, "%d: No move found.", n->index);  
  updateViews ();  
  statusBar()->message (txt);
  assert (subtree==-1 && move.cost>=n->worstCaseUpdateCost || (move.subtree->index==subtree && move.above->index==above));  
  step();  
}


void MainWindow::optimalKLStep ()
{
  int nodeIndex = QInputDialog::getInteger ("Optimal KL step for node", "Node index", -1);
  TmNode* n = treemap().getNode(nodeIndex);
  treemap().computeNonlinearEstimate();  
  if (n==NULL) return;
  TmTreemap::Move move;  
  treemap().optimalKLStep (n, INT_MAX, move);
  char txt[1000];
  if (!move.isEmpty()) sprintf (txt, "%d-->%d cost: %9.4fms-->%9.4fms", 
                                move.subtree->index, move.above->index, n->worstCaseUpdateCost*1000, move.cost*1000);
  else sprintf (txt, "No move found.");  
  updateViews ();  
  statusBar()->message (txt);
}


void MainWindow::optimizationStep ()
{
  char txt[1000];  
  if (!treemap().optimizer.optimizationQueue.empty()) sprintf (txt, "Optimizing %d ....", treemap().optimizer.optimizationQueue.front());
  else sprintf (txt, "Optimizing...");  
  statusBar()->message (txt);  
  qApp->processEvents();
  treemap().optimizer.oneKLRun ();  
  updateViews ();
  statusBar()->message (txt);
}


void MainWindow::recompute ()
{
  treemap().fullRecompute ();
  updateViews ();  
}


void MainWindow::moveSubtree (int subtree, int above)
{
  TmNode* subtreeN = treemap().getNode(subtree);
  if (subtreeN==NULL) throw runtime_error ("Could not find node ");  
  TmNode* aboveN = treemap().getNode(above);
  if (aboveN==NULL) throw runtime_error ("Could not find node ");  
  subtreeN->moveTo (aboveN);
  treemap().computeNonlinearEstimate();
  updateViews();
  statusBar()->clear();
  step();  
}


void MainWindow::moveSubtreeDialog ()
{
  int nodeIndex = QInputDialog::getInteger ("Move subtree", "Node index from", -1);
  TmNode* subtree = treemap().getNode(nodeIndex);
  if (subtree==NULL) return;
  nodeIndex = QInputDialog::getInteger ("Move subtree", "Node index to above", -1);
  TmNode* above = treemap().getNode(nodeIndex);
  if (above==NULL) return;
  subtree->moveTo (above);
  treemap().computeNonlinearEstimate();
  updateViews();
  statusBar()->clear();  
}


void MainWindow::printStatistics()
{
  char txt[1000];
  TmTreemap::SlamStatistic stat = treemap().slamStatistics ();  
  TmTreemap::TreemapStatistics stat2;
  treemap().computeStatistics (stat2, true);  
  sprintf (txt, "n=%d landmarks, m=%d measurements, p=%d poses,  %d poses marginalized, %d poses sparsified, %d nodes, %d nodes to be optimized, memory:%5.1fMB", 
           stat.n, stat.m, stat.p, stat.pMarginalized, stat.pSparsified, 
           stat2.nrOfNodes, stat2.nrOfNodesToBeOptimized,
           stat2.memory/1.0E6);
  statusBar()->message (txt);  
}


void MainWindow::step ()
{
  updateViews ();  
  qApp->processEvents();  
  while (!isStepPressed && !abortAlgorithm) qApp->processOneEvent();
  isStepPressed = autoRun;
  if (abortAlgorithm) throw SimulationContext::AlgorithmAbortedException();
  nextStep ();  
}


void MainWindow::showNodeFeatures ()
{
  TmExtendedFeatureList fl;  
  int nodeIndex = QInputDialog::getInteger ("Show features involved IN a node", "Node index", -1);
  TmNode* n = treemap().getNode(nodeIndex);
  if (n!=NULL) n->computeFeaturesInvolved (fl);
  mapView->highlightFeature (fl);  
}


void MainWindow::showRecursiveFeatures ()
{
  TmExtendedFeatureList fl;  
  int nodeIndex = QInputDialog::getInteger ("Show features involved BELOW a node", "Node index", -1);
  TmNode* n = treemap().getNode(nodeIndex);
  if (n!=NULL) treemap().computeFeaturesInvolvedBelow (n, fl);
  mapView->highlightFeature (fl);  
}


void MainWindow::goStep ()
{
  if (treemapShown<(int) sim.treemap2D.size()-1) {
    treemapShown++;
    updateViews();
  }
  else isStepPressed = true;
}


void MainWindow::storeHistoryToggled (bool on)
{
  sim.storeHistory = on;  
}


void MainWindow::autoRunToggled (bool on)
{
  autoRun = on;
  isStepPressed = on;  
}


void MainWindow::changeScale ()
{
  double scale = QInputDialog::getDouble ("Change map scale (pixel/m)", "scale", mapView->getScale());
  if (scale>=1) mapView->setScale (scale);
}


void MainWindow::changeLevel ()
{
  int levels = QInputDialog::getInteger ("Show how many levels of tree", "levels", treemapView->getDrawNrOfLevels());
  if (levels>=3) treemapView->setDrawNrOfLevels (levels);
}




void MainWindow::autoRunOff ()
{
  autoRun = isStepPressed = false;  
}


void MainWindow::abortTheAlgorithm ()
{
  abortAlgorithm = true;  
}


void MainWindow::backStep ()
{
  if (treemapShown>0) treemapShown--;
  updateViews();
}


void MainWindow::printEst ()
{
  printf ("Estimate: {");  
  for (int i=0; i<(int) treemap().feature.size(); i++) {
    char txt[5];
    featureToString (txt, i);
    if (treemap().feature[i].isDefined())
      printf ("%d%s:%2.2f, ", i, txt, treemap().feature[i].est);
  }  
  printf ("}\n");  
}

void MainWindow::assertEst ()
{
  for (int i=0; i<(int) treemap().feature.size(); i++) {
    double est = treemap().feature[i].est;    
    assert (fabs(est-(i+10))<1E-9);
  }  
}


void MainWindow::updateViews ()
{
  if (0<=treemapShown && treemapShown<(int) sim.treemap2D.size()) {
    treemapView->setTreemap (&sim.treemap2D[treemapShown]);    
    mapView->setDotsFromTreemap (sim.treemap2D[treemapShown]);
    mapView->enlargeArea ();    
  }
  else {
    treemapView->setTreemap (NULL);
    mapView->clearDots();
  }  
}


void MainWindow::ensureRobotVisible ()
{
  if (sim.isRunning()) {    
    double x, y, theta;
    treemap().robotEstimate (x, y, theta);
    mapView->ensureVisible (x, y, 0);  
  }  
}



void MainWindow::updateViewsAndMap ()
{
  updateViews ();
  ensureRobotVisible ();  
}


void MainWindow::nextStep ()
{
  if (sim.storeHistory) {    
    sim.treemap2D.push_back (sim.treemap2D.back());
    treemapShown = (int) sim.treemap2D.size()-1;
  }  
}

  

void MainWindow::fullUpdate (bool doPrintEst)
{
  treemap().assertIt();  
  treemap().computeNonlinearEstimate();
  treemap().assertIt();  
  if (doPrintEst) printEst();
  step();
}
  

void MainWindow::testVmRandom ()
{
  double sum=0, sum2=0;  
  int n = 100000;  
  VmRandom rnd;  
  for (int i=0; i<n; i++) {
    double x = rnd.gauss();
    sum += x;
    sum2 += x*x;
  }
  sum   /= n;
  sum2  /= n;  
  char txt[1000];  
  sprintf (txt, "%d Gaussians: mean %f stddev: %f", n, sum, sqrt(sum2));  
  statusBar()->message (txt);  
}


void MainWindow::testQRSpeed()
{
  double coef[4];
  if (!TmTreemap::isCompiledWithOptimization()) printf("Warning: Not compiled with optimization.");  
  TmTreemap::calibrateGaussianPerformance (300, coef, 0, "qrspeed.dat");
  printf(" Gaussian: time in microseconds %f + %f*n + %f*n^2 + %f*n^3\n", coef[0]*1E6, coef[1]*1E6, coef[2]*1E6, coef[3]*1E6);
  printf(" Normalized: %f + %f*n + %f*n^2 + n^3\n", coef[0]/coef[3], coef[1]/coef[3], coef[2]/coef[3]);  
}


void MainWindow::setArea ()
{
/*
  double loX, hiX, loY, hiY;
  sim.sim.area (loX, hiX, loY, hiY);
  double w = hiX-loX;
  double h = hiY-loY;
  loX -= 0.2*w;
  hiX += 0.2*w;
  loY -= 0.2*h;
  hiY += 0.2*h;
  mapView->setArea (loX, hiX, loY, hiY);  
*/
}




void MainWindow::runSLAM (const char* filename)
{
  try {
    abortAlgorithm = false;    
    setCaption (mainCaption+" ("+filename+")");    
    sim.loadLogFile (filename);    
    setArea ();
    int ctr = 0;     
    while (!sim.isFinished()) {
      sim.logStep ();  
      updateViewsAndMap ();
      char txt[1000];
/*      snprintf (txt, 1000, "n=%5d p=%5d (%5d mar, %5d spar) m=%6d wc=%6.3fms | %s",
               stat.n, stat.p, stat.pMarginalized, stat.pSparsified, stat.m, stat.worstCaseUpdateCost*1000,
               stat.report.c_str());      
      statusBar()->message (txt);      */
      ctr++;
      step();
    }
  }
  catch (SimulationContext::AlgorithmAbortedException err) {
    abortAlgorithm = false;    
    statusBar()->message ("algorithm aborted");    
  }  
  catch (runtime_error err) {
    fprintf (stderr, "Error in runSLAM: %s\n", err.what());
    return;    
  }
  updateViewsAndMap ();
}



