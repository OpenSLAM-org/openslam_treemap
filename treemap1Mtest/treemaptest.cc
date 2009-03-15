#include "treemapgui/tmgTreemapView.h"
#include "mainWindow.h"
#include "simulationContext.h"
#include <qapplication.h>
#include <qmainwindow.h>

//Returns the arg which contains "-batch" or -1 if there is none
int argIdx (int argc, char** argv, char* token)
{
  for (int i=1; i<argc; i++) if (strcmp(argv[i], token)==0) return i;
  return -1;  
}


int noArgIdx (int argc, char** argv)
{
  for (int i=1; i<argc; i++) if (argv[i][0]!='-') return i;
  return -1;  
}


int interactive (int argc, char** argv) 
{
  if (argc==1) {
    fprintf (stderr, "treemap1Mtest [-batch] building.bui\n\n");
  }  
  int bai = argIdx (argc, argv, "-batch");
  if (bai>=0 && bai+1<argc) {
    SimulationContext sim;
    sim.runBatchExperiment(argv[bai+1]);
    return 0;
  }
  else {    
    QApplication app (argc, argv);
    MainWindow* gui = new MainWindow;
    gui->show();
    gui->showMaximized();
    app.setMainWidget (gui);
    int loadFileIdx = noArgIdx (argc, argv);    
    if (loadFileIdx>=0) gui->runSLAM (argv[loadFileIdx]);    
    else return app.exec();
  }    
}

  
int main (int argc, char** argv)
{
  return interactive (argc, argv);  
}

#ifdef DECLARE_MAIN__
// We apparently need this to link to libf2c.a
extern "C" {
  int MAIN__ (int argc, char** argv)
{
  return main (argc, argv);  
}
}
#endif
