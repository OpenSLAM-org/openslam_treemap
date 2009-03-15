#include "tmgSlam2DView.h"

TmgSlam2DView::TmgSlam2DView (QWidget* parent, const char* name)
  :TmgMapView (parent, name)
{
}


void TmgSlam2DView::convertTmSlamDriverToDots (const TmSlamDriver2DL& treemap, int level, DotList& dot)
{
  dot.clear();  
  for (int i=0; i<treemap.landmark.size(); i++) {
    const TmSlamDriver2DL::Landmark& lm = treemap.landmark[i];      
    if (lm.featureId<0 || lm.level!=level) continue;      
    double x = treemap.feature[lm.featureId  ].est;
    double y = treemap.feature[lm.featureId+1].est;
    if (!finite(x) || !finite(y)) continue;
    char name[10];
    int ctr;      
    treemap.nameOfFeature (name, lm.featureId, ctr);
    dot.push_back (Dot (lm.featureId, x, y, name, Qt::black));    
  }
  for (int i=0; i<(int) treemap.pose.size(); i++) {
    const TmSlamDriver2DL::Pose& po = treemap.pose[i];      
    if (po.featureId<0 || po.level!=level) continue;      
    double x     = treemap.feature[po.featureId  ].est;
    double y     = treemap.feature[po.featureId+1].est;
    double theta = treemap.feature[po.featureId+2].est;
    if (finite(x) && finite(y) && finite(theta)) {
      char name[10];
      int ctr;      
      treemap.nameOfFeature (name, po.featureId, ctr);
      dot.push_back (Dot (po.featureId, x, y, name, Qt::red));    
    }    
  }  
}

  
void TmgSlam2DView::setDotsFromTreemap (const TmSlamDriver2DL& treemap, int level)
{
  DotList dot;
  convertTmSlamDriverToDots (treemap, level, dot);
  setDots (dot);
}


void TmgSlam2DView::savePng (const TmSlamDriver2DL& treemap, int level, 
			      double loX, double hiX, double loY, double hiY,
			      double scale, const char* filename)
{
  DotList dot;
  convertTmSlamDriverToDots (treemap, level, dot);
  TmgMapView::savePng (dot, loX, hiX, loY, hiY, scale, filename);  
}


void TmgSlam2DView::convertTmSlamDriverToDots (const TmSlamDriver2DP& treemap, DotList& dot)
{
  dot.clear();  
  for (int i=0; i<(int) treemap.pose2Feature.size(); i++) {
    int fid  = treemap.pose2Feature[i];      
    if (fid<0) continue;      
    double x     = treemap.feature[fid  ].est;
    double y     = treemap.feature[fid+1].est;
    double theta = treemap.feature[fid+2].est;
    if (finite(x) && finite(y) && finite(theta)) {
      char name[10];
      int ctr;
      treemap.nameOfFeature (name, fid, ctr);
      dot.push_back (Dot (fid, x, y, name, Qt::red));    
    }    
  }  
}

  
void TmgSlam2DView::setDotsFromTreemap (const TmSlamDriver2DP& treemap)
{
  DotList dot;
  convertTmSlamDriverToDots (treemap, dot);
  setDots (dot);
}

void TmgSlam2DView::savePng (const TmSlamDriver2DP& treemap, 
			      double loX, double hiX, double loY, double hiY,
			      double scale, const char* filename)
{
  DotList dot;
  convertTmSlamDriverToDots (treemap, dot);
  TmgMapView::savePng (dot, loX, hiX, loY, hiY, scale, filename);  
}


void TmgSlam2DView::highlightFeature (const TmExtendedFeatureList& fl)
{
  XycVector<int> fli;
  fli.resize (fl.size());  
  for (int i=0; i<fl.size(); i++) fli[i] = fl[i].id;
  TmgMapView::highlightFeature (fli);
}



  
