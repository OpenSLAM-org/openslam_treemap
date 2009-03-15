//!\author Udo Frese

/*!\file tmgMapView.cc Contains the implementation of class \c
  TmgMapView being a QT widget that shows a treemap (\c
  TmSlamDriver2DL) estimate as a 2D landmark map.
*/

#include "tmgMapView.h"
#include <qimage.h>
#include <qpainter.h>
#include <qcolor.h>
#include <qmessagebox.h>
#include <algorithm>
#include <stdexcept>
#include <math.h>
#include <algorithm>

using namespace std;

TmgMapView::TmgMapView (QWidget* parent, const char* name)
  :QScrollView (parent, name), dot(), doShowNames(true),
   scale(1), loX(0), hiX(0), loY(0), hiY(0),
   highlightF()
{
  viewport()->setEraseColor (Qt::white);  
}


void TmgMapView::ensureVisible (double x, double y, double margin)
{
  int height = (int) ceil((hiY-loY)*scale);  
  int xp = int((x-loX)*scale);
  int yp = height-int((y-loY)*scale);
  int mp = int(margin*scale)+5;
  QScrollView::ensureVisible (xp, yp, mp, mp);  
}


void TmgMapView::clear()
{
  dot.clear();  
  doShowNames = true;  
  scale = 1;
  loX = hiX = loY = hiY = 0;
  highlightF.clear();  
}


void TmgMapView::setDots (const DotList& dot)
{
  this->dot = dot;
  emit updateContent();  
}  


void TmgMapView::clearDots ()
{
  dot.clear();
  emit updateContent();  
}  
  
  
void TmgMapView::addDot (const Dot& newDot)
{
  dot.push_back (newDot);
  emit updateContent();    
}


void TmgMapView::setScale (double scale)
{
  this->scale = scale;
  emit updateContentSize();  
}


void TmgMapView::setDoShowNames (bool doShowNames)
{
  this->doShowNames = doShowNames;
  emit updateContent ();  
}


void TmgMapView::updateContent ()
{
  emit viewport()->update();  
}


void TmgMapView::highlightFeature (const XycVector<int>& fl)
{
  highlightF = fl;
  sort (highlightF.begin(), highlightF.end());  
  emit updateContent();  
}


void TmgMapView::updateContentSize ()
{
  resizeContents (int ((hiX-loX)*scale), int ((hiY-loY)*scale));
  emit updateGeometry();
  emit updateContent ();  
}


void TmgMapView::boundingBox (const DotList& dot, double& loX, double& hiX, double& loY, double& hiY, bool enlarge)
{
  bool first = !enlarge;  
  for (int i=0; i<dot.size(); i++) {
    const Dot& d = dot[i];
    if (finite(d.x) && finite(d.y)) {
      if (first || d.x<loX) loX = d.x;
      if (first || d.x>hiX) hiX = d.x;
      if (first || d.y<loY) loY = d.y;
      if (first || d.y>hiY) hiY = d.y;
      first = false;        
    }
  }      
}


void TmgMapView::setArea (double loX, double hiX, double loY, double hiY, bool enlarge)
{
  if (enlarge) boundingBox (dot, loX, hiX, loY, hiY, true);    
  
  if (loX!=this->loX || hiX!=this->hiX ||
      loY!=this->loY || hiY!=this->hiY) {  
    this->loX = loX;
    this->hiX = hiX;
    this->loY = loY;
    this->hiY = hiY;
    updateContentSize ();    
  }
}


void TmgMapView::enlargeArea ()
{
  double alpha = 0.05;  
  double bLoX, bHiX, bLoY, bHiY;
  boundingBox (dot, bLoX, bHiX, bLoY, bHiY, false);
  double length = max(bHiX-bLoX, bHiY-bLoY);
  if (loX <= bLoX-alpha*length &&
      bHiX+alpha*length <= hiX &&
      loY <= bLoY-alpha*length && 
      bLoY+alpha*length <= hiY) {
    // it fits
  }
  else {
    // enlarge it
    loX = bLoX - 2*alpha*length;
    hiX = bHiX + 2*alpha*length;
    loY = bLoY - 2*alpha*length;
    hiY = bHiY + 2*alpha*length;
    updateContentSize ();
  }  
}


void TmgMapView::drawContents ( QPainter * p, int clipx, int clipy, int clipw, int cliph )
{
  clipx = clipx; clipy=clipy; clipw=clipw; cliph=cliph; // to avoid warning  
  draw (p);  
}


bool TmgMapView::isElement (const XycVector<int>& list, int id)
{
  int lo=0, hi=(int) list.size()-1;
  if (hi<0) return false; // lmList empty
  if (id<list[lo]) return false;
  if (list[hi]<id) return false;    
  if (list[hi]==id) return true;
  // Invariant: lm is contained in [lo..hi]
  while (lo<hi) {
    int mid = (lo+hi)/2; // it is important that 1/2==0
    int lmMid = list[mid];
    if (id<lmMid) hi = mid;
    else lo = mid+1;
  }
  return (lo>0 && list[lo-1]==id);
}


void TmgMapView::draw (QPainter* p)
{
//  int width  = (int) ceil((hiX-loX)*scale);
  int height = (int) ceil((hiY-loY)*scale);  
  for (int i=0; i<dot.size(); i++) {
    Dot& d = dot[i];      
    if (!finite(d.x) || !finite(d.y)) continue;
    int xx = (int) round((d.x-loX)*scale);
    int yy = height-(int) round((d.y-loY)*scale);
    if (isElement (highlightF, d.id)) {
      p->setPen (Qt::blue);
      if (doShowNames) p->drawText (xx, yy, d.name);
      p->fillRect (xx-1, yy-1, 4, 4, QBrush(Qt::blue)); 
    }      
    else {
      p->setPen (Qt::black);      
      if (doShowNames) p->drawText (xx, yy, d.name);      
      p->fillRect (xx-1, yy-1, 3, 3, QBrush(d.color)); 
    }      
  }
}


void TmgMapView::savePng (const char* filename) 
{
  try {
    double loX, hiX, loY, hiY;    
    boundingBox (dot, loX, hiX, loY, hiY, false);    
    savePng (dot, loX, hiX, loY, hiY, scale, filename);  
  }
  catch (std::runtime_error& myerr) {
    QMessageBox::critical (this, "savePgm", myerr.what());
  }  
}


void TmgMapView::savePng (const DotList& dot, 
                          double loX, double hiX, double loY, double hiY,
                          double scale, const char* filename) 
{
  int width  = (int) ceil((hiX-loX)*scale);
  int height = (int) ceil((hiY-loY)*scale);  
  QImage img (width, height, 32);
  img.fill (qRgb(0xff,0xff,0xff));    
  for (int i=0; i<dot.size(); i++) {
    const Dot& d = dot[i];      
    if (!finite(d.x) || !finite(d.y)) continue;
    int xx = (int) round((d.x-loX)*scale);
    int yy = height-(int) round((d.y-loY)*scale);
    if (1<=xx && xx<width-1 && 1<=yy && yy<height-1)
      for (int ii=-1;ii<=1;ii++) for (int jj=-1;jj<=1;jj++)
        img.setPixel (xx+ii, yy+jj, d.color.rgb());      
  }
  img.save (filename, "PNG");
}

