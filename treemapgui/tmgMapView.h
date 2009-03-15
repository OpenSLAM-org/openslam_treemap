//!\author Udo Frese

/*!\file tmgMapView.h Contains the class \c TmgMapView being a QT widget that
  shows a 2D map estimate consisting of dots.
*/

#ifndef TMGMAPVIEW_H
#define TMGMAPVIEW_H

#include <xycontainer/xycVector.h>

// QT include
#include <qscrollview.h>
#include <qprinter.h>
#include <qstring.h>

//! A widget that shows a map of dots, potentially with names
class TmgMapView : public QScrollView 
{
  Q_OBJECT

 public:
  //! Std constructor
  TmgMapView (QWidget* parent=NULL, const char* name=NULL);

  //! Returns to empty state with to map shown
  void clear();  

  class Dot 
  {
  public:
    //! Number identifying the dot.
    int id;    
    //! Position
    double x, y;
    //! name
    QString name;
    //! Color
    QColor color;    

    Dot ():id(-1), x(0), y(0), name(), color() {}
    Dot (double x, double y) :id(-1), x(x), y(y), name(), color() {}
    Dot (double x, double y, char* str) :id(-1), x(x), y(y), name(str), color() {}
    Dot (int id, double x, double y, char* str) :id(id), x(x), y(y), name(str), color() {}	  
      Dot (int id, double x, double y, char* str, const QColor& color) :id(id), x(x), y(y), name(str), color(color) {}	  
  };  

  typedef XycVector<Dot> DotList;  
    
  //! Sets the list of dots shown
  void setDots (const DotList& dot);

  //! Removes all dots so nothing is shown
  void clearDots ();
  
  //! Adds one dot to the list of dots to be shown
  void addDot (const Dot& newDot);  
  
  //! Sets the area shown on the screen
  /*! If \c enlarge is true, the area is enlarged so that it covers all dots. */
  void setArea (double loX, double hiX, double loY, double hiY, bool enlarge=true); 

  //! Enlarges the area such that it contains all dots
  void enlargeArea ();  

  //! Sets the scale (pixel per unit in the estimate)
  void setScale (double scale); 

  //! Show these features as highlighted
  /*! Those features are highlighted for which \c .id is in \c fl. */
  void highlightFeature (const XycVector<int>& fl);

  //! Ensures, that \c x+/-margin, y+/-margin is visible
  void ensureVisible (double x, double y, double margin);
  
  //! Returns the current scale (pixel/m)
  double getScale () const  {return scale;}

  //! Saves \c dots with \c scale as a .png
  static void savePng (const DotList& dots, 
                       double loX, double hiX, double loY, double hiY,
                       double scale, const char* filename);  
  
  //! Saves the current map with the current settings
  void savePng (const char* filename);  

  //! Saves the current map with the current settings
  void savePng (const QString& filename) {savePng((const char*)filename);}
  



  public slots:
  //! Repaints the map
  void updateContent ();
  
  //! Resizes the scroll view
  void updateContentSize ();  
  
  //! Whether to show landmark names
  void setDoShowNames (bool doShowNames);  
  
 protected:
  //! The list of dots shown
  XycVector<Dot> dot;
  
  //! Whether to show landmark names
  bool doShowNames;  

  //! Scale: 1 in the estimate is \c scale pixel on screen
  double scale;  

  //! Area shown on screen
  double loX, hiX, loY, hiY;  


  //! Show these features as highlighted
  XycVector<int> highlightF;  


  //! Overloaded \c QScrollView function
  /*! Draws the content of this view, i.e. the building onto \c p in the region specified
    by \c clipx, \c clipy, \c clipw, \c cliph. */
  virtual void drawContents ( QPainter * p, int clipx, int clipy, int clipw, int cliph );

  //! Draws the map on \c  p
  void draw (QPainter* p);  

  static void boundingBox (const DotList& dot,
                           double& loX, double& hiX, double& loY, double& hiY, bool enlarge);  

  //! Returns, whether \c id is contained in \c list
  /*! list must be sorted. */
  static bool isElement (const XycVector<int>& list, int id);  
};


#endif
