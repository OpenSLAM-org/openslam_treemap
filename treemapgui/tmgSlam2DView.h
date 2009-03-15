//!\author Udo Frese

/*!\file tmgSlam2DView.h Contains the class \c TmgSlam2DView being a 
   QT widget that shows a treemap 2D (\c TmSlamDriver2DL or \c TmSlamDriver2DP)
   estimate as a 2D landmark map.
*/

#include "tmgMapView.h"
#include <treemap/tmSlamDriver2DL.h>
#include <treemap/tmSlamDriver2DP.h>

//! Class for showing the map estimated by a \c TmSlamDriver2DL or \c TmSlamDriver2DP object
/*! This is basically a generic \c TmgMapView plus routines for converting
  the \c TmSlamDriver2DL or \c TmSlamDriver2DP estimate to a \c XycVector<Dot>.
*/
class TmgSlam2DView : public TmgMapView
{
 public:
  //! Std constructor
  TmgSlam2DView (QWidget* parent=NULL, const char* name=NULL);

  //! Produces a dots from a \c level of a \c TmSlamDriver2DL map
  static void convertTmSlamDriverToDots (const TmSlamDriver2DL& treemap, int level, DotList& dot);  
  
  //! Show \c level of \c treemap
  /*! Does not store the treemap. */
  void setDotsFromTreemap (const TmSlamDriver2DL& treemap, int level);

  //! Saves \c level of \c treemap with \c scale as a .png
  static void savePng (const TmSlamDriver2DL& treemap, int level, 
                       double loX, double hiX, double loY, double hiY,
                       double scale, const char* filename);  

  //! Produces a dots from a \c level of a \c TmSlamDriver2DP map
  static void convertTmSlamDriverToDots (const TmSlamDriver2DP& treemap, DotList& dot);  
  
  //! Show \c level of \c treemap
  /*! Does not store the treemap. */
  void setDotsFromTreemap (const TmSlamDriver2DP& treemap);

  //! Saves \c level of \c treemap with \c scale as a .png
  static void savePng (const TmSlamDriver2DP& treemap, 
                       double loX, double hiX, double loY, double hiY,
                       double scale, const char* filename);  
  
  //! Defines, which features are highlighted
  void highlightFeature (const TmExtendedFeatureList& fl);  
};

