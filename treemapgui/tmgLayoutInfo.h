#ifndef TMGLAYOUTINFO_H
#define TMGLAYOUTINFO_H

using namespace std;

#include <xycontainer/xycVector.h>

// QT classes
#include <qrect.h>
#include <qfontmetrics.h>

// SLAM classes
class TmNode;

//! Parameter for the graph layout in \c SgTreeMapView computed from
//! the font metrics
class TmgLayoutInfo
{
public:
    //! Font metrics on which this computation is based
    QFontMetrics metric;
    
    //! Clipping rectangle
    QRect clip;
    
    //! Text height (=2 lines)
    int textHeight;

    //! Horizontal margin (between rectangle and text)
    int marginX;

    //! Vertical margin (between rectangle and text)
    int marginY;

    //! Height of a node box
    int nodeBoxHeight;

    //! Height of the space for edge lines
    int edgeSpaceHeight;

    //! Minimal height of one level (node box plus space for edges)
    /*! Use \c levelHeight() instead. It is at least \c
      baseLevelHeight. If a node in the level contains more entries than
      fit with \c baseLevelHeight, \c levelHeight(level) is
      increased accordingly. */
    int baseLevelHeight;

    /*! Height of the different levels in the tree. */
    /* Is always at least \c levelHeight but can be increased by calling
       \c increaseLevelHeight(). Do not access directly, call
       \c levelHeight() instead.
    */
    XycVector<int> effectiveLevelHeight;
    
    //! Spacing between two subtrees of level i is (i-1)*subtreeSpace
    int subtreeSpace;
    
    //! Border between the tree and the bounding box
    int border;
    
    //! Size of the insert / transfer / delete signs
    int signSize;
    
    //! Uninitialised Layout
    TmgLayoutInfo();
    
    //! A Layout parameter set for font \c metric with clipping rectangle \c clip
    TmgLayoutInfo (const QFontMetrics& metric, const QRect& clip);
    
    //! Returns the size of the node box with \c txt printed inside.
    QSize nodeBoxSize (char* txt);
    
    //! Returns the surrounding box and text box of a node
    /*! The node box is returned in \c nodeBox and the text box in
      \c textBox assuming, that the node shall be drawn at \c
      x,y (left / upper).  The text must have been already computed
      an is passed in \c txt.
    */
    void boxes (QRect& nodeBox, QRect& textBox, TmNode* node, int x, int y, char* txt);

    //! Returns the size of \c txt when printed
    /*! '[' and ']' are removed, since they are indicator for striking
         out the text between the brackets. So they are not counted.
     */
    QSize textSize (char* txt);

    //! Returns the height of tree \c level
    /*! Normally it returns \c levelHeight, except, if a larger
        value has been set by \c increaseLevelHeight. */
    int levelHeight (int level) const;

    //! Increases the height of \c level to \c height
    /*! If it is already larger, nothing is changed. */
    void increaseLevelHeight (int level, int height);
};

#endif
