#ifndef TMGDECORATION_H
#define TMGDECORATION_H

#include <xycontainer/xycVector.h>

#include <qpen.h>
#include <qbrush.h>
#include <qcolor.h>

//! QT Forward
class QPainter;
//! QT Forward
class QRect;
//! QT Forward
class QPoint;

//! SLAM Forward
class TmNode;
class TmTreemap;

//! SLAM-GUI Forward
class TmgTreemapView;

//! SLAM-GUI forward
class TmgLayoutInfo;

using namespace std;

//! Base class for decorations (additional information displayed at a
//! node in a \c TmgTreeMapView object.
class TmgDecoration
{
public:
    //! Direct access for TmgTreeMapView
    friend class TmgTreemapView;

    //! Copy constructor
    TmgDecoration (const TmgDecoration& dec);
    
    //! Destructor
    virtual ~TmgDecoration();

    //! Returns a duplicate of \c this.
    /*! This routine \em must be overloaded by each derived class,
        since the inherited implementation generates an object of the
        inherited class.  As there may be no \c TmgDecoration
        instances, the current implementation throws an error.
    */
    virtual TmgDecoration* duplicate();

    //! Draw the node decoration for \c node onto \c p.
    /*! Virtual function to be overloaded by an inheritant that draws
        the decoration for node \c node onto \c p.  \c li is an \c
        TmgTreeMapView::LayoutInfo object containing information about
        the tree layout. \c nodeBox is the box of the node (where the
        node text is written info).  \c attachPChild0, \c
        attachPChild1, \c attachPParent are the attachment points of
        the edges to child 0, child 1, and the parent
        respectively. All three lie on the border of \c nodeBox.
     */
    virtual void drawNode (QPainter* p, const TmgLayoutInfo& li, const QRect& nodeBox,
                           const QPoint& attachPChild0, const QPoint& attachPChild1, const QPoint& attachPParent);

    //! Draw the edge decoration for the edge from \c node to its
    //! parent to \c p.
    /*! Virtual function to be overloaded by an inheritant that draws
        the decoration for the edge from \c node to its parent onto
        \c p.  \c li is an \c TmgTreeMapView::LayoutInfo object
        containing information about the tree layout. \c attach is the
        attachment point of the edge line at \c node and \c
        parentAttach at its parent. */
    virtual void drawEdge (QPainter* p, const TmgLayoutInfo& li, const QPoint& attach, const QPoint& parentAttach);


    //! Returns the fill brush, this decoration defines for the
    //! background of the node rectangle
    /*! Normally the background of a node is unfilled, thus showing
        the windows background color.  A decoration can define an
        other background. In this case, it must overload this method
        to return \c true and set \c brush to the desired filling
        brush. If \c false is returned, which is the \c TmgDecoration
        implementation, than the background color is not affected.
    */
    virtual bool nodeFillBrush (QBrush& brush);

    //! Returns extra space needed by the decoration
    /*! \c left is the space needed left of the node box, \c is the
        space needed right of the node box.  The default implementation
        returns 0, 0. */
    virtual void extraWidth (QPainter* p, const TmgLayoutInfo& li, int& left, int& right);

    //! Changes the node to the corresponding node of \c treeMap base
    //! on \c node->index.
    /*! If the corresponding node does not exist in \c treeMap, the
        node is set to NULL, making the decoration effectively
        void. (\c TmgTreeMapView can cope with such decorations).
    */
    void toTreeMap (const TmTreemap& treeMap);

protected:
    //! Pointer to the \c TmgTreeMapView object this decoration belongs to.
    /*! May be \c NULL, if the decoration object is used for saving
        decorations and not currently displayed in a \c TmgTreeMapView
        widget. */
    TmgTreemapView* view;
    //! Pointer to the node to be decorated
    TmNode* node;
    //! String identifying a group of decorations.
    char* group;
    //! Whether the memory of \c group belongs to this object
    bool ownsGroup;


    //! Constructor: Only to be used by an inheriting class.
    TmgDecoration (TmNode* node, const char* group);

    //! Sets \c group to a \c group
    /*! \c group is set to a copy of \c group and \c ownsGroup set.
     */
    void setGroup (const char* group);
};

//! List of (polymorphic) decorations.
/*! We define, that the memory of the decoration belongs to the owner
    of the vector, but deallocation is not done automatically so the
    owner is responsible for that.
*/
typedef XycVector<TmgDecoration*> TmgDecorationList;

//! Decoration that puts a mark on an edge
class TmgDecorationEdgeMark : public TmgDecoration
{
public:
    //! Different types of marks
    enum MarkType {NONE=0, CIRCLE=1, DOWNWARDBRACKET=2, UPWARDBRACKET=3, CROSS=4, EDGE=5};
    

    //! Constructor for a mark on the edge from \c node to its parent.
    /*! See \c TmgDecoration::TmgDecoration \c node and \c group. \c mark
      is the type of mark, \c pen the pen (color, dashing, etc.) with which the mark
      is drawn. \c txt is written beside the mark (if not \c NULL). */
    TmgDecorationEdgeMark (TmNode* node, const char* group, MarkType mark, const QPen& pen, const char* txt=NULL);

    //! Overloaded \c TmgDecoration function
    virtual TmgDecoration* duplicate();

    //! Overloaded 'TmgDecoration' function that draws the decoration
    virtual void drawEdge (QPainter* p, const TmgLayoutInfo& li, const QPoint& attach, const QPoint& parentAttach);

    //! Overloaded destructor
    virtual ~TmgDecorationEdgeMark();

protected:
    //! Type of the mark
    MarkType mark;

    //! Pen to use in drawing the mark
    QPen pen;

    //! Text to write. Memory is owned by the object.
    const char* txt;
};


//! A node decoration redefining the background of a node
class TmgDecorationNodeBrush : public TmgDecoration
{
public:
    //! Constructs a decoration that sets the background of
    //! node \c node to \c brush. 
    /*! See \c TmgDecoration::TmgDecoration for \c group. */
    TmgDecorationNodeBrush(TmNode* node, const char* group, const QBrush& brush);

    //! Overloaded \c TmgDecoration method.
    virtual TmgDecoration* duplicate();

    //! Overloaded \c TmgDecoration method.
    virtual bool nodeFillBrush (QBrush& brush);
protected:
    QBrush _brush;
};

//! Decoration of a node: colored box
class TmgDecorationNodeMark : public TmgDecoration
{
public:
    //! Decoration with width \c width and color \c color
    TmgDecorationNodeMark (TmNode* node, const char* groupName, const QPen& color);
    
    //! Overloaded \c SgDecoration function
    virtual TmgDecoration* duplicate();

    //! Overloaded function that draws the node decoration
    virtual void drawNode (QPainter* p, const TmgLayoutInfo& li, const QRect& nodeBox,
                           const QPoint& attachPChild0, const QPoint& attachPChild1, const QPoint& attachPParent);
protected:
    QPen color;
};


#endif
