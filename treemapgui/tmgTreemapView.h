//!\author Udo Frese

/*!\file tmgTreemapViewe.h Contains the class \c TmgTreemap being a QT widget that
  shows the treemap tree with represented features.
*/

#ifndef TMGTREEMAPVIEW_H
#define TMGTREEMAPVIEW_H

#include <map>

// SLAM include
#include "treemap/tmTreemap.h"

// QT include
#include <qscrollview.h>
#include <qprinter.h>

// SLAM-GUI includes
#include "tmgLayoutInfo.h"
#include "tmgDecoration.h"


//! Widget that provides a scrollable view the tree of a \c TmTreemap object.
/*! The purpose of this widget is to display the tree of a \c TmTreemap object
  and provide some low level UI for selecting and highlighting nodes. */
class TmgTreemapView : public QScrollView
{
    Q_OBJECT

public:
    //! Standard constructor (like \c QScrollView)
    TmgTreemapView (QWidget* parent=NULL, const char* name=NULL);

    //! Destructor
    ~TmgTreemapView();
    

    //! Sets the \c SlamTreemap object to be displayed.
    /*! If \c takeTreemap is \c true, object takes ownership on \c
        *treemap and frees it in its destructor. Otherwise the object
        holds a pointer to the tree map and the user is responsible
        both for freeing the treemap and for providing that the
        treemap exists while in use by this widget.
    */
    void setTreemap (const TmTreemap* treemap, bool takeTreemap=false);
    
    //! Returns the tree map currently shown.
    const TmTreemap* getTreemap();

    //! Defines how many levels to draw maximally
    void setDrawNrOfLevels (int drawNrOfLevels);    

    //! Returns how many levels should be drawn
    int getDrawNrOfLevels () 
      {return drawNrOfLevels; }
    

    friend class TmgDecoration;
    
    //! Adds a decoration (\c decoration)
    /*! A decoration is some additional information displayed at
        a node or edge. There are different kinds of decorations
        (signs, texts, etc.). \c TmgDecoration is an abstract base
        class of all. 
        See \c TmgDecoration, \c TmgDecorationEdgeMark for concrete
        information on the different decorations. 
        The memory belongs to this object afterward.
    */
    void addDecoration (TmgDecoration* decoration);

    //! Deletes \c decoration
    void deleteDecoration (TmgDecoration* decoration);

    //! Deletes all decorations that have a group name \c groupName.
    void deleteDecorationGroup (const char* groupName);

    //! Deletes all decorations
    void deleteAllDecorations();

    //! Insert decorations describing the different regions
    //! according to \c SlamLandmarkState with base node \c r and \c n.
    /*! The decoration has group name \c name. */
    void decorateRegions (TmNode* r, TmNode* n, char* name);

    //! Prints the treemap view with the parameter specified by \c options to \c file (EPS)
    /*! The return value specifies, whether printing was successfull.
     */
    void print (const char* file, const char* options) throw (runtime_error);

    //! Saves the treemap to a .png file
    void savePng (const char* filename);    

    
public slots:
    //! updates the tree shown in the content are
    void updateContent();

    //! Sets, whether to show internal ids (\c false) or user ids (\c true)
    //! in the leaves.
    void showUserIds (bool doShowUserIds);

    //! Sets whether to show the first status line in each node
    void showStatusLine (bool doShowStatusLine);


protected:
    //! Tree map currently shown
    const TmTreemap* treemap;
    
    //! Whether this object owns \c treemap
    bool ownsTreemap;


    //! A set of decorations mapped by node index.
    typedef multimap <int, TmgDecoration*> DecorationMap;

    //! An iterator for \c DecorationMap
    typedef DecorationMap::iterator DecMI;

    //! Set of decorations mapped to the corresponding node indices
    /*! The memory belongs to this object, so before freeing the map
      all entries must be deleted. */
    DecorationMap decoration;

    //! Whether to show the user ids (\c true) or internal ids (\c false)
    //! in the leave nodes.
    bool doShowUserIds;

    //! Whether to show the first status line in each node
    bool doShowStatusLine;

    //! Whether not to show the lines and slashes
    bool onlyNode;

    //! If \c true, the BIBs are shown, as if they were CIB with a subtree following
    bool doShowBIBAsCIB;

    //! If \c true, eiliminated landmarks are shown with slashed
    /*! Otherwise they are not shown at all. */
    bool doShowEliminated;

    //! Draw only the top \c drawNrOfLevels levels of the tree
    /*! If \c -1 all levels are drawn. */
    int drawNrOfLevels;    



    //! Overloaded \c QScrollView function
    /*! Draws the content of this view, i.e. the tree onto \c p in the region specified
      by \c clipx, \c clipy, \c clipw, \c cliph. */
    virtual void drawContents ( QPainter * p, int clipx, int clipy, int clipw, int cliph );

    //! Removes \c treemap and deletes it if necessary
    void clear();
    
    //! TmgLayoutInfo belongs to \c TmgTreemapView and is thus declared as friend
    friend class TmgLayoutInfo;

    //! Layout parameter to use in \c drawContents() and subroutines.
    /*! Is computed from the actual font size in \c drawContents. */
    TmgLayoutInfo li;

    //! Bounding box of the last tree at call to \c drawContents().
    QRect lastTreeSize;
    
    //! Whether currently rendering to screen, bitmap export or postscript export
    enum {SCREEN, BITMAP, POSTSCRIPT} drawMode;

    //! Check and Modify the size of the content area if necessary
    /*! The routine looks at \c lastTreeSize and checks, whether
        it is to small or large for the current content area. If
        this is the case, the content area is enlarged and or shinked
        and \c updateGeometry() is called.
        
        The size of the content area changes in step of 1/3 of the
        visible size.
     */
    void checkContentSize();

    //! Internal subroutine for \c drawContents()
    /*! Recursively draws the subtree below \c node onto \p using the
        layout info in \c li.  The tree is aligned, so that the left
        upper corner of its bounding box is \c x, \c y. The size of the bounding
        box is returned in \c size. In \c attach the point where to connect to
        the box representing \c node is returned.

        The worst case path is shown in red by passing \c isOnWCPath.

        The layout info is derived from the font metrics by \c TmgLayoutInfo::TmgLayoutInfo.
     */
    void recursiveDraw (QPainter* p, TmNode* node, int x, int y, QSize& size, QPoint& attach, int level, bool isOnWCPath);

    //! Recursively updates the level heights in \c li considering the tree below \c node in \c level 
    void recursiveComputeLevelHeight (TmNode* node, int level);
    
    //! Like \c drawContents but with only optionally checking whether the size has changed.
    void drawAll ( QPainter * p, int clipx, int clipy, int clipw, int cliph, bool checkCS=false );

    //! Draws a single node
    /*! The node is shown as a rectangle with two lines of text. The first line contains the
        character 'R', 'I', 'O' if the RL, IIB or OIB resp. is valid. The second line contains
        a list of all landmarks represented at this node. In the case of a leaf, the landmarks
        only represented in the BIB are included but underlined. 
        If the text has already been computed it can be passed in \c txt. Otherwise it
        is computed from \c node.
    */
    void drawNode (QPainter* p, TmNode* node, int x, int y, char* txt=NULL);

    //! Draws all decorations for \c node.
    /*! \c nodeBox is the rectangle depicting the node, \c attachPChild0, \c attachPChild1
        and \c attachPParent are the attachment points for the edge to child 0, child 1
        and parent respectively (lyinng on the border of \c nodeBox).
        NOTE: It only paints the node decorations using \c TmgDecoration::drawNode, not the
        corresponding edge decorations.
    */
    void drawNodeDecorations (TmNode* node, QPainter* p, const QRect& nodeBox,
                              const QPoint& attachPChild0, const QPoint& attachPChild1, 
                              const QPoint& attachPParent);

    //! Computes the extra width needed by the decorations for node \c node
    /*! The routine proceeds by computing the maximum for the extra width of all
        decorations for node \c node. \c p is the painter used. The result is
        returned in \c left and \c right.
     */
    void decorationsExtraWidth (TmNode* node, QPainter* p, int& left, int& right);

    //! Draws a single edge onto \c p
    /*! The edge connects \c node to its parent. \c nodeAttachmentPoint is the attachment point on
        the icon for \c node. \c parentAttachmentPoint is the attachment point on the icon
        for the parent of node.
        It is possible to call this routine for the root node.
        If some edges are highlighted (see \c setTransferHighlight()) the highlight is drawn.
        if \c bold is \c true, the edge is drawn with width 2.
     */
    void drawEdge (QPainter* p, TmNode* node, const QPoint& nodeAttachmentPoint, const QPoint& parentAttachmentPoint, bool bold);

    //! Draws all decorations for the edge from \c node to its parent onto \p
    /*! \c attach is the attachent point of the edge on \c node, \c parentAttach on its parent.
        NOTE: It only paints the edge decorations using \c TmgDecoration::drawEdge, not the
        corresponding node decorations.
    */      
    void drawEdgeDecorations (TmNode* node, QPainter* p, const QPoint& attach, const QPoint& parentAttach);


    //! Returns the brush to be used for filling the box of node \c node.
    /*! If a \c TmgDecorationNodeBrush object is specified for that node,
        its brush is used. Otherwise \c defaultBrush (normally no filling) is
        returned. */
    QBrush getBackgroundBrush (TmNode* node, const QBrush& defaultBrush=QBrush());

    //! Computes the description text for \c node.
    /*! See \c drawNode() for a documentation on the text. If the text is longer
      than \c maxLen (including '\0') is is truncated. */
    void computeNodeText (char* text, int maxLen, TmNode* node);

    //! Breaks \c text by replacing ' ' with '\n' appropriately
    void lineBreak (char* text);    

    //! parses the options used for printing
    void printOptions (const char* options, QFont& font);
    

    //! Returns the extension of filename \c file.
    /*! The returned pointer points to the memory as \c file. No copy is performed.
     */
    const char* extension(const char* file);

    // Like \c p->drawText (box, txt) but strikes out text in "[]'
    void drawTextWithSlashes (QPainter* p, const QRect& box, char* txt);

    //! Copies \c file1 to \c file2 defining a bounding box
    /*! \c box is the bounding box is coordinates as used by \c QPainter. \c pageBox
        is the box of the whole page in coordinates as used by \c QPainter.
    */
    void setBoundingBox (const char* file1, const char* file2, const QRect& box, const QRect& pageBox);

    //! Prints the textual name for feature \c il[i] into \c text[ctr..]
    /*! \c ctr and i are all incremented. Since a textutal description can cover
        several features \c i can be increased by more than 1
    */
    void addFeatureText (char* txt, int& ctr, TmExtendedFeatureList& il, int& i);    
private:
    
    //! Internal routine for parsing options
    bool parse (const char* options, int& idx, const char* token);

    void checkValidity (char* txt);    
};

#endif
