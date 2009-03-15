#include "tmgDecoration.h"
#include "tmgLayoutInfo.h"

#include <qpainter.h>
#include <qpoint.h>

#include "treemap/tmNode.h"
#include "treemap/tmTreemap.h"

TmgDecoration::~TmgDecoration()
{
    if (ownsGroup && group!=NULL) delete group;
}


TmgDecoration* TmgDecoration::duplicate()
{
    assert(false);
    return NULL;
}



void TmgDecoration::drawNode (QPainter* p, const TmgLayoutInfo& li, const QRect& nodeBox,
                         const QPoint& attachPChild0, const QPoint& attachPChild1, const QPoint& attachPParent)
{
    // Empty: To be overloaded by inheritant
}


void TmgDecoration::drawEdge (QPainter* p, const TmgLayoutInfo& li, const QPoint& attach, const QPoint& parentAttach)
{
    // Empty: To be overloaded by inheritant
}


bool TmgDecoration::nodeFillBrush (QBrush& brush)
{
    // Empty: To be overloaded by inheritant
    return false;
}


void TmgDecoration::setGroup (const char* myGroup)
{
    if (ownsGroup && group!=NULL) delete group;
    group = NULL;
    if (group!=NULL) {
        group = new char[strlen(myGroup)+1];
        strcpy (group, myGroup);
        ownsGroup = true;
    }
}


TmgDecoration::TmgDecoration (TmNode* node, const char* group)
        :view (NULL), node(node), group(NULL), ownsGroup(false)
{
    setGroup (group);
}


TmgDecoration::TmgDecoration (const TmgDecoration& dec)
        :view(dec.view), node(dec.node), group(NULL), ownsGroup(false)
{
    setGroup (dec.group);
}

void TmgDecoration::extraWidth (QPainter* p, const TmgLayoutInfo& li, int& left, int& right)
{
    left = right = 0;
}


void TmgDecoration::toTreeMap (const TmTreemap& treemap)
{
  if (node!=NULL) node = treemap.getNode (node->index);
}

// ***** TmgDecorationEdgeMark
TmgDecorationEdgeMark::TmgDecorationEdgeMark (TmNode* node, const char* group, MarkType mark, const QPen& pen, const char* txt)
        : TmgDecoration (node, group), mark(mark), pen(pen), txt(txt)
{}


TmgDecoration* TmgDecorationEdgeMark::duplicate()
{
    return new TmgDecorationEdgeMark(*this);
}


void TmgDecorationEdgeMark::drawEdge (QPainter* p, const TmgLayoutInfo& li, const QPoint& attach, const QPoint& parentAttach)
{
    QPoint pC ((attach.x()+parentAttach.x())/2, (attach.y()+parentAttach.y())/2);
    QRect sR (pC.x()-li.signSize, pC.y()-li.signSize, 2*li.signSize, 2*li.signSize);
    p->save();
    p->setPen (pen);
    switch (mark) {
        case NONE:
            break;
        case CIRCLE:
            p->drawEllipse (sR);
            break;
        case DOWNWARDBRACKET:
            p->drawLine (sR.bottomLeft(), sR.topLeft());
            p->drawLine (sR.topLeft(), sR.topRight());
            p->drawLine (sR.topRight(), sR.bottomRight());
            break;
        case UPWARDBRACKET:
            p->drawLine (sR.bottomLeft(), sR.topLeft());
            p->drawLine (sR.bottomLeft(), sR.bottomRight());
            p->drawLine (sR.topRight(), sR.bottomRight());
            break;
        case CROSS:
            p->drawLine (sR.topLeft(), sR.bottomRight());
            p->drawLine (sR.bottomLeft(), sR.topRight());
            break;
        case EDGE:
            p->drawLine (attach, parentAttach);
            break;
        default:
            assert (false); // illegal mark type
    }
    if (txt!=NULL) {
        QRect textRect (sR.right()+2, sR.top(), li.metric.width(txt), sR.height());
        p->drawText (textRect, 0, txt);
    }
                        
    p->restore();
}


TmgDecorationEdgeMark::~TmgDecorationEdgeMark()
{
    if (txt!=NULL) delete txt;
}

// TmgDecorationNodeBrush

TmgDecorationNodeBrush::TmgDecorationNodeBrush (TmNode* node, const char* group, const QBrush& brush)
        :TmgDecoration (node, group), _brush(brush)
{}


TmgDecoration* TmgDecorationNodeBrush::duplicate()
{
    return new TmgDecorationNodeBrush::TmgDecorationNodeBrush(*this);
}


bool TmgDecorationNodeBrush::nodeFillBrush (QBrush& brush)
{
    brush = _brush;
    return true;
}

// TmgDecorationNodeMark



TmgDecorationNodeMark::TmgDecorationNodeMark (TmNode* node, const char* groupName, const QPen& color)
        :TmgDecoration(node, groupName), color(color)
{}

    
TmgDecoration* TmgDecorationNodeMark::duplicate()
{
    return new TmgDecorationNodeMark(*this);
}


void TmgDecorationNodeMark::drawNode (QPainter* p, const TmgLayoutInfo& li, const QRect& nodeBox,
                                     const QPoint& attachPChild0, const QPoint& attachPChild1, const QPoint& attachPParent)
{
    p->save();
    p->setPen (color);
    int w2 = (color.width()+1)/2;
    p->drawRect (nodeBox.left()-w2, nodeBox.top()-w2, nodeBox.width()+2*w2, nodeBox.height()+2*w2);
    p->restore();
}
