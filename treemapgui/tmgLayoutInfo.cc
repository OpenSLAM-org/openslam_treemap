#include "tmgLayoutInfo.h"
#include "tmgTreemapView.h"
#include "treemap/tmNode.h"

TmgLayoutInfo::TmgLayoutInfo ()
        :metric(QFont()), effectiveLevelHeight()
{}


TmgLayoutInfo::TmgLayoutInfo (const QFontMetrics& metric, const QRect& clip)
        :metric(metric), clip(clip), effectiveLevelHeight()
{
    int space = metric.width(QChar(' '));
    int mh  = metric.height();
    int ml = metric.lineSpacing();
    textHeight      =  ml + mh; // 2 lines
    marginX         = 2*space;
    marginY         = metric.lineSpacing()/4; // /2 before
    nodeBoxHeight   = textHeight + 2*marginY;
    edgeSpaceHeight = textHeight;
    baseLevelHeight     = nodeBoxHeight + edgeSpaceHeight;
    subtreeSpace    = 2*space;
    border          = 4*space;
    signSize        = 2*space;
}


QSize TmgLayoutInfo::textSize (char* txt)
{
    char buffer[1000];
    int j=0;
    for (int i=0; txt[i]!='\0' && j<999; i++) if (txt[i]!='[' && txt[i]!=']') {
        buffer[j] = txt[i];
        j++;
    }
    buffer[j]='\0';
    return metric.size (Qt::ShowPrefix, QString(buffer));
}



QSize TmgLayoutInfo::nodeBoxSize (char* txt)
{
    QSize txtSize = textSize (txt);
    return QSize (txtSize.width() + 2*marginX, txtSize.height()+2*marginY);
}


void TmgLayoutInfo::boxes (QRect& nodeBox, QRect& textBox, TmNode* node, int x, int y, char* txt)
{
    //bool isABIB = node==node->getTree()->getABIBNode();
    QSize txtSize = textSize (txt);
    int mX, mY;
    mX = marginX;
    mY = marginY;
    textBox = QRect (x+mX, y+mY, 
                     txtSize.width(), txtSize.height());
    nodeBox =  QRect (x, y,
                      txtSize.width() + 2*mX, txtSize.height()+2*mY);
}


int TmgLayoutInfo::levelHeight (int level) const
{
    if (level<(int) effectiveLevelHeight.size()) return effectiveLevelHeight[level];
    else return baseLevelHeight;
}


void TmgLayoutInfo::increaseLevelHeight (int level, int height)
{
    while (level>=(int) effectiveLevelHeight.size()) 
        effectiveLevelHeight.push_back (baseLevelHeight);
    if (effectiveLevelHeight[level]<height) effectiveLevelHeight[level] = height;
}

