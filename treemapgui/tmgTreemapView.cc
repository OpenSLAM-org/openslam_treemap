#include "tmgTreemapView.h"
#include "tmgDecoration.h"

#include <qmessagebox.h>
#include <qstrlist.h>
#include <qimage.h>
#include <qpaintdevicemetrics.h>
#include <qfontmetrics.h>
#include <qpainter.h>

#include <math.h>

TmgTreemapView::TmgTreemapView (QWidget* parent, const char* name)
        :QScrollView (parent, name), treemap(NULL), ownsTreemap(false), decoration(),
         doShowUserIds(false), doShowStatusLine(true), onlyNode(false), doShowBIBAsCIB(false),
         doShowEliminated(true), drawNrOfLevels (-1),
         drawMode(SCREEN)
{
    viewport()->setBackgroundMode( QWidget::PaletteBase );
}


TmgTreemapView::~TmgTreemapView()
{
    clear();
}


void TmgTreemapView::clear()
{
    deleteAllDecorations();
    if (ownsTreemap && treemap!=NULL) delete treemap;
    treemap = NULL;
    ownsTreemap = false;
}


void TmgTreemapView::setTreemap (const TmTreemap* treemap, bool takeTreemap)
{
  if (this->treemap!=treemap) {
    clear();
    this->treemap = treemap;
    ownsTreemap = takeTreemap && (treemap!=NULL);
  }
  emit updateContent();
}


const TmTreemap* TmgTreemapView::getTreemap()
{
    return treemap;
}


void TmgTreemapView::setDrawNrOfLevels (int drawNrOfLevels)
{
  this->drawNrOfLevels = drawNrOfLevels;
  updateContent();
}


void TmgTreemapView::showUserIds (bool doShowUserIds)
{
    doShowUserIds = doShowUserIds;
    updateContent();
}


void TmgTreemapView::showStatusLine (bool doShowStatusLine)
{
    doShowStatusLine = doShowStatusLine;
    updateContent();
}



void TmgTreemapView::updateContent()
{
    emit viewport()->update();
}


void TmgTreemapView::drawContents ( QPainter * p, int clipx, int clipy, int clipw, int cliph )
{
    drawAll (p, clipx, clipy, clipw, cliph, true);
}


void TmgTreemapView::drawAll (QPainter* p, int clipx, int clipy, int clipw, int cliph, bool checkCS )
{
    if (treemap==NULL || treemap->root==NULL) return;
    li = TmgLayoutInfo(p->fontMetrics(), QRect(clipx, clipy, clipw, cliph));
    recursiveComputeLevelHeight (treemap->root, 0);
    QSize size;
    QPoint attachment;
    int startH;
    if (onlyNode) startH = 0;
    else startH = li.edgeSpaceHeight;
    recursiveDraw (p, treemap->root, li.border, li.border+startH, size, attachment, 0, true);
    QPoint parentAttach (attachment.x(), li.border);
    if (!onlyNode) {
        drawEdge (p, treemap->root, attachment, parentAttach, true);
        drawEdgeDecorations (treemap->root, p, attachment, parentAttach);
    }
    if (checkCS) {
        lastTreeSize = QRect (li.border, li.border, size.width(), size.height()+startH);
        checkContentSize();
    }
}



void TmgTreemapView::checkContentSize()
{
    QSize newSize = QSize(contentsWidth(), contentsHeight());
    int stepX = visibleWidth()/3;
    if (stepX<50) stepX = 50;
    int stepY = visibleHeight()/3;
    if (stepY<50) stepY = 50;
    bool changed=false;
    while (lastTreeSize.right()>newSize.width()) {
        newSize.setWidth (newSize.width()+stepX);
        changed = true;
    }
    while (lastTreeSize.right()<newSize.width()-stepX) {
        newSize.setWidth (newSize.width()-stepX);
        changed =true;
    }
    while (lastTreeSize.bottom()>newSize.height()) {
        newSize.setHeight (newSize.height()+stepY);
        changed = true;
    }
    while (lastTreeSize.bottom()<newSize.height()-stepY) {
        newSize.setHeight (newSize.height()-stepY);
        changed =true;
    }
    if (changed) {
        resizeContents (newSize.width(), newSize.height());
        updateGeometry();
    }   
}


void TmgTreemapView::drawEdge (QPainter* p, TmNode* node, const QPoint& nodeAttachmentPoint, const QPoint& parentAttachmentPoint, bool isOnWCPath)
{
  if (isOnWCPath) {
    int w = 3;
    p->save();    
    p->setPen (QPen (p->pen().color(), w));
    p->drawLine (nodeAttachmentPoint, parentAttachmentPoint);
    p->restore();    
  }
  else p->drawLine (nodeAttachmentPoint, parentAttachmentPoint);
}


void TmgTreemapView::drawEdgeDecorations (TmNode* node, QPainter* p, const QPoint& attach, const QPoint& parentAttach)
{
    pair<DecMI, DecMI> range = decoration.equal_range (node->index);
    for (DecMI it=range.first; it!=range.second; it++) {
        (*it).second->drawEdge (p, li, attach, parentAttach);
    }
}


QBrush TmgTreemapView::getBackgroundBrush (TmNode* node, const QBrush& defaultBrush)
{
    QBrush brush(defaultBrush);
    pair<DecMI, DecMI> range = decoration.equal_range (node->index);
    for (DecMI it=range.first; it!=range.second; it++) {
        (*it).second->nodeFillBrush (brush);
    }
    return brush;
}


void TmgTreemapView::recursiveComputeLevelHeight (TmNode* node, int level)
{
    if (node==NULL || drawNrOfLevels>=0 && level>drawNrOfLevels) return;
    char text[1000];
    computeNodeText (text, 1000, node);
    QSize nodeBoxSize = li.nodeBoxSize(text);
    li.increaseLevelHeight (level, nodeBoxSize.height() + li.edgeSpaceHeight);
    recursiveComputeLevelHeight (node->child[0], level+1);
    recursiveComputeLevelHeight (node->child[1], level+1);    
}


void TmgTreemapView::recursiveDraw (QPainter* p, TmNode* node, int x, int y, QSize& size, QPoint& attach, int level, bool isOnWCPath)
{
    if (node==NULL || drawNrOfLevels>=0 && level>drawNrOfLevels) {
        size = QSize(0,0);
        attach = QPoint();
        return;
    }

    bool isMaxLevel = (drawNrOfLevels>=0 && level==drawNrOfLevels);    
    char text[1000];
    computeNodeText (text, 1000, node);
    QSize nodeBoxSize = li.nodeBoxSize(text);
    int leftExtra, rightExtra;

    QSize subtreeSize; // Size of children including lines to them

    int lvlH = li.levelHeight (level);
    if (node->isLeaf() || isMaxLevel) {
        decorationsExtraWidth (node, p, leftExtra, rightExtra);
        attach = QPoint(x+leftExtra+nodeBoxSize.width()/2, y);
        drawNode (p, node, x+leftExtra, y, text);
        if (doShowBIBAsCIB || !node->isLeaf() && isMaxLevel) {
            // Edges at fake CIBs
          int w5 = nodeBoxSize.width()/5;          
          QPoint attach0 (x+leftExtra+2*w5, y+nodeBoxSize.height());
          QPoint c0 (x+leftExtra+w5, y+li.levelHeight (level));
          p->drawLine (attach0, c0);
            
          QPoint attach1 (x+leftExtra+3*w5, y+nodeBoxSize.height());
          QPoint c1 (x+leftExtra+4*w5, y+li.levelHeight (level));
          p->drawLine (attach1, c1);
        }
        else lvlH = nodeBoxSize.height();
        size = QSize(nodeBoxSize.width()+leftExtra+rightExtra, lvlH);
        QRect nodeBox (x+leftExtra, y, nodeBoxSize.width(), nodeBoxSize.height());
        drawNodeDecorations (node, p, nodeBox, QPoint(), QPoint(), attach);
        return;
    }

    // First draw the children and compute their size
    QSize size0, size1;
    QPoint attach0, attach1; // Attachment point for the children
    // First draw child 0
    bool isOnWCPath0 = isOnWCPath && node->child[0]->worstCaseUpdateCost>=node->child[1]->worstCaseUpdateCost;    
    recursiveDraw (p, node->child[0], x, y+lvlH, size0, attach0, level+1, isOnWCPath0);
    // Now draw child 1
    int space = (node->getHeight()-1)*li.subtreeSpace;
    int newX = x+size0.width() + space;
    bool isOnWCPath1 = isOnWCPath && node->child[0]->worstCaseUpdateCost<=node->child[1]->worstCaseUpdateCost;    
    recursiveDraw (p, node->child[1], newX, y+lvlH, size1, attach1, level+1, isOnWCPath1);
    subtreeSize = QSize (size0.width() + space + size1.width(), 
                         max(size0.height(), size1.height()));

    // Now find where to place the node box relative to the subtrees
    int xN  = (attach0.x() + attach1.x())/2;
    int wN  = nodeBoxSize.width();
    int sW = subtreeSize.width();
    if (wN>sW || xN-wN/2<x) xN = x + wN/2; // align left
    else if (xN+wN/2>x+sW)  xN = x+sW -wN/2; // align right
    
    QRect nodeBox (xN-wN/2, y, nodeBoxSize.width(), nodeBoxSize.height());
    
    drawNode (p, node, xN-wN/2, y, text);

    attach = QPoint (xN, y);
    size = QSize (max(sW, wN), lvlH+subtreeSize.height());

    int xC0 = xN-wN/2+2*li.marginX;
    int xC1 = xN+wN-wN/2-2*li.marginX;
    if (xC0>xC1) xC0 = xC1 = xN;
    QPoint attachC0 (xC0, y+nodeBoxSize.height());
    QPoint attachC1 (xC1, y+nodeBoxSize.height());
    if (!onlyNode) {
        drawEdge (p, node->child[0], attach0, attachC0, isOnWCPath0);
        drawEdgeDecorations (node->child[0], p, attach0, attachC0);
        drawEdge (p, node->child[1], attach1, attachC1, isOnWCPath1);
        drawEdgeDecorations (node->child[1], p, attach1, attachC1);
    }
    
    drawNodeDecorations (node, p, nodeBox, attachC0, attachC1, attach);
}


void TmgTreemapView::drawTextWithSlashes (QPainter* p, const QRect& box, char* txt)
{
    if (txt==NULL) return;
    QFontMetrics metrics = p->fontMetrics();

    double perc = 0.75;
    int slashAscent  = metrics.ascent() - (int) ((1-perc)*metrics.height());
    int slashDescent = metrics.descent(); // - (int) ((1-perc)/2*metrics.height());
    
    int y = box.top() + metrics.ascent(), x=box.left();
    bool isFirstChar=true;
    enum {NOSLASH, SLASHING, REPRINTING} slashMode = NOSLASH;
    int xSlash = 0, ySlash = 0, slashPos=0;
    for (int i=0; txt[i]!='\0'; i++) {
        if (txt[i]=='[') { // Start  slashing
            if (slashMode!=NOSLASH) return;
            slashMode = SLASHING;
            xSlash = x;
            ySlash = y;
            slashPos=i;
        }
        else if (txt[i]==']') { // Stop slashing
            if (slashMode==SLASHING) {
                p->save();
                p->setPen (Qt::black);
                if (!onlyNode) 
                    p->drawLine (xSlash-1, ySlash+slashDescent, x+1, y-slashAscent);
                p->restore();
                x = xSlash;
                y = ySlash;
                i = slashPos;
                slashMode = REPRINTING;
            }
            else if (slashMode==REPRINTING) slashMode = NOSLASH;
            else return;
        }
        else if (txt[i]=='\n') { // newline
            if (slashMode!=NOSLASH) return;
            x = box.left();
            y += metrics.lineSpacing();
            isFirstChar = true;
        }
        else { // Print one char
            QChar ch(txt[i]);
            if (isFirstChar) x -= metrics.leftBearing(ch);
            isFirstChar = false;
            if (slashMode!=SLASHING) p->drawText (x, y, QString(ch));
            x += metrics.width(ch);
        }
    }
}


void TmgTreemapView::drawNode (QPainter* p, TmNode* node, int x, int y, char* txt)
{
    if (node==NULL) return;
    
    char text[1000];
    if (txt==NULL) {
        computeNodeText (text, 1000, node);
        txt=text;
    }

    QRect nodeBox, textBox;
    li.boxes (nodeBox, textBox, node, x, y, txt);

    p->save();
    if (node->isLeaf() && !doShowBIBAsCIB) {
        int w = 2;
        p->setPen (QPen (p->pen().color(), w));
    }
    p->setBrush (getBackgroundBrush(node, QBrush(Qt::white)));
    int percX = 2*200*li.marginX / nodeBox.width();
    int percY = 2*200*li.marginY / nodeBox.height();
    p->drawRoundRect (nodeBox, percX, percY);
    drawTextWithSlashes (p, textBox, txt);
    p->restore();
}


void TmgTreemapView::drawNodeDecorations (TmNode* node, QPainter* p, const QRect& nodeBox,
                                         const QPoint& attachPChild0, const QPoint& attachPChild1, 
                                         const QPoint& attachPParent)
{
    pair<DecMI, DecMI> range = decoration.equal_range (node->index);
    for (DecMI it=range.first; it!=range.second; it++) {
//        printf ("%d: %s\n", node->index, (*it).second->group);
        (*it).second->drawNode (p, li, nodeBox, attachPChild0, attachPChild1, attachPParent);
    }
}


void TmgTreemapView::decorationsExtraWidth (TmNode* node, QPainter* p, int& left, int& right)
{
    left = right = 0;
    pair<DecMI, DecMI> range = decoration.equal_range (node->index);
    for (DecMI it=range.first; it!=range.second; it++) {
        int l, r;
        (*it).second->extraWidth (p, li, l, r);
        if (l>left) left = l;
        if (r>right) right = r;
    }
}


void TmgTreemapView::addFeatureText (char* text, int& ctr, TmExtendedFeatureList& il, int& i)
{
  int n;
  treemap->nameOfFeature (text+ctr, il[i].id, n);
  int j=0;
  while (j<n && il[i+j].id==il[i].id+j) j++;
  i+=j;
  ctr += strlen (text+ctr);
} 


void TmgTreemapView::lineBreak (char* text)
{
  int width = 14;
  int i = 0;
  while (text[i]!='\0') {
    int ctr  = 0;    
    while (text[i]!='\0' && ctr<width) {
      if (text[i]!='[' && text[i]!=']') ctr++;
      if (text[i]=='\n') ctr = 0;      
      i++;
    }    
    while (text[i]!='\n' && text[i]!=' ' && text[i]!='\0') i++;
    if (text[i]==' ') {
      text[i]='\n';
      i++;
    }    
  }
}


void TmgTreemapView::computeNodeText (char* text, int maxLen, TmNode* node)
{
    text[0] = '\0';
    if (node==NULL) {
        return;
    }
    // Print flags
    char rlFlag=' ', iibFlag=' ', opFlag=' ', nlFlag=' ', mFlag=' ', iFlag=' ';
    if (node->isFlag (TmNode::IS_OPTIMIZED)) opFlag = 'O';
    if (node->isFlag (TmNode::IS_FEATURE_PASSED_VALID)) rlFlag = 'F';
    if (node->isFlag (TmNode::IS_GAUSSIAN_VALID)) iibFlag = 'G';
    if (node->linearizationPointFeature>=0) nlFlag='N';
    if (node->isFlag (TmNode::CAN_BE_MOVED)) mFlag = 'M';
    if (node->isFlag (TmNode::CAN_BE_INTEGRATED)) iFlag = 'I';    
    sprintf(text, "%3d:%c%c%c%c %5.3fms\n", node->index, iFlag, opFlag, iibFlag, mFlag, node->worstCaseUpdateCost*1000);
    

    // Now print the landmark list
    int ctr = strlen(text);


    if (node->isFlag (TmNode::IS_FEATURE_PASSED_VALID)) {      
      const TmExtendedFeatureList& rl = node->featurePassed;
      TmExtendedFeatureList il;
      node->computeFeaturesInvolved (il);
      // Run through 'il' and print all landmarks write a landmark with brackets,
      // if it is not contained in 'rl'. ('rl' is a subset of 'il')
      // This shows all landmarks involved in a node and those in brackets that
      // got eliminated there.
      int j=0;
      int n = (int) il.size();
      for (int i=0; i<n;) {  // i++ by addFeatureText
        text[ctr]=' ';        
        ctr++;
        if (ctr>maxLen-10) break;
        while (j<(int) rl.size() && rl[j].id<il[i].id) j++;
        bool useBracket  = 
          (j==(int) rl.size() || rl[j].id!=il[i].id); // Not in 'rl'
        if (useBracket) {
          if (doShowEliminated) {
            text[ctr]='[';
            ctr++;
            addFeatureText (text, ctr, il, i);
            text[ctr]=']';
            ctr++;
          }
        }
        else addFeatureText (text, ctr, il, i);
      }
      text[ctr]='\0';
    }

    checkValidity (text);    
    lineBreak (text);
}


void TmgTreemapView::addDecoration (TmgDecoration* newDecoration)
{
    if (newDecoration->node==NULL) return;
    newDecoration->view = this;
    DecorationMap::value_type val (newDecoration->node->index, newDecoration);
    decoration.insert (decoration.begin(), val);
    updateContent();
}


void TmgTreemapView::deleteDecoration (TmgDecoration* myDecoration)
{
    if (myDecoration->node==NULL) return;
    updateContent();
    pair<DecMI, DecMI> range = decoration.equal_range (myDecoration->node->index);
    for (DecMI it=range.first; it!=range.second; it++) {
        if ((*it).second==myDecoration) {
            decoration.erase (it);
            return;
        }
    }
    assert (false); // Not found
}


void TmgTreemapView::deleteDecorationGroup (const char* groupName)
{
    DecMI it=decoration.begin();
    while (it!=decoration.end()) {
        DecMI it2 = it;
        it++;
        if (strcmp((*it2).second->group, groupName)==0) {
            decoration.erase (it2);
        }
    }
    updateContent();
}

    
void TmgTreemapView::deleteAllDecorations()
{
    for (DecMI it=decoration.begin(); it!=decoration.end(); it++)
        delete (*it).second;
    decoration.clear();
    updateContent();
}


void TmgTreemapView::decorateRegions (TmNode* r, TmNode* n, char* name)
{
  /* NOT IMPLEMENTED in the new version
    if (getTreemap()==NULL) return;
    if (r==NULL || n==NULL || !n->isAncestor(r)) return;

    TmTreemap::Path path;
    
    getTreemap()->getPath (path, r, n);
    addDecoration (new TmgDecorationEdgemark (r, name, TmgDecorationEdgemark::UPWARDBRACKET, Qt::red, "D"));
    if (!r->isLeaf())
        addDecoration (new TmgDecorationEdgemark (r->child[1-path[0].whichChild()), 
                                             name, TmgDecorationEdgeMark::DOWNWARDBRACKET, Qt::red, "C"));
    for (int i=1; i<(int) path.size()-1; i++) {
        if (!path[i].node->isLeaf()) 
            addDecoration (new TmgDecorationEdgeMark (path[i].node->child[1-path[i].whichChild()), 
                                                     name, TmgDecorationEdgeMark::DOWNWARDBRACKET, Qt::red, "B"));
    }
    addDecoration (new TmgDecorationEdgeMark (n, name, TmgDecorationEdgeMark::DOWNWARDBRACKET, Qt::red, "A"));
    if (!n->isLeaf()) {
        addDecoration (new TmgDecorationEdgeMark (n->child[0], name, TmgDecorationEdgeMark::DOWNWARDBRACKET, Qt::red, "A0"));
        addDecoration (new TmgDecorationEdgeMark (n->child[1], name, TmgDecorationEdgeMark::DOWNWARDBRACKET, Qt::red, "A1"));
    }
  */
}


void TmgTreemapView::printOptions (const char* options, QFont& font)
{
    if (drawMode==POSTSCRIPT) {
        font = QFont("Helvetica", 6);
    }
    else font = this->font();
    int idx=0;
    onlyNode = false;
    doShowBIBAsCIB = false;
    doShowEliminated = true;
    while (options[idx]!='\0') {
        if (parse(options, idx, "tree")) {}
        else if (parse(options, idx, "font")) {
            int fontsize, idxincr;
            if (sscanf(options+idx,"%d %n",&fontsize, &idxincr)>=1) {
                idx += idxincr;
                font.setPointSize(fontsize);
            }
            else QMessageBox::information (this, "Print tree", "Illegal font size");
        }
        else if (parse(options, idx, "status")) {
            doShowStatusLine = true;
        }
        else if (parse(options, idx, "nostatus")) {
            doShowStatusLine = false;
        }
        else if (parse(options, idx, "noeliminated")) {
            doShowEliminated = false;
        }
        else if (parse(options, idx, "showeliminated")) {
            doShowEliminated = true;
        }
        else if (parse(options, idx, "onlynode")) {
            onlyNode = true;
        }
        else if (parse(options, idx, "BIBasCIB")) {
            doShowBIBAsCIB = true;
        }
        else {
            QMessageBox::information (this, "Print tree", "Illegal option "+QString(options+idx));
            while (options[idx]!=' ' && options[idx]!='\0') idx++;
            while (options[idx]==' ') idx++;
        }
    }
}


const char* TmgTreemapView::extension(const char* file)
{
    int i = strlen(file);
    while (i>0 && file[i]!='.') i--;
    return file+i;
}


bool TmgTreemapView::parse (const char* options, int& idx, const char* token)
{
    int n = strlen(token);
    if (strncmp(options+idx, token, n)==0 && (options[idx+n]==' ' || options[idx+n]=='\0')) {
        idx += n;
        while (options[idx]==' ') idx++;
        return true;
    }
    else return false;
}


void TmgTreemapView::setBoundingBox (const char* file1, const char* file2, const QRect& box, const QRect& pageBox)
{
    int left, right, top, bottom;
    int margin = 2;
    
    left   = box.left() - pageBox.left()-margin;
    right  = box.right() - pageBox.left()+margin;
    top    = pageBox.bottom() - box.top()+margin;
    bottom = pageBox.bottom() - box.bottom()-margin;

    char buffer[1000];
    FILE* f1 = fopen (file1, "r");
    FILE* f2 = fopen (file2, "w");
    fprintf (f2, "%%!PS-Adobe-3.0\n");
    fprintf (f2, "%%%%Creator: slamgui by Udo Frese\n");
    fprintf (f2, "%%%%BoundingBox: %d %d %d %d\n", left, bottom, right, top);
    while (!feof(f1)) {
        fgets (buffer, 1000, f1);
        if (buffer[0]=='\0') break;
        if (strncmp (buffer, "%!", 2)==0) {
            // skip line
        }
        else if (strncmp (buffer, "%%BoundingBox:", 14)==0) {
            // skip line
        }
        else fputs (buffer, f2);
    }
    fclose (f1);
    fclose (f2);
}


void TmgTreemapView::print (const char* file, const char* options) throw (runtime_error)
{
    if (file==NULL) throw runtime_error("No filename specified");
    const char* ext = extension(file);
    

    QPrinter* printer = NULL;


    bool oldDoShowStatusLine = doShowStatusLine;
    if ((strcmp(ext, ".ps")==0 || strcmp(ext, ".eps")==0)) {
        printer = new QPrinter();
        drawMode = POSTSCRIPT;
    }
    else drawMode = BITMAP;

    if (strcmp(ext, ".ps")==0) {
        printer->setOutputFileName(file);
        if (!printer->setup(this)) throw runtime_error("Could not setup printer");
    }

    QFont fnt;
    printOptions (options, fnt);

    if (drawMode==POSTSCRIPT) {
        // Do .ps and .eps via QPrinter
        if (strcmp(ext, ".eps")==0) {
            printer->setOutputFileName("/tmp/slamgui.ps");
//            printer->setPageSize (QPrinter::A4);
            printer->setFullPage (true);
        }
        QRect pageBox;
        {            
            QPainter p;
            if (!p.begin(printer)) {
                QMessageBox::information (this, "Print", "error in opening printer\n");
                throw runtime_error("error in opening printer\n");
            }
            p.setFont(fnt);
            QRect vp = p.viewport();
            drawAll (&p, vp.left(), vp.top(), vp.width(), vp.height(), true);
            p.end();
            printer->newPage();
            QPaintDeviceMetrics met (printer);
            pageBox = QRect (0, 0, met.width(), met.height());
        }
    
        delete printer;
        if (strcmp(ext,".eps")==0) {
            setBoundingBox ("/tmp/slamgui.ps", file, lastTreeSize, pageBox); 
        }
        printer = NULL;
    }
    else {
        // Do all pixel based format via QPixmap
        QStrList outputFormats = QImage::outputFormats();
        QStrListIterator of(outputFormats);
        if (ext!=NULL && ext[0]=='.') ext++; // skip the '.'
        QString extLower = QString(ext).lower();
        if (extLower=="jpg") extLower="jpeg";
        while (of.current()!=NULL) {
            QString fmt = QString(of.current()).lower();
            if (extLower==fmt) { // format supported
                break;
            }
            ++of;
        }
        if (of.current()!=NULL) {
            // We must render twice, since the first rendering sets 'lastTreeSize'
            QPixmap img(100,100);
            QPainter p (&img);
            p.setFont(fnt);
            QRect vp = p.viewport();
            drawAll (&p, vp.left(), vp.top(), vp.width(), vp.height(), true);
            p.end();

            img.resize (lastTreeSize.right(), lastTreeSize.bottom());
            img.fill (Qt::white);
            p.begin (&img);
            p.setFont(fnt);
            vp = p.viewport();
            drawAll (&p, vp.left(), vp.top(), vp.width(), vp.height(), false);
            p.end();
            img.save (file, of.current());
        }
        else {
            QMessageBox::information (this, "Print", "Format "+extLower+" not supported");
            throw runtime_error ("Format "+extLower+" not supported");
        }
    }

    doShowStatusLine = oldDoShowStatusLine;
    drawMode = SCREEN;
    updateContent();
}


void TmgTreemapView::checkValidity (char* txt)
{
  int i=0;
  while (txt[i]!='\0') i++;
}

void TmgTreemapView::savePng (const char* filename)
{
  try {
    int w = contentsWidth();      
    int h = contentsHeight();      
    
    QPixmap img (w, h);
    QPainter p(&img);
    p.fillRect (0, 0, w, h, backgroundColor());      
    drawContents (&p, 0, 0, contentsWidth(), contentsHeight());
    img.save (filename, "PNG");      
  }
  catch (const exception& error) {
    QMessageBox::critical (this, "TmgTreemapView::savePng", error.what());
  }  
}
