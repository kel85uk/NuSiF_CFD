
#include "GridView.hh"

#include <QGraphicsScene>
#include <QGraphicsSceneHoverEvent>
#include <QStyleOptionGraphicsItem>
#include <QDragEnterEvent>
#include <QWheelEvent>
#include <QMimeData>
#include <QRectF>
#include <QFontMetrics>
#include <QMenu>
#include <QAction>
#include <QDebug>
#include <QActionGroup>
#include <QSettings>
#include <QToolBar>
#include <QMessageBox>

#include <cmath>

static const double Pi = 3.14159265358979323846264338327950288419717;
static double TwoPi    = 2.0 * Pi;

// ----------------------------------- Grid Item ---------------------------------------

GridItem::GridItem ( const StaggeredGrid * g, QGraphicsItem * par )
    : QGraphicsItem( par ),
      grid( g ),
      gridSpacing( 200.0 ),
      circleWidth( 25.0 ),
      ellipseRatio( 0.4 ),
      arrowSize( 20 ),
      drawU(true), drawV(true),
      drawF(false),drawG(false),
      drawP(true),
      drawArrows(false),drawLabels(true),
      drawCells(true), drawIndices(false)
{
    precalculateBoundingBox();
    precalculateCirclePositions();

    readSettings();
}

void GridItem::precalculateBoundingBox()
{
    if(grid)
    {
        bRect.setWidth((grid->xSize()+2) * gridSpacing );
        bRect.setHeight((grid->ySize()+2) * gridSpacing);
        bRect.moveTo(0,0);
    }
    else
        bRect = QRectF();
}

void GridItem::precalculateCirclePositions()
{
    circle = QRectF (-circleWidth/2,-circleWidth/2,circleWidth,circleWidth);
    qreal smallerDim = ellipseRatio * circleWidth;
    ellipseHoriz = QRectF(-circleWidth/2,-smallerDim/2,circleWidth,smallerDim);;
    ellipseVert= QRectF (-smallerDim/2,-circleWidth/2,smallerDim,circleWidth);;

    //Move the rectangles to correct start position
    circle.translate(gridSpacing/2,gridSpacing/2);
    ellipseHoriz.translate(gridSpacing,gridSpacing/2);
    ellipseVert.translate(gridSpacing/2,gridSpacing);
}

void GridItem::setGrid(const StaggeredGrid * g)
{
    prepareGeometryChange();
    grid = g;
    precalculateBoundingBox();
    update();
}

void GridItem::saveSettings()
{
    QSettings settings;
    settings.beginGroup("GridItem");

    settings.setValue("drawU",drawU);
    settings.setValue("drawV",drawV);
    settings.setValue("drawF",drawF);
    settings.setValue("drawG",drawG);
    settings.setValue("drawP",drawP);

    settings.setValue("drawArrows",drawArrows);
    settings.setValue("drawLabels",drawLabels);
    settings.setValue("drawCells",drawCells);
    settings.setValue("drawIndices",drawIndices);
    settings.endGroup();
}

void GridItem::readSettings()
{
    QSettings settings;
    settings.beginGroup("GridItem");

    drawU=settings.value("drawU",true).toBool();
    drawV=settings.value("drawV",true).toBool();
    drawF=settings.value("drawF",false).toBool();
    drawG=settings.value("drawG",false).toBool();
    drawP=settings.value("drawP",true).toBool();

    drawArrows= settings.value("drawArrows",true).toBool();
    drawLabels= settings.value("drawLabels",false).toBool();
    drawCells = settings.value("drawCells",true).toBool();
    drawIndices=settings.value("drawIndices",false).toBool();
    settings.endGroup();
}

void GridItem::paint(QPainter *p, const QStyleOptionGraphicsItem *option, QWidget *widget)
{
    if(!grid)
        return;

    p->setClipRect(option->exposedRect);

    //Draw Cells
    if(drawCells)
    {
        p->setPen(Qt::black);
        p->setBrush(Qt::black);

        for( int yi =0; yi<= grid->ySize()+2; ++yi)
            p->drawLine(0, int( yi*gridSpacing) , int( bRect.width() ), int(yi*gridSpacing) );

        for( int xi = 0; xi<=grid->xSize()+2; ++xi)
            p->drawLine( int(xi*gridSpacing),0, int( xi*gridSpacing ), int( bRect.height() ));
    }

    if(drawU ) drawArray(p,grid->u(),ellipseHoriz,Qt::blue);
    if(drawV ) drawArray(p,grid->v(),ellipseVert, Qt::red);

    if(drawF ) drawArray(p,grid->f(),ellipseHoriz,Qt::green);
    if(drawG ) drawArray(p,grid->g(),ellipseVert, Qt::cyan);

    if(drawP ) drawArray(p,grid->p(),circle, Qt::black);

    if(drawArrows ) drawVelocityArrows(p);

}


QPointF GridItem::gridToScreenCoord( double xc, double yc )
{
    xc /= grid->dx();
    yc /= grid->dy();

    // Mirror - y coordinate goes upward in grid, and downward on canvas
    yc = grid->ySize() - yc;

    xc *= gridSpacing;
    yc *= gridSpacing;

    xc += gridSpacing;
    yc += gridSpacing;

    return QPointF(xc,yc);
}


void GridItem::drawVelocityArrows( QPainter * p )
{
   real maxVel = 0;
   for( int yi = 1; yi+1 < grid->v().ySize(); ++yi )
       for( int xi = 1; xi+1 < grid->v().xSize(); ++xi )
          maxVel = std::max( maxVel, std::fabs(grid->v()(xi,yi) ));

   for( int yi = 1; yi+1 < grid->u().ySize(); ++yi )
       for( int xi = 1; xi+1 < grid->u().xSize(); ++xi )
          maxVel = std::max( maxVel, std::fabs(grid->u()(xi,yi) ));


   real velNormalization = 2.0 * gridSpacing / maxVel;

   QPointF startPoint ( 1.5 * gridSpacing, 1.5 * gridSpacing );

   for ( int j = 0; j < grid->ySize (); ++j )
      for ( int i = 0; i < grid->xSize (); ++i )
      {
         qreal norm_u = 0.5 * ( grid->u() ( i, j )  + grid->u() ( i+1, j ) ) * velNormalization;
         qreal norm_v = 0.5 * ( grid->v() ( i, j )  + grid->v() ( i, j+1 ) ) * velNormalization;

         unsigned int inv_j = ( grid->ySize () - 1 ) - j;
         QPointF cellMidPoint = startPoint + QPointF ( i * gridSpacing, inv_j * gridSpacing );
         QPointF sourcePoint = cellMidPoint + QPointF ( -0.5 * norm_u, 0.5 * norm_v );
         QPointF destPoint = cellMidPoint + QPointF ( 0.5 * norm_u, -0.5 * norm_v );

         QLineF line ( sourcePoint, destPoint );
         if (qFuzzyCompare ( line.length (), qreal ( 0. ) ) ) return;

         // Draw the line itself
         p->setPen ( QPen ( Qt::black, 1, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin ) );
         p->drawLine ( line );

         // Draw the arrows
         double angle = ::acos ( line.dx () / line.length () );
         if (line.dy () >= 0 ) angle = TwoPi - angle;

         QPointF destArrowP1 = destPoint
                  + QPointF ( sin ( angle - Pi / 3 ) * arrowSize, cos ( angle - Pi / 3 ) * arrowSize );
         QPointF destArrowP2 = destPoint
                  + QPointF ( sin ( angle - Pi + Pi / 3 ) * arrowSize, cos ( angle - Pi + Pi / 3 ) * arrowSize );

         p->setBrush ( Qt::black );
         p->drawPolygon ( QPolygonF () << line.p2 () << destArrowP1 << destArrowP2 );

      }

}

void GridItem::drawArray(QPainter * p,const Array & a, QRectF ellipseRect, const QColor & c)
{
    //Assume ellipseRect is already at correct start position
    //ellipseRect.translate(offset);

    p->setBrush(QBrush(c));
    p->setPen(Qt::black);

    QRectF ellipseRectOriginal = ellipseRect;


    for( int j = 0; j < a.ySize(); ++j)
    {
        for( int i = 0; i < a.xSize(); ++i)
        {
            p->drawRoundedRect(ellipseRect,40,40,Qt::RelativeSize);
            //p->drawEllipse(ellipseRect);
            ellipseRect.translate(gridSpacing,0);
        }
        ellipseRect.moveLeft(ellipseRectOriginal.left());
        ellipseRect.translate(0,gridSpacing);
    }

    if(!drawLabels)
        return;


    ellipseRect = ellipseRectOriginal;

    p->save();
    p->setBrush(QColor(230,230,230));
    QFontMetrics fontMetrics = p->fontMetrics();

    for( int j = a.ySize()-1; j != -1; --j)
    {
        for( int i = 0; i < a.xSize(); ++i)
        {
            QString s;
            if(drawIndices)
                s = QString("(%1,%2)\n%3").arg(i).arg(j).arg(a(i,j));
            else
                s = QString("%1").arg(a(i,j));

            p->setPen(Qt::NoPen);

            //QRectF textRect (0,0,fontMetrics.width(s)*1.1,fontMetrics.height()*1.1);
            QRectF textRect(QPointF(0,0), fontMetrics.size(0,s));
            textRect.moveCenter(ellipseRect.center() + QPointF(0,ellipseRect.height() + textRect.height()/2));

            p->drawRoundedRect(textRect,40,40,Qt::RelativeSize);
            p->setPen(Qt::black);
            p->drawText(textRect,Qt::AlignCenter,s);

            ellipseRect.translate(gridSpacing,0);
        }
        ellipseRect.moveLeft(ellipseRectOriginal.left());
        ellipseRect.translate(0,gridSpacing);
    }

    p->restore();
}




// ----------------------------------- Grid View ---------------------------------------

#if USE_OPENGL_RENDERING
#include <QGLWidget>
#endif

GridView::GridView(QWidget * par)
    : QGraphicsView (par),
      gi(NULL)
{
    setRenderHints(QPainter::Antialiasing);

#if USE_OPENGL_RENDERING
    setViewport( new QGLWidget( QGLFormat(QGL::SampleBuffers) ) );
#endif

    scene = new QGraphicsScene(this);
    setScene(scene);

    setAcceptDrops(true);
    viewport()->setAcceptDrops(true);
    setDragMode(QGraphicsView::ScrollHandDrag);

    gi = new GridItem(NULL,0);

    scene->addItem(gi);

    //Build up context menu
    setContextMenuPolicy(Qt::CustomContextMenu);
    connect(this, SIGNAL(customContextMenuRequested(const QPoint&)),
            this, SLOT(showContextMenu(const QPoint&)));
}


void GridView::showContextMenu(const QPoint & p)
{
    if(!gi->getGrid())
        return;

    QMenu * contextMenu = new QMenu(this);
    QMenu * horizDisplayMenu = new QMenu("Horizontal Values",contextMenu);
    QMenu * verticalDisplayMenu = new QMenu("Vertical Values",contextMenu);

    // --------- Horizontal Menu ------
    QAction * showU = horizDisplayMenu->addAction("Show U");
    QAction * showF = horizDisplayMenu->addAction("Show F");
    QAction * hideHorizontal = horizDisplayMenu->addAction("Hide");

    QActionGroup * horizGroup = new QActionGroup(this);
    horizGroup->addAction(showU);
    horizGroup->addAction(showF);
    horizGroup->addAction(hideHorizontal);

    showU->setCheckable(true);
    showF->setCheckable(true);
    hideHorizontal->setCheckable(true);
    showU->setChecked(gi->getDrawU());
    showF->setChecked(gi->getDrawF());
    hideHorizontal->setChecked(!gi->getDrawU() && !gi->getDrawV());

    // --------- Vertical Menu ------
    QAction * showV = verticalDisplayMenu->addAction("Show V");
    QAction * showG = verticalDisplayMenu->addAction("Show G");
    QAction * hideVertical = verticalDisplayMenu->addAction("Hide");

    QActionGroup * verticalGroup = new QActionGroup(this);
    verticalGroup->addAction(showV);
    verticalGroup->addAction(showG);
    verticalGroup->addAction(hideVertical);

    showV->setCheckable(true);
    showG->setCheckable(true);
    hideVertical->setCheckable(true);
    showV->setChecked(gi->getDrawV());
    showG->setChecked(gi->getDrawG());
    hideVertical->setChecked(!gi->getDrawV() && !gi->getDrawG());

    // --------- Main Context Menu ------

    contextMenu->addMenu(horizDisplayMenu);
    contextMenu->addMenu(verticalDisplayMenu);
    QAction * showP = contextMenu->addAction("Show P");
    showP->setCheckable(true);
    showP->setChecked(gi->getDrawP());
    contextMenu->addSeparator();

    QAction * showLabels = contextMenu->addAction("Show Numeric Values");
    showLabels->setCheckable(true);
    showLabels->setChecked(gi->getDrawLabels());

    QAction * showCells  = contextMenu->addAction("Show Cells");
    showCells->setCheckable(true);
    showCells->setChecked(gi->getDrawCells());

    QAction * showArrows = contextMenu->addAction("Show Arrows");
    showArrows->setCheckable(true);
    showArrows->setChecked(gi->getDrawArrows());

    QAction * showCellIndices = contextMenu->addAction("Show Cell Indices");
    showCellIndices->setCheckable(true);
    showCellIndices->setChecked(gi->getDrawCellIndices());

    //Show Menu
    QAction * a = contextMenu->exec( mapToGlobal(p));

    if       (a == showU)            gi->setDrawU(true);
    else if (a == showF)            gi->setDrawF(true);
    else if (a == hideHorizontal) { gi->setDrawU(false); gi->setDrawF(false); }

    else if (a == showV)            gi->setDrawV(true);
    else if (a == showG)            gi->setDrawG(true);
    else if (a == hideVertical)   { gi->setDrawV(false); gi->setDrawG(false); }

    else if (a == showP)           gi->setDrawP(a->isChecked());
    else if (a == showLabels)      gi->setDrawLabels(a->isChecked());
    else if (a == showCells )      gi->setDrawCells(a->isChecked());
    else if (a == showArrows)      gi->setDrawArrows(a->isChecked());
    else if (a == showCellIndices) gi->setDrawCellIndices(a->isChecked());

    delete contextMenu;
}



void GridView::displayGrid(const StaggeredGrid * g )
{
    const StaggeredGrid * oldGrid = gi->getGrid();
    gi->setGrid( g );
    if( oldGrid != g )
        fitInView( gi,Qt::KeepAspectRatio );
}


void GridView::resizeEvent(QResizeEvent *)
{
    // Make scrollbars always on, because fitInView is called in resizeEvent
    // set Qt Documentation for fitInView()
    setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOn);
    setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOn);

    if(gi)
        fitInView(gi,Qt::KeepAspectRatio);
}


void GridView::scaleView(qreal scaleFactor)
{
    qreal factor = matrix().scale(scaleFactor, scaleFactor).mapRect(QRectF(0, 0, 1, 1)).width();
    if ( factor < 0.007 || factor > 1000 )
        return;

    scale(scaleFactor, scaleFactor);
}

void GridView::wheelEvent(QWheelEvent *ev)
{
    scaleView(  std::pow((double)2, ev->delta() / 240.0) );
}


