//======================================================================================================================
/*!
  Displays a StaggeredGrid in a Qt Window

  Usage:
  \code
         #include "visualization/GridView.hh"
         #include <QApplication>
         //....

         int main( int argc, char** argv )
         {
            //...

            QApplication app(argc, argv);
            GridView gridView;
            gridView.showMaximized();

            gridView.displayGrid( &yourStaggeredGrid );
            app.exec();
         }
   \endcode

 */
//======================================================================================================================


#pragma once


// Set this to 1 to switch to OpenGL rendering of QGraphicsView
#define USE_OPENGL_RENDERING 0




class QGraphicsScene;

#include <QGraphicsView>
#include <QGraphicsEllipseItem>
#include <QGraphicsRectItem>
#include <QGraphicsSceneMouseEvent>

#include <QVector>
#include <QPointF>
#include <QImage>

#include "StaggeredGrid.hh"


class GridItem: public QObject, public QGraphicsItem
{
    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)

    public:

        GridItem( const StaggeredGrid * g, QGraphicsItem * par =0 );

        virtual QRectF boundingRect() const { return bRect; }
        virtual void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);

        void setDrawU ( bool val )  { drawU = val; if(val) drawF = false; update(); saveSettings();}
        void setDrawV ( bool val )  { drawV = val; if(val) drawG = false; update(); saveSettings();}
        void setDrawF ( bool val )  { drawF = val; if(val) drawU = false; update(); saveSettings();}
        void setDrawG ( bool val )  { drawG = val; if(val) drawV = false; update(); saveSettings();}
        void setDrawP ( bool val )  { drawP = val;                        update(); saveSettings();}
        void setDrawCells      ( bool val ) { drawCells   = val; update(); saveSettings(); }
        void setDrawArrows     ( bool val ) { drawArrows  = val; update(); saveSettings(); }
        void setDrawLabels     ( bool val ) { drawLabels  = val; update(); saveSettings(); }
        void setDrawCellIndices( bool val ) { drawIndices = val; update(); saveSettings(); }


        bool getDrawU() const  { return drawU; }
        bool getDrawV() const  { return drawV; }
        bool getDrawF() const  { return drawF; }
        bool getDrawG() const  { return drawG; }
        bool getDrawP () const { return drawP; }
        bool getDrawCells()  const      { return drawCells;   }
        bool getDrawArrows() const      { return drawArrows;  }
        bool getDrawLabels() const      { return drawLabels;  }
        bool getDrawCellIndices()const  { return drawIndices; }

        void setGrid(const StaggeredGrid * g);
        const StaggeredGrid * getGrid() const { return grid; }

    protected:
        /// All drawing options are stored in QSettings
        /// when a new GridItem is created these settings are chosen as initial settings
        void saveSettings();
        void readSettings();

        void drawArray( QPainter * p, const Array & a,
                        QRectF ellipseRect, const QColor & c );

        void drawVelocityArrows( QPainter * p );

        /// Maps coordinates of a staggered grid ( f.e. particle coordinates)
        /// to coordinates which can be passed to the painter
        QPointF gridToScreenCoord( double x, double y );

        const StaggeredGrid * grid;

        /// Distance between cell borders
        qreal gridSpacing;
        /// The diameter of the circle, which represents the pressure
        qreal circleWidth;
        /// The velocities are drawn as ellipses, the ellipses have the
        /// size (circleWidth, ellipseRatio*circleWidth)
        qreal ellipseRatio;

        qreal arrowSize;

        /// Initializes ellipseHoriz,ellipseVert,circle correctly
        /// depends on gridSpacing,circleWidth and ellipseRatio
        void precalculateCirclePositions();

        void precalculateBoundingBox();

        QRectF ellipseHoriz;
        QRectF ellipseVert;
        QRectF circle;

        QRectF bRect;

        bool drawU;
        bool drawV;

        bool drawF;
        bool drawG;

        bool drawP;

        bool drawArrows;
        bool drawLabels;
        bool drawCells;
        bool drawIndices;
};





class GridView : public QGraphicsView
{
    Q_OBJECT

    public:
        GridView( QWidget * parent = 0 );
        virtual ~GridView() {}

        /// Displays given StaggeredGrid
        void displayGrid( const StaggeredGrid * g );


    protected slots:
        void showContextMenu( const QPoint & p );

    protected:
        // Zoom
        virtual void wheelEvent( QWheelEvent *event );
        virtual void scaleView ( qreal scaleFactor  );

        virtual void resizeEvent( QResizeEvent * ev );

        GridItem * gi;
        QGraphicsScene * scene;
};

