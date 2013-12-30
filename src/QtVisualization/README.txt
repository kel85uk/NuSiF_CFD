QtVisualization
---------------


What is it?
	This folder contains code for a Qt Window that can display a StaggeredGrid.
	The visualization is OPTIONAL, you do not have to use it, however it may help you for debugging your code.
	

How to build?
	Copy the "QtVisualization" folder to your "src" directory  and adapt "src/CMakeLists.txt":
	
		Example:
			find_package( Qt4 REQUIRED )
			set( QT_USE_QTOPENGL TRUE )
			include( ${QT_USE_FILE} )
			include( QtVisualization/AddQtFolders.cmake )
			
			add_qt_folders ( VISUALIZATION_FILES "QtVisualization" )
			
			include_directories( . "QtVisualization" )
			set ( OWN_SOURCES  FileReader.cc Debug.cc Array.cc main.cc SORSolver.cc  ) # and probably many more...
			
			add_executable( nusif ${OWN_SOURCES}  ${VISUALIZATION_FILES}  )
			target_link_libraries( nusif ${QT_LIBRARIES}  )
			

How to use?
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
	
