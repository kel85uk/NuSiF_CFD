
include_directories( . )

#add_executable( nusif FileReader.cc Debug.cc Array.cc main.cc )

add_executable( nusif FileReader.cc StaggeredGrid.cc Geometry2D.cc lodepng.cc GrayScaleImage.cc SORSolver.cc CGSolver.cc Debug.cc VTKWriter.cc EnhancedFluidSimulator.cc EnhancedFluidSimTest.cc)
