# Adds all files from a subfolder for compilation
# can handle ui,moc, and rc files correctly
# Usage: AddQtFolders( OUTPUT_VAR )
#        add_executable (myexec  $OUTPUT_VAR )

macro ( add_qt_folders OUT_FILES  )

    foreach( DIRECTORY ${ARGN} )
    	set( ALL_HEADER )
    	set( MOC_HEADER )	
    	set( CPP )
    	set( UI )
    	set( RES )
    	set( MOC_OUTPUT_FILES )
    	set( UI_OUTPUT_FILES )
    	set( RES_OUTPUT_FILES )
    
    	file( GLOB ALL_HEADER ${CMAKE_CURRENT_SOURCE_DIR}/${DIRECTORY}/*.hh  )
    	file( GLOB CPP        ${CMAKE_CURRENT_SOURCE_DIR}/${DIRECTORY}/*.cc  )
    	file( GLOB UI         ${CMAKE_CURRENT_SOURCE_DIR}/${DIRECTORY}/*.ui  )
    	file( GLOB RES        ${CMAKE_CURRENT_SOURCE_DIR}/${DIRECTORY}/*.qrc )
    
    	# Find files which need moc processing 
    	# (are recognized on Q_OBJECT macro inside file)
    	foreach( _current_HEADER ${ALL_HEADER} )
    		GET_FILENAME_COMPONENT(_abs_HEADER ${_current_HEADER} ABSOLUTE)
    		FILE( READ ${_abs_HEADER} _contents)
    		STRING( REGEX MATCHALL "Q_OBJECT" _match  "${_contents}" )
    
    		IF( _match)	
    			LIST( APPEND MOC_HEADER ${_current_HEADER} )
    		ENDIF (_match)
    	endforeach( _current_HEADER )
    
    	QT4_WRAP_CPP(MOC_OUTPUT_FILES ${MOC_HEADER})
    	QT4_WRAP_UI(UI_OUTPUT_FILES ${UI} )
    	QT4_ADD_RESOURCES(RES_OUTPUT_FILES ${RES})

	    LIST( APPEND ${OUT_FILES} ${CPP} ${UI_OUTPUT_FILES} ${MOC_OUTPUT_FILES} ${RES_OUTPUT_FILES} )
    endforeach( DIRECTORY)
endmacro (add_qt_folders)


