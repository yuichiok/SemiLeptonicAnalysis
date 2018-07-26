###############################################
# cmake configuration file for ForCERN
# @author Jan Engels, DESY
###############################################

SET( ForCERN_FOUND FALSE )
MARK_AS_ADVANCED( ForCERN_FOUND )

# do not store find results in cache
SET( ForCERN_INCLUDE_DIR ForCERN_INCLUDE_DIR-NOTFOUND )

FIND_PATH( ForCERN_INCLUDE_DIR
	NAMES ForCERN.h
	PATHS /exp/flc/rouene/AnalyseTop/DBD/processors/forCERN
	PATH_SUFFIXES include
	NO_DEFAULT_PATH
)
IF( NOT ForCERN_INCLUDE_DIR )
    MESSAGE( STATUS "Check for ForCERN: ${ForCERN_HOME}"
					" -- failed to find ForCERN include directory!!" )
ELSE( NOT ForCERN_INCLUDE_DIR )
    MARK_AS_ADVANCED( ForCERN_INCLUDE_DIR )
ENDIF( NOT ForCERN_INCLUDE_DIR )


# do not store find results in cache
SET( ForCERN_LIB ForCERN_LIB-NOTFOUND )

FIND_LIBRARY( ForCERN_LIB
	NAMES ForCERN
	PATHS /exp/flc/rouene/AnalyseTop/DBD/processors/forCERN
	PATH_SUFFIXES lib
	NO_DEFAULT_PATH
)
IF( NOT ForCERN_LIB )
    MESSAGE( STATUS "Check for ForCERN: ${ForCERN_HOME}"
					" -- failed to find ForCERN library!!" )
ELSE( NOT ForCERN_LIB )
    MARK_AS_ADVANCED( ForCERN_LIB )
ENDIF( NOT ForCERN_LIB )


# set variables and display results
IF( ForCERN_INCLUDE_DIR AND ForCERN_LIB )
    SET( ForCERN_FOUND TRUE )
    SET( ForCERN_INCLUDE_DIRS ${ForCERN_INCLUDE_DIR} )
    SET( ForCERN_LIBRARY_DIRS "/exp/flc/rouene/AnalyseTop/DBD/processors/forCERN/lib" )
	SET( ForCERN_LIBRARIES ${ForCERN_LIB} )
    MARK_AS_ADVANCED( ForCERN_INCLUDE_DIRS ForCERN_LIBRARY_DIRS ForCERN_LIBRARIES )
	MESSAGE( STATUS "Check for ForCERN: ${ForCERN_HOME} -- works" )
ELSE( ForCERN_INCLUDE_DIR AND ForCERN_LIB )
	IF( ForCERN_FIND_REQUIRED )
		MESSAGE( FATAL_ERROR "Check for ForCERN: ${ForCERN_HOME} -- failed!!" )
    ELSE( ForCERN_FIND_REQUIRED )
        MESSAGE( STATUS "Check for ForCERN: ${ForCERN_HOME}"
						" -- failed!! will skip this package..." )
    ENDIF( ForCERN_FIND_REQUIRED )
ENDIF( ForCERN_INCLUDE_DIR AND ForCERN_LIB )
