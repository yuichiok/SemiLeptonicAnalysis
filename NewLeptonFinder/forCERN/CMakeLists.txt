########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( ForCERN )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 2 )
SET( ${PROJECT_NAME}_VERSION_MINOR 0 )
SET( ${PROJECT_NAME}_VERSION_PATCH 1 )

# 1.0.0 for ilcsoft_v01-16 at /afs/desy
# 1.1.0 for ilcsoft_v01-16 at /afs/desy correction of problem with MC particles
# 1.2.0 for ilcsoft_v01-16 at /afs/desy gamma/gamma removal using cut on cos(theta)_PFO < 0.95
# 1.3.0 for ilcsoft_v01-16-02 at /afs/desy gamma/gamma removal using Kt jet algorithm
# 2.0.0 for ilcsoft_v01-16-02 at /afs/desy first test with new LeptonFinder
# 2.0.1 for ilcsoft_v01-16-02 at /afs/desy test with new true B jet id (17/06/13)

### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version
INCLUDE_DIRECTORIES( ${Marlin_INCLUDE_DIRS} )
LINK_LIBRARIES( ${Marlin_LIBRARIES} )
ADD_DEFINITIONS( ${Marlin_DEFINITIONS} )

# optional package

FIND_PACKAGE( ROOT )
IF( ROOT_FOUND )
    INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${ROOT_LIBRARIES} )
    ADD_DEFINITIONS( "-DMARLIN_USE_ROOT" )
    MESSAGE( STATUS "ROOT -- found" )
ELSE()
    MESSAGE( STATUS "ROOT -- not found" )
ENDIF()

FIND_PACKAGE( LCIO )
IF( LCIO_FOUND )
    INCLUDE_DIRECTORIES( ${LCIO_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${LCIO_LIBRARIES} )
    ADD_DEFINITIONS( "-DMARLIN_USE_LCIO" )
    MESSAGE( STATUS "LCIO -- found" )
ELSE()
    MESSAGE( STATUS "LCIO -- not found" )
ENDIF()

FIND_PACKAGE( CLHEP )
IF( CLHEP_FOUND )
    INCLUDE_DIRECTORIES( ${CLHEP_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${CLHEP_LIBRARIES} )
    ADD_DEFINITIONS( "-DMARLIN_USE_CLHEP" )
    MESSAGE( STATUS "CLHEP -- found" )
ELSE()
    MESSAGE( STATUS "CLHEP -- not found" )
ENDIF()

FIND_PACKAGE( MarlinUtil )
IF( MarlinUtil_FOUND )
    INCLUDE_DIRECTORIES( ${MarlinUtil_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${MarlinUtil_LIBRARIES} )
    ADD_DEFINITIONS( "-DMARLIN_USE_MarlinUtil" )
    MESSAGE( STATUS "MarlinUtil -- found" )
ELSE()
    MESSAGE( STATUS "MarlinUtil -- not found" )
ENDIF()

FIND_PACKAGE( GEAR )
IF( GEAR_FOUND )
    INCLUDE_DIRECTORIES( ${GEAR_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${GEAR_LIBRARIES} )
    ADD_DEFINITIONS( "-DMARLIN_USE_GEAR" )
    MESSAGE( STATUS "GEAR -- found" )
ELSE()
    MESSAGE( STATUS "GEAR -- not found" )
ENDIF()


### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



### LIBRARY #################################################################

# definitions to pass to the compiler
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( ./include )
#INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()

