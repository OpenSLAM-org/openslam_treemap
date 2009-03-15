PROJECT (treemap)
#SET(CMAKE_VERBOSE_MAKEFILE ON)

include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}
    ..
    ../lib
)


IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
MESSAGE ("Using no optimization for treemap")
ADD_DEFINITIONS (-DUSEPROCESSTIME -DSEGMENTATION_FAULT_ASSERT -DASSERT_LEVEL=2)
ENDIF (CMAKE_BUILD_TYPE STREQUAL "Debug")

IF (CMAKE_BUILD_TYPE STREQUAL "Release")
MESSAGE ("Using optimization for treemap")
ADD_DEFINITIONS (-DUSEPROCESSTIME -DUSESSE)
ENDIF(CMAKE_BUILD_TYPE STREQUAL "Release")

ADD_SUBDIRECTORY (../lib/xymatrix xymatrix)
ADD_SUBDIRECTORY (../lib/xymlapack xymlapack)
ADD_SUBDIRECTORY (../lib/xycontainer xycontainer)

FILE(GLOB treemap_SRC *.cc)

ADD_LIBRARY(treemap
  ${treemap_SRC}
  )
