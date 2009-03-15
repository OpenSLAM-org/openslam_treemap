PROJECT (treemapgui)
SET(CMAKE_VERBOSE_MAKEFILE ON)

INCLUDE(FindQt3)
SET (QT_MOC_EXECUTABLE ${QT_INCLUDE_DIR}/../bin/moc)
# workaround, because FindQt3 detected the qt4 moc!

include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}
    ..
    ../lib
    ${QT_INCLUDE_DIR}
)



FILE(GLOB treemapgui_SRC *.cc)
QT_WRAP_CPP(treemapgui treemapgui_MOC_SRC tmgTreemapView.h tmgMapView.h)

ADD_SUBDIRECTORY (../lib/xymatrix xymatrix)
ADD_SUBDIRECTORY (../lib/xycontainer xycontainer)
ADD_SUBDIRECTORY (../treemap treemap)


ADD_LIBRARY(treemapgui
  ${treemapgui_SRC}
  ${treemapgui_MOC_SRC}
  )
