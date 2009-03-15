PROJECT (treemap1Mtest)
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

add_definitions(${QT_DEFINITIONS})

#comment this in if the linker complains about MAIN__
#ADD_DEFINITIONS (-DDECLARE_MAIN__)

FILE(GLOB treemap1Mtest_SRC *.cc)

QT_WRAP_CPP(treemap1Mtest treemap1Mtest_MOC_SRC mainWindow.h)


ADD_SUBDIRECTORY (../lib/xymatrix xymatrix)
ADD_SUBDIRECTORY (../lib/xymlapack xymlapack)
ADD_SUBDIRECTORY (../lib/xycontainer xycontainer)
ADD_SUBDIRECTORY (../lib/vectormath vectormath)
ADD_SUBDIRECTORY (../treemap treemap)
ADD_SUBDIRECTORY (../treemapgui treemapgui)

link_directories(
    ${QT_LIBRARY_DIR}
    ${QT_QTCORE_LIBRARY}
    ${QT_QTGUI_LIBRARY_DIR}
)

ADD_EXECUTABLE(treemap1Mtest ${treemap1Mtest_SRC} ${treemap1Mtest_MOC_SRC} )


TARGET_LINK_LIBRARIES(treemap1Mtest
  treemapgui
  treemap
  xymlapack
  xymatrix
  vectormath
  lapack
  blas
  gfortran
  ${QT_LIBRARIES}
)