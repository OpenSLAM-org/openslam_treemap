PROJECT(xymatrix)

include_directories(
  ..
)

add_library(xymatrix
  xymGeneral.cc
  xymMatrixC.cc
  xymMatrixVC.cc
  xymOperations.cc
  xymVector.cc
  xymVectorV.cc
)
