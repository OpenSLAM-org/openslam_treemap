PROJECT (slamsimulator)
SET(CMAKE_VERBOSE_MAKEFILE ON)

include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}
    ..
    ../lib
)

FILE(GLOB slamsimulator_SRC *.cc)

ADD_SUBDIRECTORY (../lib/vectormath vectormath)

ADD_LIBRARY(slamsimulator
  ${slamsimulator_SRC}
  )
