include_directories(
    ${CMAKE_CURRENT_BINARY_DIR}
    ..
)

FILE(GLOB vectormath_SRC *.cc)

ADD_LIBRARY(vectormath
  ${vectormath_SRC}
  )
