# CMakeLists for the tagged spectral function libary

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add taggedSF library
add_library(taggedSF taggedSF.cpp)
target_link_libraries(taggedSF ${ROOT_LIBRARIES})
