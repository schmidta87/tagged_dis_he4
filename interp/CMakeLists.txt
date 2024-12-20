# CMakeLists for the interpolation libary

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add taggedSF library
add_library(GridInterp GridInterp.cpp)
target_link_libraries(GridInterp ${ROOT_LIBRARIES})

add_library(LineInterp LineInterp.cpp)
target_link_libraries(LineInterp ${ROOT_LIBRARIES})


