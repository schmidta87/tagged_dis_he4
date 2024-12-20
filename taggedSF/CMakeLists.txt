# CMakeLists for the tagged spectral function libary

include_directories(${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIR})

# Add taggedSF library
add_library(taggedSF taggedSF.cpp)
target_link_libraries(taggedSF ${ROOT_LIBRARIES} GridInterp LineInterp)


#add test executable
add_executable(debug_taggedSF debug_taggedSF.cpp)
target_link_libraries(debug_taggedSF ${ROOT_LIBRARIES} GridInterp LineInterp taggedSF)
