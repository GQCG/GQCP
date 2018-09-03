# src-level CMakeLists.txt


# We want to make a static library target called "gqcg", based on the project source files
add_library(${LIBRARY_NAME} ${LIBRARY_TYPE} ${PROJECT_SOURCE_FILES} ${PROJECT_INCLUDE_FILES})


# Configure the library: include headers and link to libraries
include(${CMAKE_SOURCE_DIR}/cmake/ConfigureLibrary.cmake)


# Make the install target for this library
include(${CMAKE_SOURCE_DIR}/cmake/InstallLibrary.cmake)
