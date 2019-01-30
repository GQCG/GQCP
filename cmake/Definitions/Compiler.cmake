# Set some compiler-related CMake variables


# Compiler optimizations (https://stackoverflow.com/a/41361741/7930415)
if (“${CMAKE_CXX_COMPILER_ID}” STREQUAL “Clang”)
    set(CMAKE_CXX_FLAGS_DEBUG "-g -m64 -pipe")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native -m64 -pipe")
elseif (“${CMAKE_CXX_COMPILER_ID}” STREQUAL “GNU”)
    set(CMAKE_CXX_FLAGS_DEBUG "-g -march=native -m64 -pipe -pthread")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native -m64 -pipe -pthread")
elseif (“${CMAKE_CXX_COMPILER_ID}” STREQUAL “Intel”)
    set(CMAKE_CXX_FLAGS_RELEASE "-g -xHost -m64 -pipe")
    set(CMAKE_CXX_FLAGS_RELEASE "-O2 -xHost -m64 -pipe")
endif()


# From RPATH to full installation path
# https://stackoverflow.com/questions/30398238/cmake-rpath-not-working-could-not-find-shared-object-file
set(CMAKE_INSTALL_RPATH "${LIBRARY_INSTALL_DIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
