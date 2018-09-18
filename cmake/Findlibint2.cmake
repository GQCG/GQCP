# Findlibint2.cmake. Find<package>.cmake-template from (https://cmake.org/Wiki/CMake:How_To_Find_Libraries#Writing_find_modules)

# For the moment, only a default installation of libint2 is supported, i.e. in /usr/local/libint

# Try to find libint2
# Once done, this will define
#  libint2_FOUND            libint2 is available on the system
#  libint2_INCLUDE_DIRS     the libint2 include directories + dependency include directories
#  libint2_LIBRARIES        the libint2 library + dependency libraries



# We have to find the version! Luckily, libint2 -when installed defaultly- provides a directory /usr/local/libint/x.y.z
# Try to find the libint2 installation directory. For the moment, only a default installation of libint2 is supported, i.e. in /usr/local/libint
find_path(LIBINT_PREFIX include/libint2.hpp HINTS /usr/local/libint/*/)

if("${LIBINT_PREFIX}" STREQUAL "LIBINT_PREFIX-NOTFOUND")
    message(WARNING "libint2 was not found in the default location /usr/local/libint/x.y.z")
else()
    # Set FOUND
    set(libint2_FOUND TRUE)

    # Set the INCLUDE_DIRS
    set(libint2_INCLUDE_DIRS "${libint2_INCLUDE_DIRS};${LIBINT_PREFIX}/include;${LIBINT_PREFIX}/include/libint2")

    # Set the LIBRARIES
    set(libint2_LIBRARIES "${libint2_LIBRARIES};${LIBINT_PREFIX}/lib/libint2.a")
endif()
