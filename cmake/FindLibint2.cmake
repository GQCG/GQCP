# Findlibint2.cmake. Find<package>.cmake-template from (https://cmake.org/Wiki/CMake:How_To_Find_Libraries#Writing_find_modules)
# Try to find libint2

# The following installations of libint2 are supported:
#   - /usr/local/libint
#   - any directory, but with an environment variable ${LIBINTROOT}, which contains the include, lib and share folders

# If found, this will define
#  libint2_FOUND            libint2 is available on the system
#  libint2_INCLUDE_DIRS     the libint2 include directories + dependency include directories
#  libint2_LIBRARIES        the libint2 library + dependency libraries


# We have to find the version! Luckily, libint2 -when installed defaultly- provides a directory /usr/local/libint/x.y.z
# When the user has set ${LIBINTROOT} in the enviroment, this path can also be used
find_path(LIBINT_PREFIX include/libint2.hpp HINTS /usr/local/libint/*/ ENV{LIBINTROOT} ${LIBINT_ROOT})


if("${LIBINT_PREFIX}" STREQUAL "LIBINT_PREFIX-NOTFOUND")
    message(FATAL_ERROR "libint2 was not found in the default location /usr/local/libint/x.y.z or through the environment variables")
else()
    # Set FOUND
    set(libint2_FOUND TRUE)

    # Set the INCLUDE_DIRS
    set(libint2_INCLUDE_DIRS "${libint2_INCLUDE_DIRS};${LIBINT_PREFIX}/include")

    # Set the LIBRARIES
    set(libint2_LIBRARIES "${libint2_LIBRARIES};${LIBINT_PREFIX}/lib/libint2.a")

    message(STATUS "libint2 was found at ${LIBINT_PREFIX}")
endif()
