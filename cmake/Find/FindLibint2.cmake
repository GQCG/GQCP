# FindLibint2.cmake. Find<package>.cmake-template from (https://cmake.org/Wiki/CMake:How_To_Find_Libraries#Writing_find_modules)
# Try to find Libint2

# The following installations of Libint2 are supported:
#   - /usr/local/libint
#   - any directory, but with an environment variable ${LIBINTROOT} or ${LIBINT_ROOT}, which contains the include, lib and share folders

# If found, this will define
#  Libint2_FOUND            Libint2 is available on the system
#  Libint2_INCLUDE_DIRS     the Libint2 include directories + dependency include directories
#  Libint2_LIBRARIES        the Libint2 library + dependency libraries


# We have to find the version! Luckily, Libint2 -when installed defaultly- provides a directory /usr/local/libint/x.y.z
# When the user has set ${LIBINTROOT} or ${LIBINT_ROOT} in the enviroment, this path can also be used
find_path(LIBINT_PREFIX include/libint2.hpp HINTS /usr/local/libint/*/ $ENV{LIBINTROOT} $ENV{LIBINT_ROOT} ${LIBINTROOT} ${LIBINT_ROOT})


if("${LIBINT_PREFIX}" STREQUAL "LIBINT_PREFIX-NOTFOUND")
    message(FATAL_ERROR "Libint2 was not found in the default location /usr/local/libint/x.y.z or through the environment variables")
else()
    # Set FOUND
    set(Libint2_FOUND TRUE)

    # Set the INCLUDE_DIRS
    set(Libint2_INCLUDE_DIRS "${Libint2_INCLUDE_DIRS};${LIBINT_PREFIX}/include")

    # Set the dynamic or static libraries
    if(EXISTS ${LIBINT_PREFIX}/lib/libint2.so)
        set(Libint2_LIBRARIES "${Libint2_LIBRARIES};${LIBINT_PREFIX}/lib/libint2.so")
    else()
        set(Libint2_LIBRARIES "${Libint2_LIBRARIES};${LIBINT_PREFIX}/lib/libint2.a")
    endif()

    message(STATUS "Libint2 was found at ${LIBINT_PREFIX}")
endif()
