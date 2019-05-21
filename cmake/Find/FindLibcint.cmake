# FindLibcint.cmake. Find<package>.cmake-template from (https://cmake.org/Wiki/CMake:How_To_Find_Libraries#Writing_find_modules)
# Try to find Libcint

# The following installations of Libcint are supported:
#   - /usr/local/
#   - any directory, but with an environment variable ${LIBCINTROOT} or ${LIBCINT_ROOT}, which contains the include, lib and share folders

# If found, this will define
#  Libcint_FOUND           Libcint is available on the system
#  Libcint_INCLUDE_DIRS    the Libcint include directories + dependency include directories
#  Libcint_LIBRARIES       the Libcint library + dependency libraries


# When the user has set ${LIBCINTROOT} or ${LIBCINT_ROOT} in the enviroment, this path can also be used

find_path(LIBCINT_PREFIX include/cint.h HINTS /usr/local/*/ $ENV{LIBCINTROOT} $ENV{LIBCINT_ROOT} ${LIBCINTROOT} ${LIBCINT_ROOT})


if("${LIBCINT_PREFIX}" STREQUAL "LIBCINT_PREFIX-NOTFOUND")
    message(FATAL_ERROR "Libcint was not found in the default location /usr/local or through the environment variables")

else()
    # Set FOUND
    set(Libcint_FOUND TRUE)

    # Set the INCLUDE_DIRS
    set(Libcint_INCLUDE_DIRS "${Libcint_INCLUDE_DIRS};${LIBCINT_PREFIX}/include")

    # Set the dynamic or static libraries
    if(EXISTS ${LIBCINT_PREFIX}/lib/libcint.so)
        set(Libcint_LIBRARIES "${LIBCINT_PREFIX}/lib/libcint.so")
    else()
        if(EXISTS ${LIBCINT_PREFIX}/lib/libcint.dylib)
            set(Libcint_LIBRARIES "${LIBCINT_PREFIX}/lib/libcint.dylib")
        else()
            set(Libcint_LIBRARIES "${LIBCINT_PREFIX}/lib/libcint.a")
        endif()
    endif()

    message(STATUS "Libcint was found at ${LIBCINT_PREFIX}")

endif()

if (Libcint_FOUND AND NOT TARGET libcint::libcint)
    add_library(libcint::libcint SHARED IMPORTED)
    set_target_properties(libcint::libcint
            PROPERTIES
                INTERFACE_INCLUDE_DIRECTORIES ${Libcint_INCLUDE_DIRS}
                IMPORTED_LOCATION ${Libcint_LIBRARIES}
            )
endif()
