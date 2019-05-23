
#[=======================================================================[.rst:
FindInt2
--------

Find the Int2 library (https://github.com/evaleev/libint) on the system.

Result Variables
^^^^^^^^^^^^^^^^

This module will set the following variables in your project:

``Int2_FOUND``
  System has the Int2 library installed.
``Int2_INCLUDE_DIR``
  The Int2 include directories.
``Int2_LIBRARY``
  The Int2 library.

Hints
^^^^^

``Int2_ROOT_DIR``
  Define the root directory of the Int2 installation.
#]=======================================================================]

find_path(Int2_INCLUDE_DIR libint2.hpp HINTS ${Int2_ROOT_DIR}/include)
find_library(Int2_LIBRARY int2 HINTS ${Int2_ROOT_DIR}/lib)
mark_as_advanced(Int2_INCLUDE_DIR Int2_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Int2 REQUIRED_VARS
        Int2_INCLUDE_DIR Int2_LIBRARY)

if(Int2_FOUND AND NOT TARGET Int2::Int2)
    add_library(Int2::Int2 UNKNOWN IMPORTED)
    set_target_properties(Int2::Int2 PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES ${Int2_INCLUDE_DIR}
            IMPORTED_LOCATION ${Int2_LIBRARY}
            )
endif()