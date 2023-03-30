#[=======================================================================[.rst:
FindLibcint
--------

Find the libcint library (https://github.com/sunqm/libcint) on the system.

Result Variables
^^^^^^^^^^^^^^^^

This module makes a ``Libcint::Libcint``target and will set the following variables in your project:

``Libcint_FOUND``
  System has the Libcint library installed.
``Libcint_INCLUDE_DIR``
  The Libcint include directories.
``Libcint_LIBRARY``
  The Libcint library.

Hints
^^^^^

``Libcint_ROOT_DIR``
  Define the root directory of the Libcint installation.
#]=======================================================================]
find_path(Libcint_INCLUDE_DIR cint.h HINTS ${Libcint_ROOT_DIR}/include)
find_library(Libcint_LIBRARY cint HINTS ${Libcint_ROOT_DIR}/lib)
mark_as_advanced(Libcint_INCLUDE_DIR Libcint_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Libcint REQUIRED_VARS
    Libcint_INCLUDE_DIR Libcint_LIBRARY) # sets Libcint_FOUND

if(Libcint_FOUND AND NOT TARGET Libcint::Libcint)
    add_library(Libcint::Libcint UNKNOWN IMPORTED)
    set_target_properties(Libcint::Libcint PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES ${Libcint_INCLUDE_DIR}
        IMPORTED_LOCATION ${Libcint_LIBRARY}
        INTERFACE_COMPILE_FEATURES cxx_std_11
    )
endif()
