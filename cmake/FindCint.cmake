
#[=======================================================================[.rst:
FindCint
--------

Find the Cint library (https://github.com/sunqm/libcint) on the system.

Result Variables
^^^^^^^^^^^^^^^^

This module makes a ``Cint::Cint`` target and will set the following variables in your project:

``Cint_FOUND``
  System has the Cint library installed.
``Cint_INCLUDE_DIR``
  The Cint include directories.
``Cint_LIBRARY``
  The Cint library.

Hints
^^^^^

``Cint_ROOT_DIR``
  Define the root directory of the Cint installation.
#]=======================================================================]

find_path(Cint_INCLUDE_DIR cint.h HINTS ${Cint_ROOT_DIR}/include)
find_library(Cint_LIBRARY cint HINTS ${Cint_ROOT_DIR}/lib)
mark_as_advanced(Cint_INCLUDE_DIR Cint_LIBRARY)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cint REQUIRED_VARS
        Cint_INCLUDE_DIR Cint_LIBRARY) # sets Cint_FOUND

if(Cint_FOUND AND NOT TARGET Cint::Cint)
    add_library(Cint::Cint UNKNOWN IMPORTED)
    set_target_properties(Cint::Cint PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES ${Cint_INCLUDE_DIR}
            IMPORTED_LOCATION ${Cint_LIBRARY}
            )
endif()
