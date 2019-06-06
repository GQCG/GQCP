#[=======================================================================[.rst:
FindSpectra
--------

Find the Spectra header-only library (https://spectralib.org) on the system.

Result Variables
^^^^^^^^^^^^^^^^

This module makes a ``Spectra::Spectra``target and will set the following variables in your project:

``Spectra_FOUND``
  System has the Spectra library installed.
``Spectra_INCLUDE_DIR``
  The Spectra include directories.

Hints
^^^^^

``Spectra_ROOT_DIR``
  Define the root directory of the Spectra installation.
#]=======================================================================]
find_path(Spectra_INCLUDE_DIR Spectra/SymEigsSolver.h HINTS ${CMAKE_SOURCE_DIR}/spectra/include ${Spectra_ROOT_DIR}/include)
mark_as_advanced(Spectra_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Spectra REQUIRED_VARS
        Spectra_INCLUDE_DIR) # sets Spectra_FOUND

if(Spectra_FOUND AND NOT TARGET Spectra::Spectra)
    add_library(Spectra::Spectra INTERFACE IMPORTED)
    target_link_libraries(Spectra::Spectra INTERFACE Eigen3::Eigen)
    set_target_properties(Spectra::Spectra PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES ${Spectra_INCLUDE_DIR}
            INTERFACE_COMPILE_FEATURES cxx_std_11
            )
endif()
