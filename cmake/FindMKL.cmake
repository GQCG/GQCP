
#[=======================================================================[.rst:
FindMKL
--------

Find the Intel MKL libraries https://software.intel.com/en-us/mkl on the system.

Result Variables
^^^^^^^^^^^^^^^^

This module makes a ``MKL::MKL``target and will set the following variables in your project:

``MKL_FOUND``
  System has the MKL libraries installed.
``MKL_INCLUDE_DIR``
  The MKL include directories.
``MKL_LIBRARIES``
  The MKL libraries.

Hints
^^^^^

``$ENV{MKLROOT}``
  Environment variables that defines the root directory of the MKL installation.
#]=======================================================================]
include(FindPackageHandleStandardArgs)

find_path(MKL_INCLUDE_DIR mkl.h HINTS $ENV{MKLROOT}/include)
find_library(MKL_SINGLE_DYNAMIC_LIBRARY NAMES libmkl_rt.so libmkl_rt.dylib PATHS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel64)
get_filename_component(IOMP5_ROOT $ENV{MKLROOT}/.. ABSOLUTE) # The root of IOMP installation can be one directory above the MKL root.
find_library(MKL_IOMP5 NAMES libiomp5.so libiomp5.dylib PATHS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel64 ${IOMP5_ROOT}/lib ${IOMP5_ROOT}/lib/intel64) # This library should be dynamically linked.

mark_as_advanced(MKL_INCLUDE_DIR MKL_SINGLE_DYNAMIC_LIBRARY MKL_IOMP5)

find_package_handle_standard_args(MKL REQUIRED_VARS
    MKL_INCLUDE_DIR MKL_SINGLE_DYNAMIC_LIBRARY MKL_IOMP5) # sets MKL_FOUND

if(MKL_FOUND AND NOT TARGET MKL::MKL)
    add_library(MKL::MKL INTERFACE IMPORTED)
    target_include_directories(MKL::MKL INTERFACE ${MKL_INCLUDE_DIR})
    target_link_libraries(MKL::MKL INTERFACE  ${MKL_SINGLE_DYNAMIC_LIBRARY} ${MKL_IOMP5} pthread m dl)
endif()
