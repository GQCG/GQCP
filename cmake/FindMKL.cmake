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
find_library(MKL_INTERFACE_LIBRARY NAMES libmkl_intel_lp64.a PATHS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel64)
find_library(MKL_INTEL_THREAD NAMES libmkl_intel_thread.a PATHS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel64)
find_library(MKL_CORE_LIBRARY NAMES libmkl_core.a PATHS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel64)

# get_filename_component(IOMP5_ROOT $ENV{MKLROOT}/.. ABSOLUTE) # The root of IOMP installation can be one directory above the MKL root.
# find_library(MKL_IOMP5 NAMES libiomp5.so libiomp5_db.so libiompstubs5.so PATHS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel64 ${IOMP5_ROOT}/lib ${IOMP5_ROOT}/lib/intel64 REQUIRED) # This library should be dynamically linked.
find_library(MKL_GOMP NAMES libgomp.so PATHS $ENV{MKLROOT}/lib $ENV{MKLROOT}/lib/intel64 ${GOMP_ROOT}/lib ${GOMP_ROOT}/lib/intel64 REQUIRED)

mark_as_advanced(MKL_INCLUDE_DIR MKL_INTERFACE_LIBRARY MKL_INTEL_THREAD MKL_CORE_LIBRARY MKL_GOMP)

find_package_handle_standard_args(MKL REQUIRED_VARS
  MKL_INCLUDE_DIR MKL_INTERFACE_LIBRARY MKL_INTEL_THREAD MKL_CORE_LIBRARY MKL_GOMP) # sets MKL_FOUND

if(MKL_FOUND AND NOT TARGET MKL::MKL)
  add_library(MKL::MKL INTERFACE IMPORTED)
  target_include_directories(MKL::MKL INTERFACE ${MKL_INCLUDE_DIR})

  if(NOT APPLE)
    target_link_libraries(MKL::MKL INTERFACE -Wl,--start-group ${MKL_INTERFACE_LIBRARY} ${MKL_INTEL_THREAD} ${MKL_CORE_LIBRARY} ${MKL_GOMP} -Wl,--end-group gomp pthread m dl)
  else()
    target_link_libraries(MKL::MKL INTERFACE ${MKL_INTERFACE_LIBRARY} ${MKL_INTEL_THREAD} ${MKL_CORE_LIBRARY} ${MKL_GOMP} gomp pthread m dl)
  endif()
endif()