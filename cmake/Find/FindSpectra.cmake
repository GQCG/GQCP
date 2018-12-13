# FindSpectra.cmake. Find<package>.cmake-template from (https://cmake.org/Wiki/CMake:How_To_Find_Libraries#Writing_find_modules).

# Try to find Spectra
# Once done, this will define
#  Spectra_FOUND            spectra is available on the system
#  Spectra_INCLUDE_DIRS     the spectra include directories

# Since it's a header-only library, the include directories are enough.

find_path(SPECTRA_PREFIX README.md HINTS ${CMAKE_SOURCE_DIR}/spectra /usr/local/spectra)

if("${SPECTRA_PREFIX}" STREQUAL "SPECTRA_PREFIX-NOTFOUND")
message(WARNING "Spectra was not found.")
else()
# If found, we let the user know that Spectra was found
message(STATUS "Spectra was found at ${SPECTRA_PREFIX}")

# Set Spectra_FOUND
set(Spectra_FOUND TRUE)

# Set the INCLUDE_DIRS
set(Spectra_INCLUDE_DIRS "${SPECTRA_PREFIX}/include")
endif()
