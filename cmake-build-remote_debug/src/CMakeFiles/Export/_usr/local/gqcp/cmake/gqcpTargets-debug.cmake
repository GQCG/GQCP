#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "gqcp" for configuration "Debug"
set_property(TARGET gqcp APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(gqcp PROPERTIES
  IMPORTED_LOCATION_DEBUG "/usr/local/gqcp/lib/libgqcp.so"
  IMPORTED_SONAME_DEBUG "libgqcp.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS gqcp )
list(APPEND _IMPORT_CHECK_FILES_FOR_gqcp "/usr/local/gqcp/lib/libgqcp.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
