# Parse all *.in-files


# Parse version.hpp.in and gqcp.hpp.in
configure_file(${CMAKE_SOURCE_DIR}/cmake/Parse/version.hpp.in
               ${PROJECT_INCLUDE_FOLDER}/version.hpp @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/cmake/Parse/gqcp.hpp.in
               ${PROJECT_INCLUDE_FOLDER}/gqcp.hpp @ONLY)

# Parse Config.cmake.in and ConfigVersion.cmake.in
configure_file(${CMAKE_SOURCE_DIR}/cmake/Parse/Config.cmake.in
               ${CMAKE_SOURCE_DIR}/cmake/Parse/${PROJECT_NAME}Config.cmake @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/cmake/Parse/ConfigVersion.cmake.in
               ${CMAKE_SOURCE_DIR}/cmake/Parse/${PROJECT_NAME}ConfigVersion.cmake @ONLY)
