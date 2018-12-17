# In this CMake file, we will provide the necessary commands to make the install target for this library


# The target of this project is a library called "template" (${PROJECT_NAME}).
# To specify that this target should also be exported, we add the EXPORT option. This is used in conjuction with the install(EXPORT) command below
install(TARGETS ${LIBRARY_NAME}
        EXPORT ${LIBRARY_NAME} ${EXPORT_TYPE}
        LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})


# Parse the version.hpp.in (input - .in) file (this replaces all @VAR@ by their CMake value)
configure_file(${PROJECT_INCLUDE_FOLDER}/version.hpp.in ${PROJECT_INCLUDE_FOLDER}/version.hpp @ONLY)


# Install the header files
install(DIRECTORY ${PROJECT_INCLUDE_FOLDER}/ DESTINATION ${INCLUDE_INSTALL_DIR})


# Export the target library into a ${PROJECT_NAME}Targets.cmake file.
# The file ciConfig.cmake includes this file, to be able to use this library with a find_package(template X.Y.Z) call
install(EXPORT ${LIBRARY_NAME} DESTINATION ${CMAKE_INSTALL_DIR} FILE ${PROJECT_NAME}Targets.cmake)


# Parse Config.cmake.in
configure_file(${CMAKE_SOURCE_DIR}/cmake/Config.cmake.in ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake @ONLY)


# Parse ConfigVersion.cmake.in
configure_file(${CMAKE_SOURCE_DIR}/cmake/ConfigVersion.cmake.in ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}ConfigVersion.cmake @ONLY)


# Install Config.cmake and ConfigVersion.cmake
install(FILES
            ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake
            ${CMAKE_SOURCE_DIR}/cmake/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_DIR})
