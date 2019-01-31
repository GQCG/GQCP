# 1. Configure the library

add_library(${LIBRARY_NAME} ${LIBRARY_TYPE} ${PROJECT_SOURCE_FILES} ${PROJECT_INCLUDE_FILES})


# Include this project's headers, and Spectra (see: https://github.com/arvidn/libtorrent/issues/3101#issuecomment-396195787)
target_include_directories(${LIBRARY_NAME} PUBLIC
        $<BUILD_INTERFACE:${PROJECT_INCLUDE_FOLDER}>
        $<INSTALL_INTERFACE:${PROJECT_INSTALL_INCLUDE_FOLDER}>)

# Include boost
target_include_directories(${LIBRARY_NAME} PUBLIC ${Boost_INCLUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PUBLIC ${Boost_LIBRARIES})

# Include Eigen
target_link_libraries(${LIBRARY_NAME} PUBLIC Eigen3::Eigen)

# Include Libint2
target_include_directories(${LIBRARY_NAME} PUBLIC ${Libint2_INCLUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PUBLIC ${Libint2_LIBRARIES})

# Include Spectra
target_include_directories(${LIBRARY_NAME} PRIVATE ${Spectra_INCLUDE_DIRS})

# Include MKL
if (USE_MKL)
    target_include_directories(${LIBRARY_NAME} PUBLIC ${MKL_INCLUDE_DIRS})
    target_link_libraries(${LIBRARY_NAME} PUBLIC ${MKL_LIBRARIES})
endif()



# 2. Install the library

# To specify that the library target should also be exported, we add the EXPORT option. This is used in conjuction with the install(EXPORT) command below
install(TARGETS ${LIBRARY_NAME}
        EXPORT ${LIBRARY_NAME}
        LIBRARY DESTINATION ${LIBRARY_INSTALL_DIR})


# Install the header files
install(DIRECTORY ${PROJECT_INCLUDE_FOLDER}/ DESTINATION ${INCLUDE_INSTALL_DIR})


# Export the target library into a ${PROJECT_NAME}Targets.cmake file.
# The file gqcpConfig.cmake includes this file, to be able to use this library with a find_package(template X.Y.Z) call
install(EXPORT ${LIBRARY_NAME}
        DESTINATION ${CMAKE_INSTALL_DIR}
        FILE ${PROJECT_NAME}Targets.cmake)

# Install Config.cmake and ConfigVersion.cmake
install(FILES
            ${CMAKE_SOURCE_DIR}/cmake/Parse/${PROJECT_NAME}Config.cmake
            ${CMAKE_SOURCE_DIR}/cmake/Parse/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_DIR})
