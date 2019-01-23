# Link the dependencies to the target library

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

# Include MKL (optional)
if (MKL_FOUND)
    target_include_directories(${LIBRARY_NAME} PUBLIC ${MKL_INCLUDE_DIRS})
    target_link_libraries(${LIBRARY_NAME} PUBLIC ${MKL_LIBRARIES})
endif()

# Generate documentation
if (DOXYGEN_FOUND)
    set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(docs ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            VERBATIM)
endif (DOXYGEN_FOUND)
