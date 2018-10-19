# In this CMake file, we will include the headers and link to the necessary libraries


# Include this project's headers
target_include_directories(${LIBRARY_NAME} PRIVATE ${PROJECT_INCLUDE_FOLDER})

# Include the boost headers
target_include_directories(${LIBRARY_NAME} PUBLIC ${Boost_INCLUDE_DIRS})

# Include Eigen
target_link_libraries(${LIBRARY_NAME} PUBLIC Eigen3::Eigen)

# Include libint2
target_include_directories(${LIBRARY_NAME} PUBLIC ${libint2_INCLUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PUBLIC ${libint2_LIBRARIES})

# Include Spectra
target_include_directories(${LIBRARY_NAME} PUBLIC ${spectra_INCLUDE_DIRS})

# Include numopt
target_include_directories(${LIBRARY_NAME} PUBLIC ${numopt_INCLUDE_DIRS})
target_link_libraries(${LIBRARY_NAME} PUBLIC numopt)

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
	
