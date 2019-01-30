include(${CMAKE_SOURCE_DIR}/cmake/Configure/ConfigureLibrary.cmake)

# Generate documentation
if (BUILD_DOCS)
    set(DOXYGEN_IN ${CMAKE_SOURCE_DIR}/cmake/Parse/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_BINARY_DIR}/Doxyfile)

    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)

    add_custom_target(docs ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            VERBATIM)
endif()
