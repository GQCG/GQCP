# Entry point for the CMake configuration of the gqcpy Python bindings

# Add the sources for the Python bindings
set(python_bindings_sources)

list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/gqcpy.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/version.cpp
)

add_subdirectory(include)
add_subdirectory(src)


# Configure the Python bindings
pybind11_add_module(gqcpy MODULE ${python_bindings_sources})


target_include_directories(gqcpy
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
)


target_link_libraries(gqcpy PUBLIC gqcp)
set_target_properties(gqcpy
    PROPERTIES
        SUFFIX ".so"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/gqcpy"
)
configure_file(__init__.py.in ${CMAKE_BINARY_DIR}/gqcpy/__init__.py)
configure_file(setup.py.in ${CMAKE_BINARY_DIR}/gqcpy/setup.py)
