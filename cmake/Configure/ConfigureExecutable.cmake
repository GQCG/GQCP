# Provide a function that takes the name of target executable, and includes/links the appropriate libraries to it


function(configure_executable EXECUTABLE_NAME)

    # Include this project
    target_include_directories(${EXECUTABLE_NAME} PRIVATE ${PROJECT_INCLUDE_FOLDER})
    target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${LIBRARY_NAME})


endfunction(configure_executable)
