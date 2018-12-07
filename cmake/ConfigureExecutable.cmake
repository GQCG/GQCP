# In this file, we will provide a function that takes the name of an executable, and includes/links the appropriate libraries to it.


function(configure_executable EXECUTABLE_NAME)

    target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${LIBRARY_NAME})

endfunction(configure_executable)
