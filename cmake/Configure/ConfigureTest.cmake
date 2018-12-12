# In this file, we will provide a function that takes the name of an executable, and includes/links the appropriate libraries to it.


function(configure_test TEST_NAME)

    configure_executable(${TEST_NAME})


    # At the moment, there are no extra dependencies for the unit test executables

endfunction(configure_test)
