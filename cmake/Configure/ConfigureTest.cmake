# In this file, we will provide a function that takes the name of an executable, and includes/links the appropriate libraries to it.


function(configure_test TEST_NAME)

    configure_executable(${TEST_NAME})

    # Include Spectra
    target_include_directories(${TEST_NAME} PUBLIC ${Spectra_INCLUDE_DIRS})
endfunction(configure_test)
