# In this file, we will provide a function that takes the name of an executable, and includes/links the appropriate libraries to it.


function(configure_test TEST_NAME)

    # Include this project
    target_include_directories(${TEST_NAME} PRIVATE ${PROJECT_INCLUDE_FOLDER})
    target_link_libraries(${TEST_NAME} PRIVATE ${LIBRARY_NAME})

    # Include boost
    target_include_directories(${TEST_NAME} PUBLIC ${Boost_INCLUDE_DIRS})
    target_link_libraries(${TEST_NAME} PUBLIC ${Boost_LIBRARIES})

    # Include Eigen
    target_link_libraries(${TEST_NAME} PUBLIC Eigen3::Eigen)

    # Include libint2
    target_include_directories(${TEST_NAME} PUBLIC ${Libint2_INCLUDE_DIRS})
    target_link_libraries(${TEST_NAME} PUBLIC ${Libint2_LIBRARIES})

    # Include Spectra
    target_include_directories(${TEST_NAME} PUBLIC ${spectra_INCLUDE_DIRS})

    # Include cpputil
    target_include_directories(${TEST_NAME} PRIVATE ${cpputil_INCLUDE_DIRS})
    target_link_libraries(${TEST_NAME} PRIVATE cpputil)

    # Include numopt
    target_include_directories(${TEST_NAME} PUBLIC ${numopt_INCLUDE_DIRS})
    target_link_libraries(${TEST_NAME} PUBLIC numopt)

    # Include MKL (optional)
    if (MKL_FOUND)
        target_include_directories(${TEST_NAME} PRIVATE ${MKL_INCLUDE_DIRS})
        target_link_libraries(${TEST_NAME} PRIVATE ${MKL_LIBRARIES})
    endif()

endfunction(configure_test)
