# Provide a function that takes the name of target executable, and includes/links the appropriate libraries to it


function(configure_executable EXECUTABLE_NAME)

    # Include this project
    target_include_directories(${EXECUTABLE_NAME} PRIVATE ${PROJECT_INCLUDE_FOLDER})
    target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${LIBRARY_NAME})

    # Include boost
    target_include_directories(${EXECUTABLE_NAME} PUBLIC ${Boost_INCLUDE_DIRS})
    target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${Boost_LIBRARIES})
    target_compile_definitions(${EXECUTABLE_NAME} PUBLIC -DBOOST_ALL_DYN_LINK)
    target_compile_definitions(${EXECUTABLE_NAME} PRIVATE _GLIBCXX_USE_CXX11_ABI=1)

    # Include Eigen
    target_link_libraries(${EXECUTABLE_NAME} PUBLIC Eigen3::Eigen)

    # Include libint2
    target_include_directories(${EXECUTABLE_NAME} PUBLIC ${Libint2_INCLUDE_DIRS})
    target_link_libraries(${EXECUTABLE_NAME} PUBLIC ${Libint2_LIBRARIES})

    # Include MKL (optional)
    if (MKL_FOUND)
        target_include_directories(${EXECUTABLE_NAME} PRIVATE ${MKL_INCLUDE_DIRS})
        target_link_libraries(${EXECUTABLE_NAME} PRIVATE ${MKL_LIBRARIES})
    endif()

endfunction(configure_executable)
