list(APPEND test_target_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/TransformationMatrix_test.cpp
        )

add_subdirectory(Integrals)
add_subdirectory(ScalarBasis)
add_subdirectory(SpinorBasis)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
