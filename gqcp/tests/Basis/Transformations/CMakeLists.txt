list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/JacobiRotation_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OrbitalRotationGenerators_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SimpleTransformationMatrix_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
