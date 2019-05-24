list(APPEND test_target_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/ChemicalMatrix_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ChemicalRankFourTensor_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Matrix_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ScalarFunction_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SquareMatrix_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SquareRankFourTensor_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Tensor_test.cpp
        )

add_subdirectory(optimization)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)

