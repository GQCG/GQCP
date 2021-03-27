list(APPEND test_target_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/DenseVectorizer_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ImplicitMatrixSlice_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ImplicitRankFourTensorSlice_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LeviCivitaTensor_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Matrix_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SquareMatrix_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SquareRankFourTensor_test.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Tensor_test.cpp
)

set(test_target_sources ${test_target_sources} PARENT_SCOPE)
