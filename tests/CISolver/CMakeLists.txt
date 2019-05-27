list(APPEND test_target_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/CISolver_DOCI_Davidson_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/CISolver_DOCI_Dense_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/CISolver_FCI_Davidson_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/CISolver_FCI_Dense_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/CISolver_Hubbard_Davidson_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/CISolver_Hubbard_Dense_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/CISolver_test.cpp
        )
set(test_target_sources ${test_target_sources} PARENT_SCOPE)
