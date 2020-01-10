list(APPEND python_test_target_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/DOCIRHF_test.py
        ${CMAKE_CURRENT_SOURCE_DIR}/FCI_test.py
        ${CMAKE_CURRENT_SOURCE_DIR}/FukuiDysonAnalysis_test.py
        ${CMAKE_CURRENT_SOURCE_DIR}/Hubbard_test.py
        ${CMAKE_CURRENT_SOURCE_DIR}/MullikenConstrainedFCI_test.py
        ${CMAKE_CURRENT_SOURCE_DIR}/OO-DOCI_test.py
)

set(python_test_target_sources ${python_test_target_sources} PARENT_SCOPE)
