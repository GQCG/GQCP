list(APPEND test_python
        ${CMAKE_CURRENT_SOURCE_DIR}/hubbard_binding_test.py  
        ${CMAKE_CURRENT_SOURCE_DIR}/molecule_test.py  
        ${CMAKE_CURRENT_SOURCE_DIR}/mullikenConstrained_binding_test.py  
        )
set(test_python ${test_python} PARENT_SCOPE)
