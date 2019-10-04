list(APPEND python_test_target_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/fci_binding_test.py  
        ${CMAKE_CURRENT_SOURCE_DIR}/fukui_dyson_binding_test.py  
        ${CMAKE_CURRENT_SOURCE_DIR}/hubbard_binding_test.py  
        ${CMAKE_CURRENT_SOURCE_DIR}/molecule_test.py  
        ${CMAKE_CURRENT_SOURCE_DIR}/mullikenConstrained_binding_test.py  
        ${CMAKE_CURRENT_SOURCE_DIR}/oo-doci_binding_test.py  
        )

set(python_test_target_sources ${python_test_target_sources} PARENT_SCOPE)
