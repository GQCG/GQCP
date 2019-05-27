list(APPEND test_target_sources
        ${CMAKE_CURRENT_SOURCE_DIR}/FockSpace_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/FrozenFockSpace_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/FrozenProductFockSpace_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ONV_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ProductFockSpace_test.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SelectedFockSpace_test.cpp
        )
set(test_target_sources ${test_target_sources} PARENT_SCOPE)
