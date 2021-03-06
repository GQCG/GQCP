list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/G1DM_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/G2DM_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/MixedSpinResolved2DMComponent_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Orbital1DM_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Orbital2DM_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PureSpinResolved2DMComponent_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinDensity1DM_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolved1DM_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolved2DM_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolved1DMComponent_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SpinResolved2DM_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
