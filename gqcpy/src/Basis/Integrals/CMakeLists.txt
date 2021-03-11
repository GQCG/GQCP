list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/FunctionalPrimitiveIntegralEngine_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IntegralCalculator_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/McMurchieDavidsonCoefficient_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OneElectronIntegralEngine_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
