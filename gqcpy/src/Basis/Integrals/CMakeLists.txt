list(APPEND python_bindings_sources
    ${CMAKE_CURRENT_SOURCE_DIR}/FunctionalOneElectronPrimitiveIntegralEngine_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/FunctionalTwoElectronPrimitiveIntegralEngine_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HermiteCoulombIntegral_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/IntegralCalculator_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/McMurchieDavidsonCoefficient_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/OneElectronIntegralEngine_bindings.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TwoElectronIntegralEngine_bindings.cpp
)

set(python_bindings_sources ${python_bindings_sources} PARENT_SCOPE)
