
pybind11_add_module(gqcpy MODULE bindings.cpp)
target_link_libraries(gqcpy PUBLIC gqcp)