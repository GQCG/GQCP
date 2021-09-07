# Entry point for the CMake configuration of the C++ library
add_library(gqcp SHARED "")

target_compile_options(gqcp PUBLIC -m64 -pipe)
target_compile_options(gqcp PUBLIC "$<IF:$<CONFIG:Debug>,-g,-O2>")
target_compile_options(gqcp PUBLIC "$<$<CXX_COMPILER_ID:GNU>:-pthread>")
target_link_options(gqcp PUBLIC "$<$<CXX_COMPILER_ID:GNU>:-pthread>")
target_compile_features(gqcp PUBLIC cxx_std_14)

if (OPTIMIZE_FOR_NATIVE)
    target_compile_options(gqcp PRIVATE "$<IF:$<CXX_COMPILER_ID:Intel>,-xHost,-march=native>")
endif()


add_subdirectory(src)


target_include_directories(gqcp
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/gqcp/include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/gqcp/include>
        $<INSTALL_INTERFACE:${CMAKE_PREFIX_PATH}/include/gqcp>
)

target_link_libraries(gqcp
    PUBLIC
        gqcp_version
        Eigen3::Eigen
        Boost::boost
        Int2::Int2
        cint::cint
        MKL::MKL
)
 
# Add rt library on Linux environments
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    find_library(LIB_RT rt)
    target_link_libraries(gqcp PUBLIC ${LIB_RT})
endif()

target_compile_options(gqcp PUBLIC -DEIGEN_USE_MKL_ALL -DMKL_LP64)
