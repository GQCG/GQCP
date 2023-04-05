cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(gqcp-link LANGUAGES CXX)

set(CMAKE_MODULE_PATH "${GQCP_PREFIX_PATH}/cmake/;${GQCP_INSTALL_PREFIX}/cmake")
find_package(gqcp REQUIRED)

add_executable(link link.cpp)
target_link_libraries(link PUBLIC GQCP::gqcp)
