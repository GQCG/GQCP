# In this CMake file, we will find all required packages


# Find the Boost package
find_package(Boost REQUIRED REQUIRED program_options)

# Find Eigen3
find_package(Eigen3 3.3.4 REQUIRED)

# Find libint
find_package(libint2 REQUIRED)

# Find cpputil
find_package(cpputil 1.5.1 REQUIRED)

# Find numopt
find_package(numopt 1.5.1 REQUIRED)

# Find doxygen
if(BUILD_DOCS)
    find_package(Doxygen REQUIRED dot)
endif()
