# In this CMake file, we will find all required packages


# Find the Boost package
find_package(Boost REQUIRED)

# Find Eigen3
find_package(Eigen3 3.3.4 REQUIRED)

# Find cpputil
find_package(cpputil 1.5.0 REQUIRED)

# Find libint
find_package(libint2 REQUIRED)

# Find numopt
find_package(numopt 1.4.0 REQUIRED)
