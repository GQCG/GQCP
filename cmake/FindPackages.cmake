# In this CMake file, we will find all required packages


# Find the Boost package and some required components
find_package(Boost REQUIRED COMPONENTS system thread program_options)

# Find Eigen3
find_package(Eigen3 3.3.4 REQUIRED)
