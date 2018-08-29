# gqcg 0.0.1

A C++ library that provides the glue code for our other libraries/modules.


## Dependencies


## Installation
To install this library:
1. clone the master branch, which contains the latest release

        https://github.com/GQCG/gqcg.git --branch master --single-branch
        cd gqcg

2. perform an out-of-source cmake build:

        mkdir build && cd build
        cmake -DINSTALLATION_PREFIX=prefix ..
        make && make test && sudo make install

    where
    * `prefix` is the installation prefix (defaulted to `/usr/local`) you want the library to be installed at:
        * the library `libgqcg.a` will be installed in `prefix/gqcg/lib`
        * the header files (and cmake files, see Usage) will be installed in `prefix/gqcg/include`


## Usage
Basic usage of this library can be found in the `tests` directory. If you use CMake in other projects, you can add the following CMake command to the CMakeLists.txt-file:

    find_package(gqcg 0.0.1)

CMake then provides the commands `gqcg_INCLUDE_DIRS` to be used in your `target_include_directories` and the library `gqcg` to be used in your `target_link_libraries`.
