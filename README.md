# gqcp 0.0.1
[![Build Status](https://travis-ci.org/GQCG/gqcp.svg?branch=master)](https://travis-ci.org/GQCG/gqcp)

The Ghent Quantum Chemistry Package is a C++ library for electronic structure calculations.

## Dependencies

[![Boost Dependency](https://img.shields.io/badge/Boost-1.65.1+-000000.svg)](http://www.boost.org)
[![Eigen3 Dependency](https://img.shields.io/badge/Eigen-3.3.4+-000000.svg)](http://eigen.tuxfamily.org/index.php?title=Main_Page)
[![libint2 Dependency](https://img.shields.io/badge/libint-2.3.1+-000000.svg)](https://github.com/evaleev/libint)

gqcp uses the bases packaged with libint. Please set the `LIBINT_DATA_PATH` environment variable that contains these bases. In a default installation (of e.g. version v2.3.1), the data path is given by:
```
export LIBINT_DATA_PATH=/usr/local/libint/2.3.1/share/libint/2.3.1/basis
```

## Installation
To install this library:
1. clone the master branch, which contains the latest release

        https://github.com/GQCG/gqcp.git --branch master --single-branch
        cd gqcp

2. perform an out-of-source cmake build:

        mkdir build && cd build
        cmake -DINSTALLATION_PREFIX=prefix ..
        make && make test && sudo make install

    where
    * `prefix` is the installation prefix (defaulted to `/usr/local`) you want the library to be installed at:
        * the library `libgqcp.a` will be installed in `prefix/gqcp/lib`
        * the header files (and cmake files, see Usage) will be installed in `prefix/gqcp/include`


## Usage
Basic usage of this library can be found in the `tests` directory. If you use CMake in other projects, you can add the following CMake command to the CMakeLists.txt-file:

    find_package(gqcp 0.0.1)

CMake then provides the commands `gqcp_INCLUDE_DIRS` to be used in your `target_include_directories` and the library `gqcp` to be used in your `target_link_libraries`.
