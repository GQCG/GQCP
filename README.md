# gqcp 0.1.0
[![Build Status](https://travis-ci.org/GQCG/gqcp.svg?branch=master)](https://travis-ci.org/GQCG/gqcp)

The Ghent Quantum Chemistry Package is a C++ library for electronic structure calculations.


## Dependencies

[![Boost Dependency](https://img.shields.io/badge/Boost-1.65.1+-000000.svg)](http://www.boost.org)
[![Eigen3 Dependency](https://img.shields.io/badge/Eigen-3.3.4+-000000.svg)](http://eigen.tuxfamily.org/index.php?title=Main_Page)
[![libint2 Dependency](https://img.shields.io/badge/libint-2.3.1+-000000.svg)](https://github.com/evaleev/libint)

[![cpputil Dependency](https://img.shields.io/badge/cpputil-1.5.0+-blue.svg)](https://github.com/GQCG/cpputil)
[![numopt Dependency](https://img.shields.io/badge/numopt-1.5.0+-blue.svg)](https://github.com/GQCG/numopt)

As gqcp uses the bassisets packaged with libint, please set the `LIBINT_DATA_PATH` environment variable to the folder that contains these bases. In a default installation (of e.g. version v2.3.1), the data path is given by:

    export LIBINT_DATA_PATH=/usr/local/libint/2.3.1/share/libint/2.3.1/basis


## Installation
To install this library:
1. clone the master branch, which contains the latest release

        https://github.com/GQCG/gqcp.git --branch master --single-branch
        cd gqcp

2. perform an out-of-source build:

        mkdir build && cd build
        cmake -DINSTALLATION_PREFIX=prefix ..
        make && make test && sudo make install

    where
    * `prefix` is the installation prefix (defaulted to `/usr/local`) you want the library to be installed at:
        * the library `libgqcp.a` will be installed in `prefix/gqcp/lib`
        * the header files (and cmake files, see Usage) will be installed in `prefix/gqcp/include`


## Usage
Basic usage of this library can be found in the `tests` directory, and common use cases are explained in the corresponding [Wiki page](https://github.com/GQCG/gqcp/wiki/Common-use-cases). If you use CMake in other projects, you can add the following CMake command in your project:

    find_package(gqcp 0.1.0)

CMake then provides the commands `gqcp_INCLUDE_DIRS` to be used in your `target_include_directories` and the library `gqcp` to be used in your `target_link_libraries`.


## Documentation
Documentation can be built with Doxygen (in which case Graphviz is required for UML generation). If you pass the following option to CMake

    cmake -DBUILD_DOCS=ON ..

a custom `docs` target will be configured. After the documentation is compiled through

    make docs

the HTML documentation can be found in the `docs/html` directory inside the `build` directory. Navigating the documentation is easiest if you start with the `index.html` file.
