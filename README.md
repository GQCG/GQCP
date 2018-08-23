# temp vx.x.x
[![Build Status](https://travis-ci.org/GQCG/ci.svg?branch=master)](https://travis-ci.org/GQCG/ci)

A C++ library for performing configuration interaction (CI) calculations.

## Dependencies
[![Boost Dependency](https://img.shields.io/badge/Boost-1.65.1+-000000.svg)](http://www.boost.org)
[![Eigen3 Dependency](https://img.shields.io/badge/Eigen-3.3.4+-000000.svg)](http://eigen.tuxfamily.org/index.php?title=Main_Page)
[![libint2 Dependency](https://img.shields.io/badge/libint-2.3.1+-000000.svg)](https://github.com/evaleev/libint)
[![Spectra Dependency](https://img.shields.io/badge/Spectra-0.6.1+-000000.svg)](https://github.com/yixuan/spectra/)


[![cpputil Dependency](https://img.shields.io/badge/cpputil-1.5.0+-blue.svg)](https://github.com/GQCG/cpputil)
[![bmqc Dependency](https://img.shields.io/badge/bmqc-1.2.2+-blue.svg)](https://github.com/GQCG/bmqc)
[![libwint Dependency](https://img.shields.io/badge/libwint-3.0.0+-blue.svg)](https://github.com/GQCG/libwint)
[![hf Dependency](https://img.shields.io/badge/hf-3.0.0+-blue.svg)](https://github.com/GQCG/hf)
[![numopt Dependency](https://img.shields.io/badge/numopt-1.4.0+-blue.svg)](https://github.com/GQCG/numopt)


## Installation
To install this library:
1. clone the master branch, which contains the latest release

        https://github.com/GQCG/ci.git --branch master --single-branch
        cd ci

2. perform an out-of-source cmake build:

        mkdir build && cd build
        cmake -DINSTALLATION_PREFIX=prefix ..
        make && make test && sudo make install

    where
    * `prefix` is the installation prefix (defaulted to `/usr/local`) you want the library to be installed at:
        * the library `libci.a` will be installed in `prefix/ci/lib`
        * the header files (and cmake files, see Usage) will be installed in `prefix/ci/include`


## Usage
Basic usage of this library can be found in the `tests` directory. If you use CMake in other projects, you can add the following CMake command to the CMakeLists.txt-file:

    find_package(ci 1.2.1)

CMake then provides the commands `ci_INCLUDE_DIRS` to be used in your `target_include_directories` and the library `ci` to be used in your `target_link_libraries`.
