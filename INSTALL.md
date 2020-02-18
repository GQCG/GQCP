# Installation

## Prerequisites

Before installing GQCP, please make sure the following dependencies are available on your system:

[![Boost Dependency](https://img.shields.io/badge/Boost-<=1.69-000000.svg)](http://www.boost.org)
[![Eigen3 Dependency](https://img.shields.io/badge/Eigen-3.3.4+-000000.svg)](http://eigen.tuxfamily.org/index.php?title=Main_Page)
[![libint2 Dependency](https://img.shields.io/badge/libint-2.3.1+-000000.svg)](https://github.com/evaleev/libint)
[![libcint Dependency](https://img.shields.io/badge/gqcg_libcint-develop-000000.svg)](https://github.com/GQCG/libcint/tree/develop)

Note that we offer Conda packages for these installation requirements:

```bash
    conda env create -f environment.yml
    conda activate gqcg_dev
```

If you use your own installation of libint, please set the `LIBINT_DATA_PATH` environment variable to the folder that contains these bases. In a default installation (of e.g. version v2.3.1), the data path is given by:

```bash
   export LIBINT_DATA_PATH=/usr/local/libint/2.3.1/share/libint/2.3.1/basis
```

## CMake out-of-source build

For a default CMake build, the steps are the following:

1. clone the develop branch, which contains the latest release

    ```bash
        git clone https://github.com/GQCG/GQCP.git --branch develop --single-branch --recurse-submodules
        cd GQCP
    ```

2. perform an out-of-source build:

    ```bash
        mkdir build && cd build
        cmake .. (CMake options)
        make && make test && sudo make install
    ```

The possible CMake options are listed below.

## CMake options

In general, please set and pass the following options to the `cmake ..` command:

* `-DCMAKE_PREFIX_PATH=prefix_path`, with `prefix_path` the path to those libraries and includes that are not in default locations, but are grouped together.
   For instance, setting the prefix_path to `/usr/local` ensures that the folders `cmake`, `lib`, `lib64` and `include` can be found.
   You should set the prefix to libInt2, libCint, Eigen and Intel MKL.
   If you have chosen to use conda for some dependencies, `prefix_path` should be set to e.g. `/anaconda3/envs/.../`

* `-DCMAKE_INSTALL_PREFIX=prefix`, with `prefix` (defaulted to `/usr/local`) the installation prefix you want the library to be installed it. This option controls where the library is installed:
    * the header files will be installed in `prefix/include`
    * the compiled library will be installed in `prefix/lib`
    * drivers (optional) and benchmarks (optional) will be installed in `prefix/bin`
    * CMake target files will be installed in `prefix/cmake`

    We note that setting `CMAKE_INSTALL_PREFIX=~/.local` is preferred as this is also makes sure that the installed Python modules can be found automatically.


For this library, there are several extra options and configuration arguments you can pass to the `cmake ..` command:

* `-DBUILD_TESTS=TRUE` specifies that tests should be built and run.

* `-DBUILD_BENCHMARKS=TRUE` makes sure CMake adds the benchmark executables as targets. This uses [Google benchmark](https://github.com/google/benchmark), so make sure you have this installed if you wish to proceed with benchmarking on your system.

* `-DBUILD_DOCS=TRUE` specifies that documentation should be built using Doxygen, in which case Graphviz is required for UML generation. A custom `docs` target will then be configured by CMake, so that

        make docs

    compiles the documentation. After compilation, the HTML documentation can be found in the `docs/html` directory inside your out-of-source `build` directory. Navigating the documentation is easiest if you start with the `index.html` file.

* `-DBUILD_PYTHON_BINDINGS=TRUE` makes sure that selected pieces of the GQCP library can be called from Python. This uses [PyBind11](https://github.com/pybind/pybind11), so make sure you have this installed if you wish to use GQCPY on your system.

* `-DPYTHON_EXECUTABLE=python_executable` with `python_executable` the path to your preferred Python executable.
