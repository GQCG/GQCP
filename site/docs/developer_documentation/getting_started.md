---
id: developer_getting_started
title: Getting started
sidebar_label: Getting started
---

Before delving into GQCP's structure under the hood, let's make sure you're set up for development to GQCP.


## Cloning the repository

For GQCP, we use the [git flow](https://www.atlassian.com/git/tutorials/comparing-workflows/gitflow-workflow) workflow. All features are therefore based on the `develop` branch, which contains the latest developments.

```bash
git clone https://github.com/GQCG/GQCP.git --branch develop --single-branch --recurse-submodules
cd GQCP
```


## Installing dependencies

Please make sure the following dependencies are available on your system:

[![Boost Dependency](https://img.shields.io/badge/Boost-<=1.69-000000.svg)](http://www.boost.org)
[![Eigen3 Dependency](https://img.shields.io/badge/Eigen-3.3.4+-000000.svg)](http://eigen.tuxfamily.org/index.php?title=Main_Page)
[![libint2 Dependency](https://img.shields.io/badge/libint-2.3.1+-000000.svg)](https://github.com/evaleev/libint)
[![libcint Dependency](https://img.shields.io/badge/gqcg_libcint-develop-000000.svg)](https://github.com/GQCG/libcint/tree/develop)

You may install these manually, but please note that we offer a conda environment which contains these dependencies from the start. In the root directory of this repository, create the `gqcp_dev` conda environment from the `environment.yml`-file that we provide.
```bash
conda env create -f environment.yml
conda activate gqcp_dev
```

Since we (still) depend on Libint's basissets, after installation, you'll have to set the `LIBINT_DATA_PATH` environment variable to the folder that contains the libint bases. In a default installation (of Libint's version v2.3.1), the data path is given by:

```bash
export LIBINT_DATA_PATH=${CONDA_PREFIX}/share/libint/2.3.1/basis
```

You will have to either export this environment variable every time you activate the `gqcp` environment or (better) put this export in your .bashrc or (preferred) [add this environment variable to your virtual environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#saving-environment-variables).




### CMake out-of-source build

Right now, you're set up for development! Start a feature branch, edit some source files and create a pull request so we can merge your changes.

We have continuous integration set up, but in order to locally compile the library and run the test cases, we can use CMake to perform an out-of-source build. Out-of-source means that we're creating temporary folder in the root directory of this repository, usually just called `build`. Afterwards, we let CMake do its thing and then we run the test suite and possibly finalize with the installation of GQCP.
```bash
mkdir build && cd build
cmake .. (CMake options)
make && make test && sudo make install
```

Here `(CMAKE options)` still has to filled in by the options that GQCP supports.

* Specify the C compiler using `-DCMAKE_C_COMPILER=cc`, with `cc` the C-compiler that you would like to use.

* Specify the C++ compiler using `-DCMAKE_CXX_COMPILER=cxx`, with `cxx` the C++-compiler that you would like to use.

* In order to let CMake find those libraries and include headers that are not in default locations, please use the `-DCMAKE_PREFIX_PATH=prefix_path` options, with `prefix_path` the path to those libraries and includes (grouped together).

   If you have chosen to use conda for some dependencies, `prefix_path` should be set to e.g. `/anaconda3/envs/.../`

   If you use a custom installation of GQCP's dependencies, you should provide this option in order to let CMake find GQCP's dependencies (libInt2, libCint, Eigen and Intel MKL), which therefore should be installed in the subfolders of `prefix_path`. For instance, providing `-DCMAKE_PREFIX_PATH=/usr/local` that any subfolders may be found (e.g. the `cmake`, `lib`, `lib64` and `include` folders).

* In order to control where the library should be installed, you may specify `-DCMAKE_INSTALL_PREFIX=prefix`, with `prefix` the installation prefix you want the library to be installed it. When omitted, `prefix` defaults to `/usr/local`.
    * the header files will be installed in `prefix/include`
    * the compiled library will be installed in `prefix/lib`
    * drivers (optional) and benchmarks (optional) will be installed in `prefix/bin`
    * CMake target files will be installed in `prefix/cmake`

    We should note that setting `CMAKE_INSTALL_PREFIX=~/.local` is preferred as this is also makes sure that the installed Python modules can be found automatically.

* `-DBUILD_TESTS=TRUE` specifies that tests should be built and run.

* `-DBUILD_BENCHMARKS=TRUE` makes sure CMake adds the benchmark executables as targets. This uses [Google benchmark](https://github.com/google/benchmark), so make sure you have this installed if you wish to proceed with benchmarking on your system.

* `-DBUILD_DOCS=TRUE` specifies that the API documentation for the C++ library should be built using Doxygen, in which case Graphviz is required for UML generation. A custom `docs` target will then be configured by CMake, so that

        make docs

    compiles the documentation. After compilation, the HTML documentation can be found in the `docs/html` directory inside your out-of-source `build` directory. Navigating the documentation is easiest if you start with the `index.html` file.

* `-DBUILD_PYTHON_BINDINGS=TRUE` makes sure that selected pieces of the GQCP library can be called from Python. This uses [pybind11](https://github.com/pybind/pybind11), so make sure you have this installed if you wish to use `gqcpy` on your system.

* `-DPYTHON_EXECUTABLE=python_executable` with `python_executable` the path to your preferred Python executable.

* `-DPYTHON_LIBRARY=python_library` with `python_library` the path to the libraries that support your preferred Python executable.


As you can see, there are a lot of options that can (and should) be passed to CMake. For quick reference, here's a command that should work most of the time. In your out-of-source build directory, initialize CMake, make all targets, run all tests and install the library:

```bash
cmake .. -DCMAKE_PREFIX_PATH=${CONDA_PREFIX} \
    -DCMAKE_INSTALL_PREFIX=~/.local \
    -DBUILD_TESTS=TRUE \ 
    -DBUILD_PYTHON_BINDINGS=TRUE \
    -DPYTHON_EXECUTABLE=${CONDA_PREFIX}/bin/python \ 
    -DPYTHON_LIBRARY=${CONDA_PREFIX}/lib/libpython3.8.a
make -j 4 && make test && make install
```

`make -j 4` will make sure that 4 targets will be built at once. In order to increase compilation speed, you may increase this number, related to the number of cores you have available on your machine.
