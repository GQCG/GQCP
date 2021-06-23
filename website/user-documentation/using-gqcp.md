
# Using GQCP(y)

After you have installed GQCP, you can either interface directly with GQCP [through C++](#using-gqcp-through-cpp) or use its [Python bindings](#using-gqcp-through-python-bindings-gqcpy).

## Using GQCP through CPP

### Hello GQCP: source file

Given the following XYZ file

```bash
3

O          0.00000       -0.07579        0.00000
H          0.86681        0.60144        0.00000
H         -0.86681        0.60144        0.00000
```

a prototypical `Hello GQCP` driver `driver.cpp` is given by

```C++
#include <iostream>
#include <gqcp.hpp>

int main(int argc, char * argv [] ) {
    const auto water_molecule = GQCP::Molecule::ReadXYZ("../water.xyz");  // creates a neutral molecule
    std::cout << "Hello GQCP!" << std::endl;
    std::cout << water_molecule << std::endl;
    return 0;
}
```

### Hello GQCP: build stage

We have packaged `GQCP` into its own CMake modules, such that it can be easily found

```bash
find_package(gqcp REQUIRED)
```
and linked

```bash
target_link_libraries(${EXTERNAL_PROJECT_TARGET} PUBLIC GQCP::gqcp)
```

Given the `driver.cpp` source file detailed above, the following `CMakeLists.txt` should allow you to compile that code

```bash
cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(driver VERSION 0.1.0 LANGUAGES CXX)

find_package(gqcp REQUIRED)

add_executable(driver driver.cpp)
target_link_libraries(driver PUBLIC GQCP::gqcp)
```

given that you pass the location of GQCP during the CMake build procedure. On the Docker platform, this translates into

```bash
mkdir build && cd build
cmake .. -DCMAKE_PREFIX_PATH=/usr/local/miniconda3/
make -j${NCPUS}
```

with `${NCPUS}` the number of CPUs that can be reserved for compilation. After compilation, you can find the `driver.x` executable in the `build` folder.

## Using GQCP through Python bindings: GQCPy

In contrast to using `GQCP` through C++, using the provided Python bindings is as easy as

```python
import gqcpy
import numpy as np
```
