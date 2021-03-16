# Using GQCP in an external package

External packages can be linked to `gqcp` using CMake. In order to find and load `gqcp` and its dependencies, CMake must look into the appropriate search path where all `cmake` files of the dependencies are installed.

```bash
list(APPEND CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH})
```

When Cmake tries to find the `gqcp` package, it will look into the given directory to find and load the package.

```bash
find_package(gqcp REQUIRED)
```

Once `gqcp` is found, the external project can be linked to use its functionalities.

```bash
target_link_libraries(${EXTERNAL_LIBRARY} PUBLIC GQCP::gqcp)
```