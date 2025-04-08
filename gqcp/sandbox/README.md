The compiler has been set up such that any `xxx.cpp` file inside this `sandbox` folder is compiled into an executable in `build/gqcp/sandbox/` once it has been added to `sandbox_target_sources` inside `sandbox/CmakeLists.txt`, line 4.
To compile, just run
```bash
cmake .. -BUILD_SANDBOX=TRUE && make -j 4
```
If GQCP has already been compiled before inside the `build/` directory (using the [command specified in the website developer documentation](https://gqcg.github.io/GQCP/developer-documentation/getting-started.html#cmake-options-quick-reference)), only the _changed_ files will be recompiled. This means, if you only make changes to files inside `sandbox/`, then compilation will be very fast, and you will not have to wait long to try out a small test.