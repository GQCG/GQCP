# Troubleshooting guide

Can't build or use GQCP? Compiler errors? Other faults? You're at the right address. Here, we've compiled a list issues that have been encountered in the past.

## (HPC) EMT-instruction failure
If you encounter an EMT-instruction failure when using the Intel compilers on the HPC clusters, reduce the number of processors passed to the `make` command, e.g. replace `make -j 6` by `make -j 3`, or even `make -j 1`.

## (HPC) Slow compilation
Store and compile the code on `$VSC_SCRATCH` to achieve optimal performance. 

Still have a problem even after reading through the documentation and the troubleshooting guide? Create an [issue](https://github.com/GQCG/GQCP/issues/new/choose) and help us help you.
