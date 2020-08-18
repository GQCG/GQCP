---
id: documentation
title: Documentation overview
sidebar_label: Overview
---


## Documentation overview

An overview of the C++ documentation can be found [here](cpp_documentation.md), while the documentation for the Python bindings can be found [here](python_documentation.md)


## General structure of GQCP

GQCP has two main components: the C++ library `gqcp` and its associated Python bindings `gqcpy`. This is reflected in the folder structure of the GitHub repository, where we have, among other directories, the `gqcp` and `gqcpy` directories. As is common, the C++ library has `include`, `src` and `tests` folder, while the folder that contains the Python bindings only has a `src` folder and a `tests` folder.

Inside each of these folders, we have added the following structure:
- __Basis__: for everything related to spinors, spinor basis and transformations
- __Mathematical__: for general mathematical utilities, not directly related to quantum chemistry
- __Molecule__: for everything related to molecules and nuclei
- __ONVBasis__: for everything related to ONVs and bases for Fock (sub)spaces
- __Operator__: for everything related to first- and second-quantized operators
- __Processing__: for collecting everything that happens _after_ the determination of the optimal values of the electronic structure model's parameters
- __QCMethod__: for the determination of the optimal parameters of an electronic structure model
- __Utilities__: for collecting general utilities that do not belong elsewhere
