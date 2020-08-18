---
id: documentation
title: Documentation overview
sidebar_label: Overview
---


## Why?

At GQCG, we strongly believe that every complex problem can be solved. Calculating the electronic structure is one of them. From the combination of our passion for quantum chemistry and open-source software, GQCP saw its first light in 2017 and has been growing ever since. 

But why does our quantum chemistry community need yet another electronic structure code? While other code bases might have better performance than GQCP, GQCP is modern at its core. GQCP is natively written in C++, so we have access to the most modern software techniques and compilers. GQCPY is its Python-bounded counterpart, so we fully embrace Python's present role as data manipulating language. Gone are the days of providing input files or writing executables, with GQCPY and Jupyter Notebooks you can naturally adapt a work flow that is both playful and systematic at the same time.


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


## Documentation overview

An overview of the C++ documentation can be found [here](cpp_documentation.md), while the documentation for the Python bindings can be found [here](python_documentation.md)