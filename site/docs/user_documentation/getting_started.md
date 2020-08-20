---
id: user_getting_started
title: Getting started
sidebar_label: Getting started
---

Before discussing GQCP's capabilities, let's make sure you're all set up to use GQCP.


## Installation through conda

The quickest way to install GQCP is through conda. The following commands create a conda environment called `gqcp` and installs GQCP and its dependencies in it.

```bash
conda create --name gqcp
source activate gqcp
conda install -c gqcg -c intel -c conda-forge gqcp
```


Since we (still) depend on Libint's basissets, after installation, you'll have to set the `LIBINT_DATA_PATH` environment variable to the folder that contains the libint bases. In a default installation (of Libint's version v2.3.1), the data path is given by:

```bash
export LIBINT_DATA_PATH=${CONDA_PREFIX}/share/libint/2.3.1/basis
```

You will have to either export this environment variable every time you activate the `gqcp` environment or (better) put this export in your .bashrc or (preferred) [add this environment variable to your virtual environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#saving-environment-variables).
