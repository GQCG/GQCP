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


Since we (still) depend on Libint's basissets, after installation, you'll have to set the `LIBINT_DATA_PATH` environment variable to the folder that contains the libint bases. In a default installation (of Libint's version v2.4.2), the data path is given by:

```bash
export LIBINT_DATA_PATH=${CONDA_PREFIX}/share/libint/2.4.2/basis
```

You will have to either export this environment variable every time you activate the `gqcp` environment or (better) put this export in your .bashrc or (preferred) [add this environment variable to your virtual environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#saving-environment-variables).

## Installation through Docker

First [install Docker](https://docs.docker.com/get-docker/). Then pull our `GQCP` image to the infrastructure in question

```bash
docker pull gqcg/gqcp
```

If you bind the correct ports to the container based on that image

```bash
docker run -p 8888:8888 -it gqcg/gqcp
```

you can set up a Jupyter notebook in the Docker container

```bash
jupyter notebook --ip=0.0.0.0 --no-browser --allow-root
```

and access that Jupyter notebook in your local browser by looking in the terminal output for the link that begins with `127.0.0.1` 

```bash
To access the notebook, open this file in a browser:
        file:///root/.local/share/jupyter/runtime/nbserver-24-open.html
Or copy and paste one of these URLs:
        http://a696331630e9:8888/xxxxxx
     or http://127.0.0.1:8888/?xxxxxx
```

and copying that link in your browser.

