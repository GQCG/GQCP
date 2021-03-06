# Getting started


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

## Installation through Singularity

For HPC systems, [singularity](https://sylabs.io/docs/) offers a more secure fork of Docker. Singularity converts Docker images to Singularity images on the fly. For the [UGent HPC](https://www.ugent.be/hpc/en) this translates into first making sure [that you are able to download and use Singularity images](https://vlaams-supercomputing-centrum-vscdocumentation.readthedocs-hosted.com/en/latest/software/singularity.html)

```bash
mkdir $VSC_SCRATCH/containers
mkdir $VSC_SCRATCH/containers/cache
mkdir $VSC_SCRATCH/containers/tmp
```

then you can pull and convert the Docker image in the scratch space

```bash
SINGULARITY_CACHEDIR=$VSC_SCRATCH/containers/cache SINGULARITY_TMPDIR=$VSC_SCRATCH/containers/tmp SINGULARITY_PULLFOLDER=$VSC_SCRATCH/containers singularity pull docker://gqcg/gqcp
```


> **Note**: The `SINGULARITY_*` variables used above can also be set user wide in your `.bashrc`. 

After converting the Docker image to a Singularity file you can shell into the resulting `*.sif`

```bash
singularity shell $VSC_SCRATCH/containers/gqcp_latest.sif
```

You can launch a Jupyter notebook with the following command in the Singularity container

```bash
jupyter notebook --ip=0.0.0.0 --no-browser --allow-root
```

Note the node number of the cluster on which you are running the container, e.g.:

```
http://node2473.golett.os:8888/?token=a3bf3b8bbb14618586ff6716bfa8f51886ff14d0f5e2c458
 or http://127.0.0.1:8888/?token=a3bf3b8bbb14618586ff6716bfa8f51886ff14d0f5e2c458
```

and set up an SSH tunnel in a separate terminal, such that you can access that notebook from your local browser

```bash
ssh -J vsc{xxxxx}@login.hpc.ugent.be -L 8888:localhost:8888 vsc{xxxxx}@node{nnnn}.{cluster}.os
```

where `{xxxxx}` is your vsc number, `{nnnn}` is the node number mentioned above and `{cluster}` is the cluster on which your Jupyter notebook is running. After this tunnel has been set up, you can access your notebook through your local browser by using the `127.0.0.1:8888` link mentioned above:

```bash
http://127.0.0.1:8888/?token=a3bf3b8bbb14618586ff6716bfa8f51886ff14d0f5e2c458
```
