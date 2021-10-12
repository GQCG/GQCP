FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y wget \
    build-essential \
    clang \
    git \ 
    gdb \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && apt-get autoremove -y

ENV PATH="/usr/local/miniconda3/bin:${PATH}"
ARG PATH="/usr/local/miniconda3/bin:${PATH}"
RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -p /usr/local/miniconda3 -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda --version

ENV LIBINT_DATA_PATH=/usr/local/miniconda3/share/libint/2.6.0/basis
ARG LIBINT_DATA_PATH=/usr/local/miniconda3/share/libint/2.6.0/basis
RUN conda install -c conda-forge -c intel -c gqcg openmp=8.0.1 cmake=3.20.2 boost-cpp=1.75.0 eigen=3.3.9 blas=2.15 pybind11=2.6.2 benchmark numpy=1.20.3 jupyter=1.0.0 libint=2.6.0 cint=3.0.17 mkl=2020.4 mkl-include=2020.4 mkl-static=2020.4 intel-openmp=2020.3 doxygen=1.9.1 graphviz=2.47.2 jupyter-book=0.10.2 libblas=*=*mkl

RUN ldconfig

ENTRYPOINT bash
