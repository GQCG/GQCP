FROM ubuntu:18.04

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
    && mkdir /usr/local/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda --version

ENV LIBINT_DATA_PATH=/usr/local/miniconda3/share/libint/2.3.1/basis
ARG LIBINT_DATA_PATH=/usr/local/miniconda3/share/libint/2.3.1/basis
RUN conda install -c conda-forge openmp=8.0.1 cmake=3.14 boost-cpp=1.69.0 eigen=3.3.7 blas=2.15 mkl=2019.0 pybind11=2.4.3 benchmark numpy=1.18.1 jupyter
RUN conda install -c gqcg libint=2.3.1 cint=3.0.17
RUN conda install -c intel mkl-include=2019.0 mkl-static=2019.0 intel-openmp=2019.0

RUN ldconfig

ENTRYPOINT bash
