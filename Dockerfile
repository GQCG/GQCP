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
    && bash Miniconda3-latest-Linux-x86_64.sh -p /usr/local/miniconda3 -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 
RUN conda --version

ENV LIBINT_DATA_PATH=/usr/local/miniconda3/share/libint/2.4.2/basis
ARG LIBINT_DATA_PATH=/usr/local/miniconda3/share/libint/2.4.2/basis
RUN conda install -c conda-forge -c intel -c gqcg openmp=8.0.1 cmake=3.14 boost-cpp=1.69.0 eigen=3.3.7 blas=2.15 mkl=2019.0 pybind11=2.4.3 benchmark numpy=1.18.1 jupyter matplotlib pandas h5py libint=2.4.2 cint=3.0.17 mkl-include=2019.0 mkl-static=2019.0 intel-openmp=2019.0 

COPY . .
RUN mkdir build && cd build && cmake .. \
    -DCMAKE_PREFIX_PATH=/usr/local/miniconda3 \
    -DCMAKE_INSTALL_PREFIX=/usr/local/miniconda3 \
    -DBUILD_TESTS=TRUE \
    -DBUILD_PYTHON_BINDINGS=TRUE 
RUN cd build && make -j2 VERBOSE=1 && env CTEST_OUTPUT_ON_FAILURE=1 make test && make install

RUN ldconfig

ENTRYPOINT bash
