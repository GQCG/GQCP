FROM gqcg/gqcp-dev:latest

COPY . .
RUN mkdir build && cd build && cmake .. \
    -DCMAKE_PREFIX_PATH=/usr/local/miniconda3 \
    -DCMAKE_INSTALL_PREFIX=/usr/local/miniconda3 \
    -DBUILD_TESTS=TRUE \
    -DBUILD_PYTHON_BINDINGS=TRUE 
RUN cd build && make -j2 VERBOSE=1 && env CTEST_OUTPUT_ON_FAILURE=1 make test && make install

RUN ldconfig

ENTRYPOINT bash
