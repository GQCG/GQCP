#!/usr/bin/env bash

CC="ccache $CC"
CXX="ccache $CXX"

mkdir build_release
cd build_release
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DPROJECT_INSTALL_DIR=${PREFIX} -DLIBRARY_TYPE=SHARED -DUSE_MKL=ON -DBUILD_DOCS=ON ..
make && make install
