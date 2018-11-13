#!/usr/bin/env bash

mkdir build_release
cd build_release
cmake -DCMAKE_INSTALL_PREFIX=${PREFIX} -DPROJECT_INSTALL_DIR=${PREFIX} -DLIBRARY_TYPE=SHARED -DBUILD_DOCS=ON ..
make && make test && make install
