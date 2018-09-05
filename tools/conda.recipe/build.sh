#!/usr/bin/env bash

mkdir build_for_release
cd build_for_release
cmake ..
make && make test && sudo make install
