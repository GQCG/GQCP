conda activate gqcg_dev

cmake .. -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DCMAKE_INSTALL_PREFIX=~/.local -DBUILD_TESTS=TRUE -DBUILD_PYTHON_BINDINGS=TRUE

make -j2 
make -j2 gqcp
hercompileer enkel die test die mij interesseren
make test / make tests
make install
