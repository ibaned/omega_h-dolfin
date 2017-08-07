#!/bin/bash -ex
cmake $HOME/src/omega_h \
-DCMAKE_INSTALL_PREFIX:PATH=$HOME/install/gcc/omega_h-dolfin \
-DCMAKE_CXX_COMPILER:FILEPATH=$HOME/install/gcc/mpich/bin/mpicxx \
-DBUILD_SHARED_LIBS:BOOL=ON \
-DOmega_h_USE_MPI:BOOL=ON \
-DOmega_h_USE_ZLIB:BOOL=ON \
-DOmega_h_USE_DOLFIN:BOOL=ON \
-DDOLFIN_PREFIX=$HOME/install/gcc/dolfin \
-DBOOST_ROOT=$HOME/install/gcc/boost \
-DBUILD_TESTING:BOOL=ON \
2>&1 | tee config_log
