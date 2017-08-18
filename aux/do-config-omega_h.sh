#!/bin/bash -ex
cmake $HOME/omega_h \
-DCMAKE_INSTALL_PREFIX:PATH=$HOME/omega_h/install \
-DCMAKE_CXX_COMPILER:FILEPATH=`which mpicxx` \
-DCMAKE_C_COMPILER:FILEPATH=`which mpicc` \
-DBUILD_SHARED_LIBS:BOOL=ON \
-DOmega_h_USE_MPI:BOOL=ON \
-DOmega_h_USE_ZLIB:BOOL=ON \
-DOmega_h_USE_DOLFIN:BOOL=ON \
2>&1 | tee config_log
