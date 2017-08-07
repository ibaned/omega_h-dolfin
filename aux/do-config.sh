#!/bin/bash -ex
cmake $HOME/src/dolfin-poisson \
-DCMAKE_CXX_COMPILER=$HOME/install/gcc/mpich/bin/mpicxx \
-DCMAKE_PREFIX_PATH="$HOME/install/gcc/omega_h-dolfin" \
-DBOOST_ROOT=$HOME/install/gcc/boost \
2>&1 | tee config_log
