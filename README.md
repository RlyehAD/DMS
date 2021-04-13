
               Welcome to the official version of DMS!

DMS (Deductive Multiscale Simulator) an open-source multiscale all-atom 
simulator that is based on the GROMACS molecular dynamics source code.

DMS is free software, distributed under the GNU General
Public License, version 2.1. DMS uses a substantial portion of the source
code of the molecular dynamics package GROMACS. If you want to distribute 
a modified version or use part of DMS in your own program, remember that 
the entire project must be licensed according to the requirements of the LGPL 
v2.1 license under which you received this copy of DMS. We request that it 
must clearly be labeled as derived work. It should not use the name 
"official DMS", and make sure support questions are directed to you instead 
of the DMS developers.

                               * * * * *
## Acknowledgement

The development of DMS is mainly funded by academic research grants. 
To help us fund development, we humbly ask that you cite the DMS papers:

* Multiscale Factorization Method for Simulating Mesoscopic Systems with 
  Atomic Precision
  A. A. Mansour, P. J. Ortoleva 
  J. Chem. Theory Comput., 2014, 10, 518
  DOI: 10.1021/ct400615a

* Energy transfer between a nanosystem and its host fluid: A multiscale 
  factorization approach
  Y. V. Sereda, J. M. Espinosa-Duran, P. J. Ortoleva
  J. Chem. Phys., 2014, 140, pp 074102
  DOI: 10.1063/1.4864200


## Compilation with conda on Linux
WARNING: Your conda enviroments may not be located in /miniconda3/envs. Please use 'conda info --envs' to check.
```
conda create -n dms -c conda-forge gcc_linux-64
conda activate dms
conda install -c conda-forge gxx_linux-64
conda install -c conda-forge cmake hdf5 petsc fftw boost
git clone git@github.com:CTCNano/DMS.git
cd DMS
mkdir build
cd build

PETSC_LIBRARIES=/home/USER_NAME/miniconda3/envs/dms/lib/libpetsc.so PETSC_INCLUDES=/home/USER_NAME/miniconda3/envs/dms/include   
cmake .. -DCMAKE_PREFIX_PATH=/miniconda3/envs/dms/ -DCMAKE_INSTALL_PREFIX=INSTALLATION_ADDRESS/dms   
-DCMAKE_C_COMPILER=/home/USER_NAME/miniconda3/envs/dms/bin/x86_64-conda_cos6-linux-gnu-gcc   
-DBUILD_SHARED_LIBS=ON -DGMX_DEFAULT_SUFFIX=OFF -DGMX_BINARY_SUFFIX=_dms -DGMX_LIBS_SUFFIX=_dms   
-DGMX_BUILD_MDRUN_ONLY=ON -DGMX_GPU=OFF -DGMX_MPI=ON -DGMX_OPENMP=OFF

make
make install
```
