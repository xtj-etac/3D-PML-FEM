BLAS_DIR     =/thfs1/software/openblas/0.3.12-gcc9.3.0-openmp
LAPACK_DIR   =/thfs1/software/lapack/3.8.0-gcc9.3.0
PARMETIS_DIR =/thfs1/home/test653/shutian/YHMAX/parmetis-4.0.3
METIS_DIR    =/thfs1/software/spack/opt/linux-ubuntu20.04-aarch64/gcc-8.3.0/metis-5.1.0-65tfhw7
SUPERLU_DIST =/thfs1/home/test653/shutian/YHMAX-dist/superlu_dist-8.0.0
EIGEN_DIR     =/thfs1/home/test653/shutian/YHMAX/eigen-git-mirror-master
PETSC_DIR    =/thfs1/home/test653/shutian/YHMAX-dist/petsc-3.20.5



LIBS   += -lpetsc
LIBS   += -lopenblas
LIBS   += -llapack
LIBS   += -lparmetis -lmetis -lsuperlu_dist
LIBS   += -lgfortran -fopenmp -lgcc_s -ldl -lstdc++ -lX11 -lm  -lmpifort -lmpi

###############################################################
#LDFLAGS += -L${MPI_DIR}/lib
LDFLAGS += -L${BLAS_DIR}/lib
LDFLAGS += -L${LAPACK_DIR}/lib
LDFLAGS += -L${METIS_DIR}/lib
LDFLAGS += -L${PARMETIS_DIR}/build/_install/lib
LDFLAGS += -L${PETSC_DIR}/arch-linux-c-debug/lib -L${PETSC_DIR}/arch-linux-c-opt/lib
LDFLAGS += -L${SUPERLU_DIST}/_install/lib


###############################################################
#INCLUDES += -I${MPI_DIR}/include
INCLUDES += -I${METIS_DIR}/include
INCLUDES += -I${PARMETIS_DIR}/build/_install/include
INCLUDES += -I${PETSC_DIR}/include -I${PETSC_DIR}/arch-linux-c-opt/include
INCLUDES += -I${SUPERLU_DIST}/_install/include 
INCLUDES += -I${BLAS_DIR}/include/openblas
INCLUDES += -I${EIGEN_DIR}/_install/include

#LOADOPTS     =  -Wl,-rpath -Wl,-L/thfs1/software/mpich/mpi-n-gcc9.3.0/lib -lmpi -Wl,--enable-new-dtags -Wl,-rpath -Wl,/thfs1/software/mpich/mpi-n-gcc9.3.0/lib -Wl,--enable-new-dtags   

CC           = /thfs1/software/mpich/mpi-n-gcc9.3.0/bin/mpif90
CFLAGS 	     = -std=c99 -O3 -fopenmp -DNDEBUG -fPIC
##-fPIC
CXX          = /thfs1/software/mpich/mpi-n-gcc9.3.0/bin/mpicxx
CXXFLAGS     = -std=c++11  -O3 -fopenmp -DNDEBUG 
