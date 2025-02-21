# 3D-PML-FEM
## title:
Parallel finite element forward modeling of 3-D magnetotelluric conductivity and permeability anisotropy with coupled PML boundary conditions  

This code is for the coupled PML (Perfectly Matched Layer) three-dimensional magnetotelluric anisotropic parallel finite element forward modeling, which includes a single-frequency point parallel test case for model 3 (Figure 11 in the paper). The 'data' directory contains the input files, and the 'outfile' directory contains the output files, which output phase and apparent resistivity.The test time output containing a test case is in the log-1.out file. Please note that due to time constraints, this code is only an initial version.  
## Name of the code/library:
3D-PML-FEM     
## Contact:   
e-mail xiaotiaojie@nudt.edu.cn    
## Hardware requirements:   
Tianhe-2 Supercomputer   



# Getting Started
## Installation
Clone  
$ git clone https://github.com/  
$ cd PML

## Configure your environment
1、You will need to install the following external libraries:  
gcc9.3.0、mpich、superlu_dist-8.0.0、petsc-3.20.5 

2、Then you need to modify the make.inc file: 
BLAS_DIR     =/Your path/openmp  
LAPACK_DIR   =/Your path/gcc9.3.0  
PARMETIS_DIR =/Your path/parmetis-4.0.3  
METIS_DIR    =/Your path/metis  
SUPERLU_DIST =/Your path/superlu_dist-8.0.0  
EIGEN_DIR     =/Your path/eigen-git-mirror-master  
PETSC_DIR    =/Your path/petsc-3.20.5  

3、Compile your code    
$ make

4、Run your code  
$ sh mybatch.sh X1 X2 (X1 represents the number of nodes used at run time and X2 represents the total number of processes used at run time)  

$notice    
You need to select the mpi run command based on your environment, as shown below:  
mpiexe -n numpro(Number of processes) ./main

# Citation  
If you use this code for your research, please cite the paper.

