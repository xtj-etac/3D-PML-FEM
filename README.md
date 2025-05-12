# 3D-PML-FEM
## Author:  
Tiaojie Xiao;Tian Shu et al.
## Date:   
2024-11
## Language： 
C++  
If you use this code in your research, please do not use it for commercial purposes and please provide the source.

## Title:
Parallel finite element forward modeling of 3-D magnetotelluric conductivity and permeability anisotropy with coupled PML boundary conditions  

## Contact:   
e-mail xiaotiaojie@nudt.edu.cn    
## Hardware requirements:   
Tianhe-2 Supercomputer    

## Brief Introduction  
This code is for the coupled PML (Perfectly Matched Layer) three-dimensional magnetotelluric anisotropic parallel finite element forward modeling, which includes a single-frequency point parallel test case for model 3 (Figure 11 in the paper).  

This library contains the original data of all the experiments in the paper. 
There are 21 directorys in repository.The LibraryFile file contains the library files that need to be installed. The other directories all contain a complete example, code, experimental data, and running log information (except for the directory of Model 1).
Note that in the experiments of Model 1, we compare our results with those computed using the MATLAB code provided by Yang et al. to verify correctness.Due to MATLAB software licensing restrictions, we only provide the experimental data in the Model 1 directory. Note that due to the same code structure, we only provide standardized and detailed code comments in the model2_PML_single_air file.

The directory naming convention is as follows (note that some fields are not mandatory):  
name:      model_boundary_process_model type_main author of the program_permeability type  
            |       |        |        |                  |                     |
            |       |        |        |                  |                     |  
Field ID:   1       2        3        4                  5                     6


 

# Getting Started
## Installation
Clone  
$ git clone https://github.com/xtj-etac/3D-PML-FEM.git  
$ cd 3D-PML-FEM

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

