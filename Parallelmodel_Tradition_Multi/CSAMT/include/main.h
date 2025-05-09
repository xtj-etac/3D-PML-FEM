#include<iostream>
#include<stdio.h>
#include"mpi.h"
#include<omp.h>
#include<math.h>
#include<unistd.h>
#include<time.h>
#include<fstream>
#include"eigen3/Eigen/Sparse"
#include<eigen3/Eigen/Dense>
#include<complex.h>
#include"./fmodel.h"
#include"./v0Assembly.h"
#include"./v1Assembly.h"
#include"./v2Assembly.h"
#include"./rhoAndpha.h"
#include <petscmat.h>
#include <petscsys.h>
//#include"../include/v0Assembly.h"
//#include"../include/v1Assembly.h"
//#include"../include/v2Assembly.h"
//#include"../include/rhoAndpha.h"
#include<vector>
#include<array>
#include<lapacke.h>
#include<petsc.h>
