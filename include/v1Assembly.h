#ifndef V1ASSEMBLY_H
#define V1ASSEMBLY_H
#include<cmath>
#include <petsc.h>
extern "C"
{
    #include"superlu_zdefs.h"
}
using namespace std;

class v1Assembly
{
public:
    double w,mu0;
    int polarization;
    int m_len;
    int datasetSize;
    typedef Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> MatrixX;

    MatrixX  EpxA,EpyA,EpzA;
    MatrixX  linshi;
    Eigen::MatrixXd HzA;

    vector<vector<complex<double>>> Exx1;
    vector<vector<complex<double>>> Eyy1;
    vector<vector<complex<double>>> Ezz1;
    vector<vector<complex<double>>> Hxx1;
    vector<vector<complex<double>>> Hyy1;
    vector<vector<complex<double>>> Hzz1;
   
public:
    v1Assembly(Fmodel fmodel,double f);
    ~v1Assembly();

    void fDirBdaries(int groupid);
    void importdata(string f1,string f2,MatrixX &E);
    void creatv1(v0Assembly v0Assembly,Fmodel fmodel,MPI_Comm curComm);

};


#endif
