#ifndef V2ASSEMBLY_H
#define V2ASSEMBLY_H
#include<cmath>
#include <petsc.h>
//#include"slu_zdefs.h"
extern "C"
{
    #include"superlu_zdefs.h"
}
using namespace std;

class v2Assembly
{
public:
    double w,mu0;
    int polarization;
    typedef Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> MatrixX;
    int m_len;
    MatrixX  EpxB,EpyB,EpzB;
    MatrixX  linshi;
    Eigen::MatrixXd HzB;
    //Eigen::MatrixXd HzB;
    //int m_len;
    //double *m_f;
    vector<vector<complex<double>>> Exx2;
    vector<vector<complex<double>>> Eyy2;
    vector<vector<complex<double>>> Ezz2;
    vector<vector<complex<double>>> Hxx2;
    vector<vector<complex<double>>> Hyy2;
    vector<vector<complex<double>>> Hzz2;
   
public:
    v2Assembly(Fmodel fmodel,double f);
    ~v2Assembly();

    void fDirBdaries(int groupid);
    void importdata(string f1,string f2,MatrixX &E);
    void creatv2(v0Assembly v0Assembly,Fmodel fmodel,MPI_Comm curComm);

};


#endif
