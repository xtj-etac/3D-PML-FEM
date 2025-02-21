#ifndef V0ASSEMBLY_H
#define V0ASSEMBLY_H
#include <petscmat.h>
using namespace std;

class v0Assembly
{
public:
    double w,mu0;
    Mat v0;
    Mat K3;
    Mat K4;

public:
    v0Assembly(double f);
    ~v0Assembly();
    void creatv0(Fmodel fmodel,MPI_Comm row_comm);
};
#endif
