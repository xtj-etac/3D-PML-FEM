#ifndef FMODEL_H
#define FMODEL_H

using namespace std;

#define N_MODEL_X 20  
#define N_MODEL_Y 20  
#define N_MODEL_Z 20   
#define NPML_X  6
#define NPML_Y  6
#define NPML_Z  6
#define Nair NPML_Z+4 
#define main_dx 100
#define main_dy 100
#define main_dz 100
#define PML_DX 5e-2
#define PML_DY 5e-2
#define PML_DZ 5e-2
#define find_kappa 3
#define find_sigma 0.84
#define find_alpha 2.08


class Fmodel
{
public:
    double epsilo0,eps0,miu0,omega;
    int NN,m;
    int NX,NY,NZ;
    double sigma_max,kappa_max,alpha_max;
    double *A_X=nullptr,*B_Y=nullptr,*C_Z=nullptr;
    int NE,NL;

    double *xekappa=nullptr,*xesigma=nullptr,*xealpha=nullptr;
    double *yekappa=nullptr,*yesigma=nullptr,*yealpha=nullptr;
    double *zekappa=nullptr,*zesigma=nullptr,*zealpha=nullptr;
    complex<double> *Sex_origin=nullptr,*Sey_origin=nullptr,*Sez_origin=nullptr;

    vector<vector<int>> ME;
    vector<array<double,3>> rho,rho_b,ms;
    vector<double> alpha_S,alpha_D,alpha_L;
    vector<double> alpha_b_S,alpha_b_D,alpha_b_L;
    vector<double> ms_S,ms_D,ms_L;

public:
    Fmodel(double f);
    ~Fmodel();
    
    void settMatrix();
    void init_Fmodel();
   
    void freeFmodel();
};
#endif
