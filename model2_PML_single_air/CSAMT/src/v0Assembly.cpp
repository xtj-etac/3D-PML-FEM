//Matrix element calculation and matrix assembly
#include "../include/main.h"
#include "../include/fmodel.h"
using namespace std;
typedef Eigen::SparseMatrix<complex<double>> MatrixX;
v0Assembly::v0Assembly(double f)
{
    w=2*M_PI*f, mu0=(4*pow(10,-7))*M_PI;
}

v0Assembly::~v0Assembly()
{
}

void v0Assembly::creatv0(Fmodel fmodel,MPI_Comm curComm)
{
    PetscLogDouble start,end1 ;
    PetscTime(&start);
    int NX=fmodel.NX;
    int NY=fmodel.NY;
    int NZ=fmodel.NZ;
    int NE=fmodel.NE;
    int NL=fmodel.NL;
    
    MPI_Request request;
    MPI_Status status;

    //Obtain the process ID and the size of the process set from the communication domain.
    int curRank, curSize;
    MPI_Comm_size(curComm, &curSize);
    MPI_Comm_rank(curComm, &curRank);
    //Create matrices v0, k3, and k4
    MatCreate(curComm, &v0);
    MatSetSizes(v0, PETSC_DECIDE, PETSC_DECIDE, NL,NL);
    MatSetFromOptions(v0);
    MatSetUp(v0);

    MatCreate(curComm, &K3);
    MatSetSizes(K3, PETSC_DECIDE, PETSC_DECIDE,  NL,NL);
    MatSetFromOptions(K3);
    MatSetUp(K3);

    MatCreate(curComm, &K4);
    MatSetSizes(K4, PETSC_DECIDE, PETSC_DECIDE,  NL,NL);
    MatSetFromOptions(K4);
    MatSetUp(K4);

    //Allocate the matrix unit calculations to different processes.
    int slide = floor((double)NE / curSize);
    int st_idx, ed_idx;
    st_idx = curRank*slide+1; 
    if(curRank != curSize-1)
        ed_idx = (curRank+1)*slide;
    else
        ed_idx = NE;

    for(int h=st_idx;h<=ed_idx;h++)
    {
            //Calculate the anisotropic conductivity based on the main shaft resistivity and Euler rotation angle
            Eigen::Matrix3d sigma_primary;
            sigma_primary<<1/fmodel.rho[h][0], 0, 0, 0 ,1/fmodel.rho[h][1],0, 0 ,0 ,1/fmodel.rho[h][2];

            Eigen::Matrix3d RzfuAlpha_S,RzfuAlpha_D,RzfuAlpha_L,RzAlpha_L,RzAlpha_S,RzAlpha_D,sigma_tensor;
            RzfuAlpha_S<<cos(-fmodel.alpha_S[h]),sin(-fmodel.alpha_S[h]),0,-sin(-fmodel.alpha_S[h]),cos(-fmodel.alpha_S[h]),0,0,0,1;
            RzfuAlpha_D<<1,0,0,0,cos(-fmodel.alpha_D[h]),sin(-fmodel.alpha_D[h]),0,-sin(-fmodel.alpha_D[h]),cos(-fmodel.alpha_D[h]);
            RzfuAlpha_L<<cos(-fmodel.alpha_L[h]),sin(-fmodel.alpha_L[h]),0,-sin(-fmodel.alpha_L[h]),cos(-fmodel.alpha_L[h]),0,0,0,1;
            RzAlpha_L<<cos(fmodel.alpha_L[h]),sin(fmodel.alpha_L[h]),0,-sin(fmodel.alpha_L[h]),cos(fmodel.alpha_L[h]),0,0,0,1;
            RzAlpha_D<<1,0,0,0,cos(fmodel.alpha_D[h]),sin(fmodel.alpha_D[h]),0,-sin(fmodel.alpha_D[h]),cos(fmodel.alpha_D[h]);
            RzAlpha_S<<cos(fmodel.alpha_S[h]),sin(fmodel.alpha_S[h]),0,-sin(fmodel.alpha_S[h]),cos(fmodel.alpha_S[h]),0,0,0,1;
            sigma_tensor=RzfuAlpha_S*RzfuAlpha_D*RzfuAlpha_L*sigma_primary*RzAlpha_L*RzAlpha_D*RzAlpha_S;

            double sigma_xx,sigma_xy,sigma_xz,sigma_yy,sigma_yx,sigma_yz,sigma_zz,sigma_zx,sigma_zy;
            sigma_xx=sigma_tensor(0,0);
            sigma_xy=sigma_tensor(0,1);
            sigma_xz=sigma_tensor(0,2);
            sigma_yx=sigma_tensor(1,0);
            sigma_yy=sigma_tensor(1,1);
            sigma_yz=sigma_tensor(1,2);
            sigma_zx=sigma_tensor(2,0);
            sigma_zy=sigma_tensor(2,1);
            sigma_zz=sigma_tensor(2,2);
           
            //Calculate the anisotropic conductivity of the background medium 
            // based on the principal axis resistivity of the medium and the Euler rotation angle.
            Eigen::Matrix3d sigma_b_primary;
            sigma_b_primary<<1/fmodel.rho_b[h][0], 0, 0, 0 ,1/fmodel.rho_b[h][1],0, 0 ,0 ,1/fmodel.rho_b[h][2];

            Eigen::Matrix3d RzfuAlpha_b_S,RzfuAlpha_b_D,RzfuAlpha_b_L,RzAlpha_b_L,RzAlpha_b_S,RzAlpha_b_D,sigma_b_tensor;
            RzfuAlpha_b_S<<cos(-fmodel.alpha_b_S[h]),sin(-fmodel.alpha_b_S[h]),0,-sin(-fmodel.alpha_b_S[h]),cos(-fmodel.alpha_b_S[h]),0,0,0,1;
            RzfuAlpha_b_D<<1,0,0,0,cos(-fmodel.alpha_b_D[h]),sin(-fmodel.alpha_b_D[h]),0,-sin(-fmodel.alpha_b_D[h]),cos(-fmodel.alpha_b_D[h]);
            RzfuAlpha_b_L<<cos(-fmodel.alpha_b_L[h]),sin(-fmodel.alpha_b_L[h]),0,-sin(-fmodel.alpha_b_L[h]),cos(-fmodel.alpha_b_L[h]),0,0,0,1;
            RzAlpha_b_L<<cos(fmodel.alpha_b_L[h]),sin(fmodel.alpha_b_L[h]),0,-sin(fmodel.alpha_b_L[h]),cos(fmodel.alpha_b_L[h]),0,0,0,1;
            RzAlpha_b_D<<1,0,0,0,cos(fmodel.alpha_b_D[h]),sin(fmodel.alpha_b_D[h]),0,-sin(fmodel.alpha_b_D[h]),cos(fmodel.alpha_b_D[h]);
            RzAlpha_b_S<<cos(fmodel.alpha_b_S[h]),sin(fmodel.alpha_b_S[h]),0,-sin(fmodel.alpha_b_S[h]),cos(fmodel.alpha_b_S[h]),0,0,0,1;
            sigma_b_tensor=RzfuAlpha_b_S*RzfuAlpha_b_D*RzfuAlpha_b_L*sigma_b_primary*RzAlpha_b_L*RzAlpha_b_D*RzAlpha_b_S;

            double sigma_b_xx,sigma_b_xy,sigma_b_xz,sigma_b_yy,sigma_b_yx,sigma_b_yz,sigma_b_zz,sigma_b_zx,sigma_b_zy;
            sigma_b_xx=sigma_b_tensor(0,0);
            sigma_b_xy=sigma_b_tensor(0,1);
            sigma_b_xz=sigma_b_tensor(0,2);
            sigma_b_yx=sigma_b_tensor(1,0);
            sigma_b_yy=sigma_b_tensor(1,1);
            sigma_b_yz=sigma_b_tensor(1,2);
            sigma_b_zx=sigma_b_tensor(2,0);
            sigma_b_zy=sigma_b_tensor(2,1);
            sigma_b_zz=sigma_b_tensor(2,2);

            //Calculate the anisotropic relative permeability based on the main shaft magnetization rate and Euler rotation angle
            Eigen::Matrix3d ms_primary;
            ms_primary<<fmodel.ms[h][0], 0, 0, 0 ,fmodel.ms[h][1],0, 0 ,0 ,fmodel.ms[h][2];

            Eigen::Matrix3d RzfuMs_S,RxfuMs_D,RxfuMs_L,RzMs_L,RzMs_S,RzMs_D,ms_tensor;
            RzfuMs_S<<cos(-fmodel.ms_S[h]),sin(-fmodel.ms_S[h]),0,-sin(-fmodel.ms_S[h]),cos(-fmodel.ms_S[h]),0,0,0,1;
            RxfuMs_D<<1,0,0,0,cos(-fmodel.ms_D[h]),sin(-fmodel.ms_D[h]),0,-sin(-fmodel.ms_D[h]),cos(-fmodel.ms_D[h]);
            RxfuMs_L<<cos(-fmodel.ms_L[h]),sin(-fmodel.ms_L[h]),0,-sin(-fmodel.ms_L[h]),cos(-fmodel.ms_L[h]),0,0,0,1;
            RzMs_L<<cos(fmodel.ms_L[h]),sin(fmodel.ms_L[h]),0,-sin(fmodel.ms_L[h]),cos(fmodel.ms_L[h]),0,0,0,1;
            RzMs_D<<1,0,0,0,cos(fmodel.ms_D[h]),sin(fmodel.ms_D[h]),0,-sin(fmodel.ms_D[h]),cos(fmodel.ms_D[h]);
            RzMs_S<<cos(fmodel.ms_S[h]),sin(fmodel.ms_S[h]),0,-sin(fmodel.ms_S[h]),cos(fmodel.ms_S[h]),0,0,0,1;
            ms_tensor=RzfuMs_S*RxfuMs_D*RxfuMs_L*ms_primary*RzMs_L*RzMs_D*RzMs_S;

            Eigen::Matrix3d mur_tensor,eye;
            eye.setIdentity();
            mur_tensor=eye+ms_tensor; 
            double mur_xx,mur_xy,mur_xz,mur_yy,mur_yx,mur_yz,mur_zz,mur_zx,mur_zy;
            mur_xx=mur_tensor(0,0);
            mur_xy=mur_tensor(0,1);
            mur_xz=mur_tensor(0,2);
            mur_yx=mur_tensor(1,0);
            mur_yy=mur_tensor(1,1);
            mur_yz=mur_tensor(1,2);
            mur_zx=mur_tensor(2,0);
            mur_zy=mur_tensor(2,1);
            mur_zz=mur_tensor(2,2);

            Eigen::Matrix3d mur_inv;
            mur_inv=mur_tensor.inverse();
            Eigen::Matrix3d gama;
            gama=eye-mur_inv;
            
    Eigen::Matrix4d Pe_xyPeT_yx,Pe_xzPeT_yx,Pe_xyPeT_zx,Pe_xzPeT_zx;
    Eigen::Matrix4d Pe_yxPeT_yx,Pe_yzPeT_yx,Pe_yxPeT_zx,Pe_yzPeT_zx;
    Eigen::Matrix4d Pe_xyPeT_xy,Pe_xzPeT_xy,Pe_xyPeT_zy,Pe_xzPeT_zy;
    Eigen::Matrix4d Pe_zxPeT_yx,Pe_zyPeT_yx,Pe_zxPeT_zx,Pe_zyPeT_zx;
    Eigen::Matrix4d Pe_yxPeT_xy,Pe_yzPeT_xy,Pe_yxPeT_zy,Pe_yzPeT_zy;
    Eigen::Matrix4d Pe_zxPeT_xy,Pe_zyPeT_xy,Pe_zxPeT_zy,Pe_zyPeT_zy;
    Eigen::Matrix4d Pe_xyPeT_xz,Pe_xzPeT_xz,Pe_xyPeT_yz,Pe_xzPeT_yz;
    Eigen::Matrix4d Pe_yxPeT_xz,Pe_yzPeT_xz,Pe_yxPeT_yz,Pe_yzPeT_yz;
    Eigen::Matrix4d Pe_zxPeT_xz, Pe_zyPeT_xz, Pe_zxPeT_yz, Pe_zyPeT_yz;
    Eigen::Matrix4d K2exx,K2exy,K2exz,K2eyx,K2eyy,K2eyz,K2ezx,K2ezy,K2ezz;

    Pe_xyPeT_yx << 2, 1, -2, -1,
              1, 2, -1, -2,
             -2, -1, 2, 1,
             -1, -2, 1, 2; 

   Pe_xzPeT_yx<< -1, -1,  1,  1,
                  1,  1, -1, -1,
                  -1, -1,  1,  1,
                   1,  1, -1, -1;   
        
    Pe_xyPeT_zx<< -1,  1, -1,  1,
                  -1,  1, -1,  1,
                  1, -1,  1, -1,
                  1, -1,  1, -1;   
          

    Pe_xzPeT_zx<<2, -2,  1, -1,
                   -2,  2, -1,  1,
                    1, -1,  2, -2,
                   -1,  1, -2,  2;

    
        //  K21
    Pe_yxPeT_yx<<-1, -1,  1,  1,
                 1,  1, -1, -1,
                 -1, -1,  1,  1,
                 1,  1, -1, -1;   
        
    Pe_yzPeT_yx<<1,  1, -1, -1,
                  1,  1, -1, -1,
                  -1, -1,  1,  1,
                  -1, -1,  1,  1;  
        
    Pe_yxPeT_zx<<1, -1,  1, -1,
                  -1,  1, -1,  1,
                 1, -1,  1, -1,
                  -1,  1, -1,  1;  
          
    Pe_yzPeT_zx<<-2,  2, -1,  1,
                  -1,  1, -2,  2,
                           2, -2,  1, -1,
                           1, -1,  2, -2; 
        //  K31
    
    Pe_zxPeT_yx<< 1,  1, -1, -1,
                           1,  1, -1, -1,
                          -1, -1,  1,  1,
                          -1, -1,  1,  1;          
        
    Pe_zyPeT_yx << -2, -1,  2,  1,
                           2,  1, -2, -1,
                          -1, -2,  1,  2,
                           1,  2, -1, -2;     
    Pe_zxPeT_zx<< -1,  1, -1,  1,
                              -1,  1, -1,  1,
                               1, -1,  1, -1,
                               1, -1,  1, -1;          
    Pe_zyPeT_zx<< 1, -1,  1, -1,
                          -1,  1, -1,  1,
                           1, -1,  1, -1,
                          -1,  1, -1,  1;
       //  K12
   
    Pe_xyPeT_xy<<-1,  1, -1,  1,
                              -1,  1, -1,  1,
                               1, -1,  1, -1,
                               1, -1,  1, -1;    
    Pe_xzPeT_xy<< 1, -1,  1, -1,
                          -1,  1, -1,  1,
                           1, -1,  1, -1,
                          -1,  1, -1,  1;       
    Pe_xyPeT_zy<< 1,  1, -1, -1,
                           1,  1, -1, -1,
                          -1, -1,  1,  1,
                          -1, -1,  1,  1;          
    Pe_xzPeT_zy<< -2, -1,  2,  1,
                           2,  1, -2, -1,
                          -1, -2,  1,  2,
                           1,  2, -1, -2;
        //  K22
    
    Pe_yxPeT_xy<<  2, -2,  1, -1,
                              -2,  2, -1,  1,
                               1, -1,  2, -2,
                              -1,  1, -2,  2;
    Pe_yzPeT_xy<< -1,  1, -1,  1,
                          -1,  1, -1,  1,
                           1, -1,  1, -1,
                           1, -1,  1, -1;
    Pe_yxPeT_zy<< -1, -1,  1,  1,
                           1,  1, -1, -1,
                          -1, -1,  1,  1,
                       1,  1, -1, -1;
    Pe_yzPeT_zy<<  2,  1, -2, -1,
                               1,  2, -1, -2,
                              -2, -1,  2,  1,
                              -1, -2,  1,  2;
       // K32
    
    Pe_zxPeT_xy<< -2,  2, -1,  1,
                          -1,  1, -2,  2,
                           2, -2,  1, -1,
                           1, -1,  2, -2;
    Pe_zyPeT_xy<<  1, -1,  1, -1,
                          -1,  1, -1,  1,
                           1, -1,  1, -1,
                          -1,  1, -1,  1;
    Pe_zxPeT_zy<<  1,  1, -1, -1,
                           1,  1, -1, -1,
                          -1, -1,  1,  1,
                          -1, -1,  1,  1;
    Pe_zyPeT_zy<< -1, -1,  1,  1,
                               1,  1, -1, -1,
                              -1, -1,  1,  1,
                               1,  1, -1, -1;
       //K13
    
    Pe_xyPeT_xz<<1,  1, -1, -1,
                           1,  1, -1, -1,
                          -1, -1,  1,  1,
                          -1, -1,  1,  1;
    Pe_xzPeT_xz<< -1, -1,  1,  1,
                               1,  1, -1, -1,
                              -1, -1,  1,  1,
                               1,  1, -1, -1;
    Pe_xyPeT_yz<< -2,  2, -1,  1,
                          -1,  1, -2,  2,
                           2, -2,  1, -1,
                           1, -1,  2, -2;
    Pe_xzPeT_yz<<  1, -1,  1, -1,
                          -1,  1, -1,  1,
                           1, -1,  1, -1,
                          -1,  1, -1,  1;
        
        //  K23
    
    Pe_yxPeT_xz<< -2, -1,  2,  1,
                           2,  1, -2, -1,
                          -1, -2,  1,  2,
                           1,  2, -1, -2;
    Pe_yzPeT_xz<<  1,  1, -1, -1,
                           1,  1, -1, -1,
                          -1, -1,  1,  1,
                          -1, -1,  1,  1;    
    Pe_yxPeT_yz<<  1, -1,  1, -1,
                          -1,  1, -1,  1,
                           1, -1,  1, -1,
                          -1,  1, -1,  1;           
    Pe_yzPeT_yz<< -1,  1, -1,  1,
                              -1,  1, -1,  1,
                               1, -1,  1, -1,
                               1, -1,  1, -1;  
    
        //  K33   
        
    Pe_zxPeT_xz<<  2,  1, -2, -1,
                               1,  2, -1, -2,
                              -2, -1,  2,  1,
                              -1, -2,  1,  2;        
    Pe_zyPeT_xz<< -1, -1,  1,  1,
                           1,  1, -1, -1,
                          -1, -1,  1,  1,
                           1,  1, -1, -1;
    Pe_zxPeT_yz<< -1,  1, -1,  1,
                          -1,  1, -1,  1,
                           1, -1,  1, -1,
                           1, -1,  1, -1;          
    Pe_zyPeT_yz<<  2, -2,  1, -1,
                              -2,  2, -1,  1,
                               1, -1,  2, -2,
                              -1,  1, -2,  2;  
    
    // K2e   
    K2exx<<4, 2, 2, 1,
               2, 4, 1, 2,
               2, 1, 4, 2,
               1, 2, 2, 4; 
    K2exy<<2, 1, 2, 1,
           2, 1, 2, 1,
           1, 2, 1, 2,
           1, 2, 1, 2;
    K2exz<<2,2,1,1,
           1,1,2,2,
           2,2, 1, 1,
           1, 1, 2, 2;
    K2eyx<<2,2, 1, 1,
           1, 1, 2, 2,
           2, 2, 1, 1,
           1, 1, 2, 2;
    K2eyy<<4, 2, 2, 1, 
               2, 4, 1, 2,
               2, 1, 4, 2,
               1, 2, 2, 4; 
    K2eyz<<2, 1, 2, 1,
           2, 1, 2, 1,
           1, 2, 1, 2,
           1, 2, 1, 2;
    K2ezx<<2, 1, 2, 1,
           2, 1, 2, 1,
           1, 2, 1, 2,
           1, 2, 1, 2;   
    K2ezy<<2, 2, 1, 1,
           1, 1, 2, 2,
           2, 2, 1, 1,
           1, 1, 2, 2;
    K2ezz<<4, 2, 2, 1,
               2, 4, 1, 2,
               2, 1, 4, 2,
               1, 2, 2, 4;
        int a_xn,b_yn,c_zn;
        a_xn=(h-1)%NX+1;
        b_yn=((h-a_xn)%(NX*NY))/NX+1;
        c_zn=floor((h-1)/(NX*NY)+1);
        
        double a,b,c;
        a=fmodel.A_X[a_xn-1]; 
        b=fmodel.B_Y[b_yn-1];
        c=fmodel.C_Z[c_zn-1];
        
        Pe_xyPeT_yx=a*b/6/c*Pe_xyPeT_yx;
        Pe_xzPeT_yx=a/4*Pe_xzPeT_yx;
        Pe_xyPeT_zx=a/4*Pe_xyPeT_zx;
        Pe_xzPeT_zx=a*c/6/b*Pe_xzPeT_zx;

        Pe_yxPeT_yx=a*b/4/c*Pe_yxPeT_yx;
        Pe_yzPeT_yx=b/4*Pe_yzPeT_yx;
        Pe_yxPeT_zx=a/4*Pe_yxPeT_zx; 
        Pe_yzPeT_zx=c/6*Pe_yzPeT_zx;

        Pe_zxPeT_yx=a/4*Pe_zxPeT_yx;
        Pe_zyPeT_yx=b/6*Pe_zyPeT_yx;
        Pe_zxPeT_zx=a*c/4/b*Pe_zxPeT_zx;
        Pe_zyPeT_zx=c/4*Pe_zyPeT_zx;

        Pe_xyPeT_xy=a*b/4/c*Pe_xyPeT_xy;
        Pe_xzPeT_xy=a/4*Pe_xzPeT_xy;
        Pe_xyPeT_zy=b/4*Pe_xyPeT_zy;
        Pe_xzPeT_zy=c/6*Pe_xzPeT_zy;

        Pe_yxPeT_xy=a*b/6/c*Pe_yxPeT_xy;
        Pe_yzPeT_xy=b/4*Pe_yzPeT_xy;
        Pe_yxPeT_zy=b/4*Pe_yxPeT_zy;
        Pe_yzPeT_zy=b*c/6/a*Pe_yzPeT_zy;

        Pe_zxPeT_xy=a/6*Pe_zxPeT_xy;
        Pe_zyPeT_xy=b/4*Pe_zyPeT_xy;
        Pe_zxPeT_zy=c/4*Pe_zxPeT_zy;
        Pe_zyPeT_zy=b*c/4/a*Pe_zyPeT_zy;

        Pe_xyPeT_xz=a/4*Pe_xyPeT_xz;
        Pe_xzPeT_xz=a*c/4/b*Pe_xzPeT_xz;
        Pe_xyPeT_yz=b/6*Pe_xyPeT_yz;
        Pe_xzPeT_yz=c/4*Pe_xzPeT_yz;

        Pe_yxPeT_xz=a/6*Pe_yxPeT_xz;
        Pe_yzPeT_xz=c/4*Pe_yzPeT_xz;
        Pe_yxPeT_yz=b/4*Pe_yxPeT_yz;
        Pe_yzPeT_yz=b*c/4/a*Pe_yzPeT_yz;

        Pe_zxPeT_xz=a*c/6/b*Pe_zxPeT_xz;
        Pe_zyPeT_xz=c/4*Pe_zyPeT_xz;
        Pe_zxPeT_yz=c/4*Pe_zxPeT_yz;
        Pe_zyPeT_yz=b*c/6/a*Pe_zyPeT_yz;

        K2exx=w*mu0*a*b*c/36*K2exx;
        K2exy=w*mu0*a*b*c/24*K2exy;
        K2exz=w*mu0*a*b*c/24*K2exz;
        K2eyx=w*mu0*a*b*c/24*K2eyx;
        K2eyy=w*mu0*a*b*c/36*K2eyy;
        K2eyz=w*mu0*a*b*c/24*K2eyz;
        K2ezx=w*mu0*a*b*c/24*K2ezx;
        K2ezy=w*mu0*a*b*c/24*K2ezy;
        K2ezz=w*mu0*a*b*c/36*K2ezz;


        Eigen::Matrix<complex<double>,4,4> K1exx,K1eyx,K1ezx,K1exy,K1eyy,K1ezy,K1exz,K1eyz,K1ezz;

        K1exx=mur_inv(1,1)*fmodel.Sey_origin[h]/fmodel.Sex_origin[h]/fmodel.Sez_origin[h]*Pe_xyPeT_yx
            +mur_inv(2,1)*Pe_xzPeT_yx+mur_inv(1,2)*Pe_xyPeT_zx
            +mur_inv(2,2)*fmodel.Sez_origin[h]/fmodel.Sex_origin[h]/fmodel.Sey_origin[h]*Pe_xzPeT_zx;

        K1eyx=mur_inv(0,1)*Pe_yxPeT_yx+mur_inv(2,1)*Pe_yzPeT_yx+mur_inv(0,2)*Pe_yxPeT_zx
            +mur_inv(2,2)*fmodel.Sez_origin[h]/fmodel.Sex_origin[h]/fmodel.Sey_origin[h]*Pe_yzPeT_zx;

        K1ezx=mur_inv(0,1)*Pe_zxPeT_yx+mur_inv(1,1)*fmodel.Sey_origin[h]/fmodel.Sex_origin[h]/fmodel.Sez_origin[h]*Pe_zyPeT_yx+
            mur_inv(0,2)*Pe_zxPeT_zx+mur_inv(1,2)*Pe_zyPeT_zx;
            
        K1exy=mur_inv(1,0)*Pe_xyPeT_xy+mur_inv(2,0)*Pe_xzPeT_xy+
            mur_inv(1,2)*Pe_xyPeT_zy+mur_inv(2,2)*fmodel.Sez_origin[h]/fmodel.Sex_origin[h]/fmodel.Sey_origin[h]*Pe_xzPeT_zy;

        K1eyy=mur_inv(0,0)*fmodel.Sex_origin[h]/fmodel.Sey_origin[h]/fmodel.Sez_origin[h]*Pe_yxPeT_xy+mur_inv(2,0)*Pe_yzPeT_xy+
            mur_inv(0,2)*Pe_yxPeT_zy+mur_inv(2,2)*fmodel.Sez_origin[h]/fmodel.Sex_origin[h]/fmodel.Sey_origin[h]*Pe_yzPeT_zy;

        K1ezy=mur_inv(0,0)*fmodel.Sex_origin[h]/fmodel.Sey_origin[h]/fmodel.Sez_origin[h]*Pe_zxPeT_xy+mur_inv(1,0)*Pe_zyPeT_xy+
            mur_inv(0,2)*Pe_zxPeT_zy+mur_inv(1,2)*Pe_zyPeT_zy;

        K1exz=mur_inv(1,0)*Pe_xyPeT_xz+mur_inv(2,0)*Pe_xzPeT_xz+
            mur_inv(1,1)*fmodel.Sey_origin[h]/fmodel.Sex_origin[h]/fmodel.Sez_origin[h]*Pe_xyPeT_yz+mur_inv(2,1)*Pe_xzPeT_yz;

        K1eyz=mur_inv(0,0)*fmodel.Sex_origin[h]/fmodel.Sey_origin[h]/fmodel.Sez_origin[h]*Pe_yxPeT_xz+mur_inv(2,0)*Pe_yzPeT_xz+
            mur_inv(0,1)*Pe_yxPeT_yz+mur_inv(2,1)*Pe_yzPeT_yz;

        K1ezz=mur_inv(0,0)*fmodel.Sex_origin[h]/fmodel.Sey_origin[h]/fmodel.Sez_origin[h]*Pe_zxPeT_xz+mur_inv(1,0)*Pe_zyPeT_xz+
            mur_inv(0,1)*Pe_zxPeT_yz+mur_inv(1,1)*fmodel.Sey_origin[h]/fmodel.Sex_origin[h]/fmodel.Sez_origin[h]*Pe_zyPeT_yz;
    
        Eigen::Matrix<complex<double>,12,12> K1e;
        K1e<<K1exx, K1exy, K1exz,
            K1eyx, K1eyy, K1eyz,
            K1ezx, K1ezy, K1ezz;

        Eigen::Matrix<complex<double>,4,4> K4exx,K4eyx,K4ezx,K4exy,K4eyy,K4ezy,K4exz,K4eyz,K4ezz;

        K4exx=gama(1,1)*fmodel.Sey_origin[h]/fmodel.Sex_origin[h]/fmodel.Sez_origin[h]*Pe_xyPeT_yx+gama(2,1)*Pe_xzPeT_yx+
              gama(1,2)*Pe_xyPeT_zx+gama(2,2)*fmodel.Sez_origin[h]/fmodel.Sex_origin[h]/fmodel.Sey_origin[h]*Pe_xzPeT_zx;
        K4eyx=gama(0,1)*Pe_yxPeT_yx+gama(2,1)*Pe_yzPeT_yx+
              gama(0,2)*Pe_yxPeT_zx+gama(2,2)*fmodel.Sez_origin[h]/fmodel.Sex_origin[h]/fmodel.Sey_origin[h]*Pe_yzPeT_zx;
        K4ezx=gama(0,1)*Pe_zxPeT_yx+gama(1,1)*fmodel.Sey_origin[h]/fmodel.Sex_origin[h]/fmodel.Sez_origin[h]*Pe_zyPeT_yx+
              gama(0,2)*Pe_zxPeT_zx+gama(1,2)*Pe_zyPeT_zx;
        K4exy=gama(1,0)*Pe_xyPeT_xy+gama(2,0)*Pe_xzPeT_xy+
              gama(1,2)*Pe_xyPeT_zy+gama(2,2)*fmodel.Sez_origin[h]/fmodel.Sex_origin[h]/fmodel.Sey_origin[h]*Pe_xzPeT_zy;
        K4eyy=gama(0,0)*fmodel.Sex_origin[h]/fmodel.Sey_origin[h]/fmodel.Sez_origin[h]*Pe_yxPeT_xy+gama(2,0)*Pe_yzPeT_xy+
              gama(0,2)*Pe_yxPeT_zy+gama(2,2)*fmodel.Sez_origin[h]/fmodel.Sex_origin[h]/fmodel.Sey_origin[h]*Pe_yzPeT_zy;
        K4ezy=gama(0,0)*fmodel.Sex_origin[h]/fmodel.Sey_origin[h]/fmodel.Sez_origin[h]*Pe_zxPeT_xy+gama(1,0)*Pe_zyPeT_xy+
              gama(0,2)*Pe_zxPeT_zy+gama(1,2)*Pe_zyPeT_zy;
        K4exz=gama(1,0)*Pe_xyPeT_xz+gama(2,0)*Pe_xzPeT_xz+
              gama(1,1)*fmodel.Sey_origin[h]/fmodel.Sex_origin[h]/fmodel.Sez_origin[h]*Pe_xyPeT_yz+gama(2,1)*Pe_xzPeT_yz;
        K4eyz=gama(0,0)*fmodel.Sex_origin[h]/fmodel.Sey_origin[h]/fmodel.Sez_origin[h]*Pe_yxPeT_xz+gama(2,0)*Pe_yzPeT_xz+
              gama(0,1)*Pe_yxPeT_yz+gama(2,1)*Pe_yzPeT_yz;
        K4ezz=gama(0,0)*fmodel.Sex_origin[h]/fmodel.Sey_origin[h]/fmodel.Sez_origin[h]*Pe_zxPeT_xz+gama(1,0)*Pe_zyPeT_xz+
              gama(0,1)*Pe_zxPeT_yz+gama(1,1)*fmodel.Sey_origin[h]/fmodel.Sex_origin[h]/fmodel.Sez_origin[h]*Pe_zyPeT_yz;

        Eigen::Matrix<complex<double>,12,12> K4e;
        K4e<<K4exx,K4exy,K4exz,
             K4eyx,K4eyy,K4eyz,
             K4ezx,K4ezy,K4ezz;
        complex<double> mysqrt(0,-1);

        Eigen::Matrix<complex<double>,12,12> K2e,K3e;
        K2e<<sigma_xx*mysqrt*K2exx*fmodel.Sez_origin[h]*fmodel.Sey_origin[h]/fmodel.Sex_origin[h],sigma_xy*mysqrt*K2exy,sigma_xz*mysqrt*K2exz,
             sigma_yx*mysqrt*K2eyx,sigma_yy*mysqrt*K2eyy*fmodel.Sex_origin[h]*fmodel.Sez_origin[h]/fmodel.Sey_origin[h],sigma_yz*mysqrt*K2eyz,
             sigma_zx*mysqrt*K2ezx,sigma_zy*mysqrt*K2ezy,sigma_zz*mysqrt*K2ezz*fmodel.Sex_origin[h]*fmodel.Sey_origin[h]/fmodel.Sez_origin[h];

        K3e<<(sigma_xx-sigma_b_xx)*mysqrt*K2exx*fmodel.Sez_origin[h]*fmodel.Sey_origin[h]/fmodel.Sex_origin[h],(sigma_xy-sigma_b_xy)*mysqrt*K2exy,(sigma_xz-sigma_b_xz)*mysqrt*K2exz,
             (sigma_yx-sigma_b_yx)*mysqrt*K2eyx,(sigma_yy-sigma_b_yy)*mysqrt*K2eyy*fmodel.Sex_origin[h]*fmodel.Sez_origin[h]/fmodel.Sey_origin[h],(sigma_yz-sigma_b_yz)*mysqrt*K2eyz,
             (sigma_zx-sigma_b_zx)*mysqrt*K2ezx,(sigma_zy-sigma_b_zy)*mysqrt*K2ezy,(sigma_zz-sigma_b_zz)*mysqrt*K2ezz*fmodel.Sex_origin[h]*fmodel.Sey_origin[h]/fmodel.Sez_origin[h];
        
        for(int j=1;j<=12;j++)
        { 
            for(int k=1;k<=12;k++)
            { 
                int NJ,NK;
                NJ=fmodel.ME[j][h]-1;
                NK=fmodel.ME[k][h]-1;
                //v0=K1+K2
                MatSetValue(v0, NJ,NK,K1e(j-1,k-1), ADD_VALUES);
                MatSetValue(v0, NJ,NK,K2e(j-1,k-1), ADD_VALUES);
                
                MatSetValue(K3, NJ,NK,-K3e(j-1,k-1), ADD_VALUES);
                MatSetValue(K4, NJ,NK,K4e(j-1,k-1), ADD_VALUES);
            } 
        }
    }
    //compress v0\k3\k4;
    MatAssemblyBegin(v0, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(v0, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(K3, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K3, MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(K4, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(K4, MAT_FINAL_ASSEMBLY);
    PetscTime(&end1);
    PetscPrintf(curComm,"_v0_time=%f\n",end1-start);
    
}
 
