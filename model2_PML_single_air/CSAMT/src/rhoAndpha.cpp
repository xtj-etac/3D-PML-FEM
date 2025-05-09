//Post-processing involves calculating the apparent resistivity and phase based on the formulas 51 and 52 from the paper.
#include "../include/main.h"

using namespace std;

rhoAndpha::rhoAndpha(int NX)
{
    Zxx.resize(NX); // Create NX sub-vectors
    for(auto& row : Zxx) {
        row.resize(NX); // Each sub-vector is initialized with NX elements.
    }

    Zxy.resize(NX); 
    for(auto& row : Zxy) {
        row.resize(NX); 
    }

    Zyx.resize(NX); 
    for(auto& row : Zyx) {
        row.resize(NX); 
    }

    Zyy.resize(NX); 
    for(auto& row : Zyy) {
        row.resize(NX); 
    }

    PhaseXX.resize(NX); 
    for(auto& row : PhaseXX) {
        row.resize(NX); 
    }

    PhaseXY.resize(NX); 
    for(auto& row :PhaseXY) {
        row.resize(NX); 
    }

    PhaseYX.resize(NX); 
    for(auto& row : PhaseYX) {
        row.resize(NX); 
    }

    PhaseYY.resize(NX); 
    for(auto& row : PhaseYY) {
        row.resize(NX); 
    }

    ResXX.resize(NX);
    for(auto& row : ResXX) {
        row.resize(NX); 
    }

    ResXY.resize(NX);
    for(auto& row :ResXY) {
        row.resize(NX);
    }

    ResYX.resize(NX);
    for(auto& row : ResYX) {
        row.resize(NX);
    }

    ResYY.resize(NX);
    for(auto& row : ResYY) {
        row.resize(NX);
    }
}

rhoAndpha::~rhoAndpha()
{
}
void rhoAndpha::rho(v1Assembly v1,v2Assembly v2,Fmodel fmodel,int group_id)
{
    int NX=fmodel.NX;
    int NY=fmodel.NX;
    double w=fmodel.omega;
    double mu=fmodel.miu0;
    complex<double> mysqrt(0,1);
    string s=to_string(group_id);
    cout<<"s="<<s<<endl;
    string f1="../outfile/"+s+"PhaseXY.txt";
    char *fPXY=(char *)f1.c_str();
    string f2="../outfile/"+s+"PhaseYX.txt";
    char *fPYX=(char *)f2.c_str();
    string f3="../outfile/"+s+"RseXY.txt";
    char *fRXY=(char *)f3.c_str();
    string f4="../outfile/"+s+"RseYX.txt";
    char *fRYX=(char *)f4.c_str();
    FILE* fpw_PhaseXY = fopen(fPXY,"w");
    FILE* fpw_PhaseYX = fopen(fPYX,"w");
    FILE* fpw_RseXY = fopen(fRXY,"w");
    FILE* fpw_RseYX = fopen(fRYX,"w");
    cout<<"start for循环"<<endl;
    for(int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {
            
            complex<double> Temp=v1.Hxx1[i][j]*v2.Hyy2[i][j]-v1.Hyy1[i][j]*v2.Hxx2[i][j];
            Zxx[i][j]=(v1.Exx1[i][j]*v2.Hyy2[i][j]-v2.Exx2[i][j]*v1.Hyy1[i][j])/Temp;
            Zxy[i][j]=(v2.Exx2[i][j]*v1.Hxx1[i][j]-v1.Exx1[i][j]*v2.Hxx2[i][j])/Temp;
            Zyx[i][j]=(v1.Eyy1[i][j]*v2.Hyy2[i][j]-v2.Eyy2[i][j]*v1.Hyy1[i][j])/Temp;
            Zyy[i][j]=(v2.Eyy2[i][j]*v1.Hxx1[i][j]-v1.Eyy1[i][j]*v2.Hxx2[i][j])/Temp;
            PhaseXX[i][j]=-atan((Zxx[i][j].imag())/(Zxx[i][j].real()))*180/M_PI;  
            PhaseXY[i][j]=-atan((Zxy[i][j].imag())/(Zxy[i][j].real()))*180/M_PI;
            fprintf(fpw_PhaseXY,"%e\n",PhaseXY[i][j]);
            PhaseYX[i][j]=-atan((Zyx[i][j].imag())/(Zyx[i][j].real()))*180/M_PI;
            fprintf(fpw_PhaseYX,"%e \n",PhaseYX[i][j]);
            PhaseYY[i][j]=-atan((Zyy[i][j].imag())/(Zyy[i][j].real()))*180/M_PI;
            ResXY[i][j]=abs(((Zxy[i][j])*(Zxy[i][j]))*mysqrt/(w*mu));
            fprintf(fpw_RseXY,"%e \n",ResXY[i][j]);
            ResYX[i][j]=abs(((Zyx[i][j])*(Zyx[i][j]))*mysqrt/(w*mu));
            fprintf(fpw_RseYX,"%e \n",ResYX[i][j]);
            ResXX[i][j]=abs((Zxx[i][j])*(Zxx[i][j])*mysqrt/(w*mu));
            ResYY[i][j]=abs((Zyy[i][j])*(Zyy[i][j])*mysqrt/(w*mu));        
        }
        fprintf(fpw_PhaseXY,"\n");
        fprintf(fpw_PhaseYX,"\n");
        fprintf(fpw_RseXY,"\n");
        fprintf(fpw_RseYX,"\n");
    }
    int midx = floor(NX/2)-1, midy = floor(NY/2)-1;
    printf("PhaseXX mid: %f\n", PhaseXX[midx][midy]);
    printf("PhaseXY mid: %f\n", PhaseXY[midx][midy]);
    printf("PhaseYX mid: %f\n", PhaseYX[midx][midy]);
    printf("PhaseYY mid: %f\n", PhaseYY[midx][midy]);

    printf("RseXX mid: %f\n", ResXX[midx][midy]);
    printf("RseXY mid: %f\n", ResXY[midx][midy]);
    printf("RseYX mid: %f\n", ResYX[midx][midy]);
    printf("RseYY mid: %f\n", ResYY[midx][midy]);

}
