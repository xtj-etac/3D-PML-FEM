/*
 * @Author: Shutian Shutian@653.com
 * @Date: 2023-11-14 19:07:09
 * @LastEditors: Shutian Shutian@653.com
 * @LastEditTime: 2023-12-13 17:32:48
 * @FilePath: /st/project2/matlabC_v3/HYMAX-APP/codeMAX/rhoAndpha.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "../include/main.h"
//#include "../include/fmodel.h"

using namespace std;

rhoAndpha::rhoAndpha(int NX)
{
    Zxx.resize(NX); // 创建NX个子vector
    for(auto& row : Zxx) {
        row.resize(NX); // 每个子vector初始化为NX个元素
    }

    Zxy.resize(NX); // 创建NX个子vector
    for(auto& row : Zxy) {
        row.resize(NX); // 每个子vector初始化为NX个元素
    }

    Zyx.resize(NX); // 创建NX个子vector
    for(auto& row : Zyx) {
        row.resize(NX); // 每个子vector初始化为NX个元素
    }

    Zyy.resize(NX); // 创建NX个子vector
    for(auto& row : Zyy) {
        row.resize(NX); // 每个子vector初始化为NX个元素
    }

    PhaseXX.resize(NX); // 创建NX个子vector
    for(auto& row : PhaseXX) {
        row.resize(NX); // 每个子vector初始化为NX个元素
    }

    PhaseXY.resize(NX); // 创建NX个子vector
    for(auto& row :PhaseXY) {
        row.resize(NX); // 每个子vector初始化为NX个元素
    }

    PhaseYX.resize(NX); // 创建NX个子vector
    for(auto& row : PhaseYX) {
        row.resize(NX); // 每个子vector初始化为NX个元素
    }

    PhaseYY.resize(NX); // 创建NX个子vector
    for(auto& row : PhaseYY) {
        row.resize(NX); // 每个子vector初始化为NX个元素
    }

    ResXX.resize(NX); // 创建NX个子vector
    for(auto& row : ResXX) {
        row.resize(NX); // 每个子vector初始化为NX个元素
    }

    ResXY.resize(NX); // 创建NX个子vector
    for(auto& row :ResXY) {
        row.resize(NX); // 每个子vector初始化为NX个元素
    }

    ResYX.resize(NX); // 创建NX个子vector
    for(auto& row : ResYX) {
        row.resize(NX); // 每个子vector初始化为NX个元素
    }

    ResYY.resize(NX); // 创建NX个子vector
    for(auto& row : ResYY) {
        row.resize(NX); // 每个子vector初始化为NX个元素
    }
}

rhoAndpha::~rhoAndpha()
{
}
void rhoAndpha::rho(v1Assembly v1,v2Assembly v2,Fmodel fmodel,int group_id)
{
    int NX=fmodel.NX;
    int NY=fmodel.NY;
    double w=fmodel.omega;
    double mu=fmodel.miu0;
    cout<<"start Zxx"<<endl;
    complex<double> mysqrt(0,1);
    string s=to_string(group_id);//0/1
    string f1="../outfile/"+s+"PhaseXY.txt";
    char *fPXY=(char *)f1.c_str();
    string f2="../outfile/"+s+"PhaseYX.txt";
    char *fPYX=(char *)f2.c_str();
    string f3="../outfile/"+s+"RseXY.txt";
    char *fRXY=(char *)f3.c_str();
    string f4="../outfile/"+s+"RseYX.txt";
    char *fRYX=(char *)f4.c_str();
   // FILE* fpw_PhaseXX = fopen("../outfile/PhaseXX.txt","w");
    FILE* fpw_PhaseXY = fopen(fPXY,"w");
    FILE* fpw_PhaseYX = fopen(fPYX,"w");
   // FILE* fpw_PhaseYY = fopen("../outfile/PhaseYY.txt","w");
    FILE* fpw_RseXY = fopen(fRXY,"w");
    FILE* fpw_RseYX = fopen(fRYX,"w");
    // cout<<"start 2"<<endl;
    for(int i=0;i<NX;i++)
    {
        for(int j=0;j<NY;j++)
        {
            
            complex<double> Temp=v1.Hxx1[i][j]*v2.Hyy2[i][j]-v1.Hyy1[i][j]*v2.Hxx2[i][j];
            Zxx[i][j]=(v1.Exx1[i][j]*v2.Hyy2[i][j]-v2.Exx2[i][j]*v1.Hyy1[i][j])/Temp;
            Zxy[i][j]=(v2.Exx2[i][j]*v1.Hxx1[i][j]-v1.Exx1[i][j]*v2.Hxx2[i][j])/Temp;
            Zyx[i][j]=(v1.Eyy1[i][j]*v2.Hyy2[i][j]-v2.Eyy2[i][j]*v1.Hyy1[i][j])/Temp;
            Zyy[i][j]=(v2.Eyy2[i][j]*v1.Hxx1[i][j]-v1.Eyy1[i][j]*v2.Hxx2[i][j])/Temp;
            //printf("Zxx[i][j]=%e,%e\n",Zxx[i][j].real(),Zxx[i][j].imag());
            PhaseXX[i][j]=-atan((Zxx[i][j].imag())/(Zxx[i][j].real()))*180/M_PI;  
           // fprintf(fpw_PhaseXX,"%e\n",i,j,PhaseXX[i][j]);
            //printf("PhaseXX=%f\n",PhaseXX[i][j]);
            PhaseXY[i][j]=-atan((Zxy[i][j].imag())/(Zxy[i][j].real()))*180/M_PI;
            
            fprintf(fpw_PhaseXY,"%e\n",PhaseXY[i][j]);

            PhaseYX[i][j]=-atan((Zyx[i][j].imag())/(Zyx[i][j].real()))*180/M_PI;
            fprintf(fpw_PhaseYX,"%e \n",PhaseYX[i][j]);
            PhaseYY[i][j]=-atan((Zyy[i][j].imag())/(Zyy[i][j].real()))*180/M_PI;
            //fprintf(fpw_PhaseYY,"(%d,%d)=%e\n",i,j,PhaseYY[i][j]);

            ResXY[i][j]=abs(((Zxy[i][j])*(Zxy[i][j]))*mysqrt/(w*mu));
            fprintf(fpw_RseXY,"%e \n",ResXY[i][j]);
            ResYX[i][j]=abs(((Zyx[i][j])*(Zyx[i][j]))*mysqrt/(w*mu));
            fprintf(fpw_RseYX,"%e \n",ResYX[i][j]);
            ResXX[i][j]=abs((Zxx[i][j])*(Zxx[i][j])*mysqrt/(w*mu));
            //fprintf(fpw_RseXX,"%e \n",ResXX[i][j]);
            ResYY[i][j]=abs((Zyy[i][j])*(Zyy[i][j])*mysqrt/(w*mu));
            //fprintf(fpw_RseYY,"%e \n",ResYY[i][j]);

            
        }
        fprintf(fpw_PhaseXY,"\n");
        fprintf(fpw_PhaseYX,"\n");
        fprintf(fpw_RseXY,"\n");
        fprintf(fpw_RseYX,"\n");
        //printf("第i=%d行\n",i);
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
