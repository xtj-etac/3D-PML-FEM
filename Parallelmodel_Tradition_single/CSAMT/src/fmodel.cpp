#include "../include/main.h"
using namespace std;

Fmodel::Fmodel(double f)
{
    epsilo0=1.0/(36*M_PI)*1e-9;
    eps0=1/(36*M_PI)*1e-9;
    miu0=M_PI*4e-7;
    omega=2*M_PI*f;
    NX=N_MODEL_X+2*Grid_extend_X;
    NY=N_MODEL_Y+2*Grid_extend_Y;
    NZ=N_MODEL_Z+2*Grid_extend_Z;
    NE=NX*NY*NZ; //模型的六面体单元总数（单元总数）
    NL=NX*(NY+1)*(NZ+1)+(NX+1)*NY*(NZ+1)+(NX+1)*(NY+1)*NZ;
    NN=0;
    m=1;
    sigma_max = find_sigma * miu0 * sqrt(miu0/eps0)/PML_DZ;
    kappa_max = find_kappa;
    alpha_max = find_alpha;
    //printf("epsilo0=%e ,eps0=%e ,miu0=%e, omega=%e\n",epsilo0,eps0,miu0,omega);
    //printf("NE=%d, NL=%d\n",NE,NL);
    //printf("sigma_max=%e ,kappa_max=%e ,alpha_max=%e\n",sigma_max,kappa_max,alpha_max);
    
    init_Fmodel();

    
}

void Fmodel::init_Fmodel()
{

    //int lenABC=46;
    int lenABC=N_MODEL_X+2*Grid_extend_X;
    A_X=new double[lenABC];
    B_Y=new double[lenABC];
    C_Z=new double[lenABC];
    
    double extend_DX[Grid_extend_X]={125,156,195,244,305,381,476,596,745,931,1160,1450,1820,2270,2840,3550,4440,5550,6940,8670,10800,13600,16900};

    for(int i=0;i<lenABC;i++)
    {
        if(i<Grid_extend_X)
        {
            A_X[i]=extend_DX[Grid_extend_X-i-1];
            B_Y[i]=extend_DX[Grid_extend_X-i-1];
            C_Z[i]=extend_DX[Grid_extend_X-i-1];
        }else if(i>=Grid_extend_X&&i<Grid_extend_X+N_MODEL_X)
        {
            A_X[i]=100;
            B_Y[i]=100;
            C_Z[i]=100;
        }else
        {
            A_X[i]=extend_DX[i-(Grid_extend_X+N_MODEL_X)];
            B_Y[i]=extend_DX[i-(Grid_extend_X+N_MODEL_X)];
            C_Z[i]=extend_DX[i-(Grid_extend_X+N_MODEL_X)];
        }
            
    }
   

 /*   cout << "Values A_X:" << endl;
    for (auto it = 0; it <lenABC; ++it) {
        std::cout << A_X[it] << " ";
    }
    cout << endl;

    cout << "Values B_Y:" << endl;
    for (auto it = 0; it <lenABC; ++it) {
        std::cout << B_Y[it] << " ";
    }
    cout << endl;

    cout << "Values C_Z:" << endl;
    for (auto it = 0; it <lenABC; ++it) {
        std::cout << C_Z[it] << " ";
    }
    cout << endl;*/
}

Fmodel::~Fmodel(){
   /* delete []A_X;
    delete []B_Y;
    delete []C_Z;
    delete []xekappa;
    delete []xesigma;
    delete []xealpha;
    delete []yekappa;
    delete []yesigma;
    delete []yealpha;
    delete []zekappa;
    delete []zesigma;
    delete []zealpha;
    delete []Sex_origin;
    delete []Sey_origin;
    delete []Sez_origin;*/

    //freeFmodel();
}

//设置h、rho、alpha_S、alpha_D,alpha_L的初始值
void Fmodel::settMatrix()
{
    int m_size=NE+1;
    ME.resize(13, std::vector<int>(m_size, 0.0));
    rho.resize(m_size, array<double, 3>{0.0, 0.0, 0.0});
    rho_b.resize(m_size, array<double, 3>{0.0, 0.0, 0.0});
    ms.resize(m_size, array<double, 3>{0.0, 0.0, 0.0});

    ms_S.resize(m_size,0);
    ms_D.resize(m_size,0);
    ms_L.resize(m_size,0);

    alpha_S.resize(m_size,0);
    alpha_D.resize(m_size,0);
    alpha_L.resize(m_size,0);

    alpha_b_S.resize(m_size,0);
    alpha_b_D.resize(m_size,0);
    alpha_b_L.resize(m_size,0);

    int h=(Nair+NN)*NX*NY;
    int h2=NZ*NX*NY;
    //cout<<"h="<<h<<endl;

   // #pragma omp parallel for
    for(int i=1;i<=h;i++)
    {
        for(int j=0;j<3;j++)
        {
            rho[i][j]=1e12;
            rho_b[i][j]=1e12;
        }
    }      
   // #pragma omp parallel for
    for(int i=h+1;i<=h2;i++)//h-h2 
    {
        for(int j=0;j<3;j++)
        {
            rho[i][j]=100.0;
            rho_b[i][j]=100.0;
        }
    }    
 
    int tmp=0;
    //#pragma omp parallel for num_threads(4) lastprivate(tmp)
    for(int L=NZ/2-1;L<=NZ/2+2;L++)
    {
        for(int M=NY/2-1;M<=NY/2+2;M++)
        {
            for(int N=NX/2-3;N<=NX/2+4;N++)
            {
                tmp=(L-1)*NX*NY+(M-1)*NX+N;
                rho[tmp][0]=70.0;
                rho[tmp][1]=30.0;
                rho[tmp][2]=50.0;
                alpha_S[tmp]=45;
                alpha_D[tmp]=30;
                alpha_L[tmp]=10;
                ms[tmp][0]=0.3;
                ms[tmp][1]=0.2;
                ms[tmp][2]=0.1;
                ms_S[tmp]=40;
                ms_D[tmp]=20;
                ms_L[tmp]=30;
            }
        }
    } 

    xekappa = new double[m_size];
    xesigma = new double[m_size];
    xealpha = new double[m_size];
    yekappa = new double[m_size];
    yesigma = new double[m_size];
    yealpha = new double[m_size];
    zekappa = new double[m_size];
    zesigma = new double[m_size];
    zealpha = new double[m_size];
    Sex_origin = new complex<double>[m_size];
    Sey_origin = new complex<double>[m_size];
    Sez_origin = new complex<double>[m_size];

    //初始化某些值为1：
    for(int i=0;i<m_size;i++)
    {
        xekappa[i]=1.0;
        yekappa[i]=1.0;
        zekappa[i]=1.0;
        Sex_origin[i]=complex<double>(1.0,0);
        Sey_origin[i]=complex<double>(1.0,0);
        Sez_origin[i]=complex<double>(1.0,0);
    }
    
 /*   for(int k=1;k<=NZ;k++)
    {
        for(int j=1;j<=NY;j++)
        {
            for(int i=1;i<=NPML_X;i++)
            {
                int h=(k-1)*NX*NY+(j-1)*NX+i;
                double distance=(NPML_X+0.5-i)/NPML_X;
                xesigma[h]=sigma_max*(exp(m*distance)-1)/sqrt(omega*eps0/100);
                xekappa[h]=1+kappa_max*(exp(m*distance)-1);
                xealpha[h]=alpha_max*(exp(m*(1-distance))-1);
                Sex_origin[h]=xekappa[h]+(sqrt(2)*xesigma[h]/(xealpha[h]-complex<double>(0, 1)));
            }

            for(int i=NPML_X+N_MODEL_X+1;i<=NX;i++)
            {
                int h=(k-1)*NX*NY+(j-1)*NX+i;
                double distance=(i-0.5-(NPML_X+N_MODEL_X))/NPML_X;
                xesigma[h]=sigma_max*(exp(m*distance)-1)/sqrt(omega*eps0/100);
                xekappa[h]=1+kappa_max*(exp(m*distance)-1);
                xealpha[h]=alpha_max*(exp(m*(1-distance))-1);
                Sex_origin[h]=xekappa[h]+(sqrt(2)*xesigma[h]/(xealpha[h]-complex<double>(0, 1)));
            }
            
        }
    }

    for(int k=1;k<=NZ;k++)
    {
        for(int i=1;i<=NX;i++)
        {
            for(int j=1;j<=NPML_Y;j++)
            {
                int h=(k-1)*NX*NY+(j-1)*NX+i;
                double distance=(NPML_Y+0.5-j)/NPML_Y;
                yesigma[h]=sigma_max*(exp(m*distance)-1)/sqrt(omega*eps0/100);
                yekappa[h]=1+kappa_max*(exp(m*distance)-1);
                yealpha[h]=alpha_max*(exp(m*(1-distance))-1);
                Sey_origin[h]=yekappa[h]+(sqrt(2)*yesigma[h]/(yealpha[h]-complex<double>(0, 1)));
            }
   
            for(int j=NPML_Y+N_MODEL_Y+1;j<=NY;j++)
            {
                int h=(k-1)*NX*NY+(j-1)*NX+i;
                double distance=(j-0.5-(NPML_Y+N_MODEL_Y))/NPML_Y;
                yesigma[h]=sigma_max*(exp(m*distance)-1)/sqrt(omega*eps0/100);
                yekappa[h]=1+kappa_max*(exp(m*distance)-1);
                yealpha[h]=alpha_max*(exp(m*(1-distance))-1);
                Sey_origin[h]=yekappa[h]+(sqrt(2)*yesigma[h]/(yealpha[h]-complex<double>(0, 1)));
            }
            
        }
    }


    for(int j=1;j<=NY;j++)
    {
        for(int i=1;i<=NX;i++)
        {
            for(int k=1;k<=NPML_Z;k++)
            {
                int h=(k-1)*NX*NY+(j-1)*NX+i;
                double distance=(NPML_Z+0.5-k)/NPML_Z;
                zesigma[h]=sigma_max*(exp(m*distance)-1)/sqrt(omega*eps0/100);
                zekappa[h]=1+kappa_max*(exp(m*distance)-1);
                zealpha[h]=alpha_max*(exp(m*(1-distance))-1);
                Sez_origin[h]=zekappa[h]+(sqrt(2)*zesigma[h]/(zealpha[h]-complex<double>(0, 1)));
                
            }

            for(int k=NPML_Z+N_MODEL_Z+1;k<=NZ;k++)
            {
                int h=(k-1)*NX*NY+(j-1)*NX+i;
                double distance=(k-0.5-(NPML_Z+N_MODEL_Z))/NPML_Y;
                zesigma[h]=sigma_max*(exp(m*distance)-1)/sqrt(omega*eps0/100);
                zekappa[h]=1+kappa_max*(exp(m*distance)-1);
                zealpha[h]=alpha_max*(exp(m*(1-distance))-1);
                Sez_origin[h]=zekappa[h]+(sqrt(2)*zesigma[h]/(zealpha[h]-complex<double>(0, 1)));
            }
        }
    }
*/


        //omp_set_max_active_levels(1);//设置支持嵌套并行
       // #pragma omp parallel for num_threads(10)
        for(int L=1;L<=NZ;L++)
        {
          //  #pragma omp parallel for num_threads(10)// lastprivate(ME)
            for(int M=1;M<=NY;M++)
            {
                for(int N=1;N<=NX;N++)
                {
                    long tmp=(L-1)*NY*NX+(M-1)*NX+N;
                    long IL1=(L-1)*(NX*(NY+1)+NY*(NX+1)+(NX+1)*(NY+1));
                    ME[1][tmp]=IL1+(M-1)*(2*NX+1)+N;
                    ME[2][tmp]=IL1+(M-1)*(2*NX+1)+N+2*NX+1;
                    ME[3][tmp]=IL1+(M-1)*(2*NX+1)+N+NX*(NY+1)+NY*(NX+1)+(NX+1)*(NY+1);
                    ME[4][tmp]=IL1+(M-1)*(2*NX+1)+N+NX*(NY+1)+NY*(NX+1)+(NX+1)*(NY+1)+2*NX+1;
                    ME[5][tmp]=IL1+(M-1)*(2*NX+1)+N+NX;
                    ME[6][tmp]=IL1+(M-1)*(2*NX+1)+N+NX+NX*(NY+1)+NY*(NX+1)+(NX+1)*(NY+1);
                    ME[7][tmp]=IL1+(M-1)*(2*NX+1)+N+NX+1;
                    ME[8][tmp]=IL1+(M-1)*(2*NX+1)+N+NX+1+NX*(NY+1)+NY*(NX+1)+(NX+1)*(NY+1);
                    ME[9][tmp]=IL1+NX*(NY+1)+NY*(NX+1)+(M-1)*(NX+1)+N;
                    ME[10][tmp]=IL1+NX*(NY+1)+NY*(NX+1)+(M-1)*(NX+1)+N+1;
                    ME[11][tmp]=IL1+NX*(NY+1)+NY*(NX+1)+(M-1)*(NX+1)+N+1+NX;
                    ME[12][tmp]=IL1+NX*(NY+1)+NY*(NX+1)+(M-1)*(NX+1)+N+2+NX;
                }
            }
        }  
    
}

void Fmodel::freeFmodel()
{
    ME.clear();
    ME.shrink_to_fit();
    rho.clear();
    rho.shrink_to_fit();
    rho_b.clear();
    ms.clear();
    ms.shrink_to_fit();
    alpha_b_D.clear();
    alpha_b_D.shrink_to_fit();
    alpha_D.clear();
    alpha_D.shrink_to_fit();
    alpha_b_L.clear();
    alpha_b_L.shrink_to_fit();
    alpha_L.clear();
    alpha_L.shrink_to_fit();
    alpha_b_S.clear();
    alpha_b_S.shrink_to_fit();
    alpha_S.clear();
    alpha_S.shrink_to_fit();
    ms_D.clear();
    ms_D.shrink_to_fit();
    ms_L.clear();
    ms_L.shrink_to_fit();
    ms_S.clear();
    ms_S.shrink_to_fit();
}
/*void Fmodel::printME(int (*ME)[34201])
{
        FILE *outFile=fopen("../outfile/ME.txt","w");
        for(int i=1;i<9;i++)
        {
            for(int j=1;j<34201;j++)
            {
                fprintf(outFile,"%d \n",ME[i][j]);
                fflush(outFile);
                //if(i==2&&j==)
            }   
            //fprintf(outFile,"\n");
        }
        //printf("%ld\n",ME[4][25]);
}

void Fmodel::printRho(double (*rho)[3],double *alpha_S,Config config)
{ 
            FILE *outFile3=fopen("../outfile/alpha_S.txt","w");
            for(int i=1;i<config.h2+1;i++)
            { 
                if(alpha_S[i]!=0.0)
                {
                    fprintf(outFile3,"i=%d %f\n",i,alpha_S[i]);
                    fflush(outFile3);
                }
            }
            FILE *outFile2=fopen("../outfile/rho2.txt","w");
            for(int i=1;i<config.h2+1;i++)
            {
                    fprintf(outFile2,"%f  %f  %f\n",rho[i][0],rho[i][1],rho[i][2]);
                    fflush(outFile2);
            }
} */
