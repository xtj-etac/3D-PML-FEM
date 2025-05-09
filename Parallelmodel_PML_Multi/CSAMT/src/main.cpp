#include "../include/main.h"
using namespace std;
// 获取进程的内存使用情况
long getMemoryUsage() {
    std::ifstream file("/proc/self/status");
    std::string line;
    long vmsize = 0, rss = 0;

    while (std::getline(file, line)) {
        if (line.find("VmSize:") != std::string::npos) {
            sscanf(line.c_str(), "VmSize: %ld kB", &vmsize);
        } else if (line.find("VmRSS:") != std::string::npos) {
            sscanf(line.c_str(), "VmRSS: %ld kB", &rss);
        }
    }

    return rss * 1024; // 返回 RSS 内存使用量（字节）
}

int main(int argc,char **argv)
{ 
    // 打印初始内存使用情况
    //std::cout << "Initial Memory Usage: " << getMemoryUsage() / 1024.0 << " KB" << std::endl;
    PetscErrorCode ierr;
    PetscMPIInt rank,size;
    //MPI_Init(&argc,&argv);
    PetscFunctionBeginUser;
    ierr = PetscInitialize(&argc, &argv, NULL, "Program help message or NULL");
    if (ierr) {
        PetscPrintf(PETSC_COMM_WORLD, "Error during initialization: %d\n", ierr);
        return ierr;
    }

    PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
    PetscPrintf(PETSC_COMM_WORLD, "Initial Memory Usage:%f KB\n",getMemoryUsage() / 1024.0);
    clock_t start,finish;
    //double f[16]={0.01,0.1,0.2,0.5,1,2,5,10,20,50,100,100,200,500,1000,2000};
    double f[8]={0.125,0.25,0.375,0.5,0.625,0.75,0.875,1};
    //double f[2]={10,10};//
    PetscLogDouble starttime,endtime,star1,end1,Alltime;
    int f_size=sizeof(f)/sizeof(f[0]);//
    MPI_Status status;
    //
    PetscMPIInt group_id;//
    group_id=rank%f_size;//color=0,1
    MPI_Comm row_comm;
    PetscCallMPI(MPI_Comm_split(PETSC_COMM_WORLD,group_id,rank,&row_comm));//{0,1,2,3}{4,5,6,7}
    //
    PetscMPIInt row_rank,row_size;
    PetscCallMPI(MPI_Comm_rank(row_comm,&row_rank));//
    PetscCallMPI(MPI_Comm_size(row_comm,&row_size));
    PetscTime(&star1);
    int FreqStep=0;
    if(size==1)//进程数为0时，Freq进位为1;
    {
        FreqStep=1;
    }else{//否则按照频点个数为步长跳转
        FreqStep=f_size;
    }
    for(int Freq=group_id;Freq<f_size;Freq+=FreqStep)
    {
        PetscPrintf(row_comm,"Freq=%d  rank=%d f[Freq]=%f\n",Freq,row_rank,f[Freq]);
        PetscPrintf(row_comm,"start fmodel\n");
        double curFreq=f[Freq];
        Fmodel fmodel(curFreq); 
        fmodel.settMatrix();

        //2.v0Assembly.m
        PetscPrintf(row_comm,"开始计算v0\n");
        v0Assembly v0A(curFreq);
        v0A.creatv0(fmodel,row_comm);
       //3.v1Assembly.m
        v1Assembly v1(fmodel,curFreq);
        v1.fDirBdaries(Freq);
        v1.creatv1(v0A,fmodel,row_comm);
      //4.v2Assembly.m
        v2Assembly v2(fmodel,curFreq);
        v2.fDirBdaries(Freq);
        v2.creatv2(v0A,fmodel,row_comm);
        //5.rhoAndpha.m
        if(row_rank==0)
        {
            cout<<"start startrho fmodel.NX="<<fmodel.NX<<endl;
            double startrho=MPI_Wtime();
            rhoAndpha rho_alpha(fmodel.NX);
            rho_alpha.rho(v1,v2,fmodel,Freq);
            double endrho=MPI_Wtime();
        } 
    }
    // 打印分配内存后的内存使用情况
   
    PetscCallMPI(MPI_Comm_free(&row_comm));//
    PetscTime(&end1);//=MPI_Wtime();  
    Alltime=end1-star1;
    MPI_Barrier(PETSC_COMM_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "===fmodel time(MPI):%f seconds\n",Alltime);

    ierr = PetscFinalize();
    if (ierr) {
        PetscPrintf(PETSC_COMM_WORLD, "Error during finalization: %d\n", ierr);
        return ierr;
    }
    std::cout << "Memory Usage after allocation: " << getMemoryUsage() / 1024.0 << " KB" << std::endl;
    return 0;
}
