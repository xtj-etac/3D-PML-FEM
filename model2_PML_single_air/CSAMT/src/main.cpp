#include "../include/main.h"
using namespace std;
// Obtain the memory usage of the process
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

    return rss * 1024; // Return RSS memory usage (bytes)
}

int main(int argc,char **argv)
{ 
    // Print the initial memory usage
    PetscErrorCode ierr;
    PetscMPIInt rank,size;
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
    double f[1]={0.125};//test frequency
    PetscLogDouble starttime,endtime,start1,end1,Alltime;
    int f_size=sizeof(f)/sizeof(f[0]);
    MPI_Status status;
    PetscMPIInt group_id;
    group_id=rank%f_size;
    MPI_Comm row_comm;
    //Dividing the communication domain
    PetscCallMPI(MPI_Comm_split(PETSC_COMM_WORLD,group_id,rank,&row_comm));
    PetscMPIInt row_rank,row_size;
    PetscCallMPI(MPI_Comm_rank(row_comm,&row_rank));
    PetscCallMPI(MPI_Comm_size(row_comm,&row_size));
    PetscTime(&start1);
    int FreqStep=0;
    if(size==1)//When the process count is 1, the Freq rollover occurs with a value of 1.
    {
        FreqStep=1;
    }else{//Otherwise, the jump will be made based on the number of frequency points as the step size.
        FreqStep=f_size;
    }
    for(int Freq=group_id;Freq<f_size;Freq+=FreqStep)
    {
        PetscPrintf(row_comm,"Freq=%d  rank=%d f[Freq]=%f\n",Freq,row_rank,f[Freq]);
        PetscPrintf(row_comm,"start fmodel\n");
        PetscLogDouble startfmodel,endtfmodel ;
        PetscTime(&startfmodel);
        double curFreq=f[Freq];
        Fmodel fmodel(curFreq); 
        fmodel.settMatrix();
        PetscTime(&endtfmodel);
        PetscPrintf(PETSC_COMM_WORLD,"_fmodel_time=%f\n",endtfmodel-startfmodel);

        //Call the matrix element calculation and matrix assembly code
        PetscPrintf(row_comm,"开始计算v0\n");
        v0Assembly v0A(curFreq);
        v0A.creatv0(fmodel,row_comm);
        //source A
        v1Assembly v1(fmodel,curFreq);
        v1.fDirBdaries(Freq);
        v1.creatv1(v0A,fmodel,row_comm);
        //source B
        v2Assembly v2(fmodel,curFreq);
        v2.fDirBdaries(Freq);
        v2.creatv2(v0A,fmodel,row_comm);
        //Post-processing, calculation of apparent resistivity and phase
        PetscLogDouble start_t,end_t ;
        PetscTime(&start_t);
        if(row_rank==0)
        {
            cout<<"start startrho fmodel.NX="<<fmodel.NX<<endl;
            double startrho=MPI_Wtime();
            rhoAndpha rho_alpha(fmodel.NX);
            rho_alpha.rho(v1,v2,fmodel,Freq);
            double endrho=MPI_Wtime();
        } 
        PetscTime(&end_t);
        PetscPrintf(PETSC_COMM_WORLD,"Post-processing_time=%f\n",end_t-start_t);
    }

    PetscCallMPI(MPI_Comm_free(&row_comm));
    PetscTime(&end1); 
    Alltime=end1-start1;
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
