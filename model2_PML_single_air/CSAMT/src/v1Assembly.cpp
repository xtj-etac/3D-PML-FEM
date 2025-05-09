//Solve for the electric field and magnetic field of source A
#include "../include/main.h"
#include "../include/fmodel.h"
using namespace std;

v1Assembly::v1Assembly(Fmodel fmodel,double f)
{
    w=2*M_PI*f;
    mu0=(4e-7)*M_PI;
    m_len=fmodel.NZ+1;
    datasetSize=N_MODEL_Z+1;
    EpxA.resize(m_len,1);
    EpyA.resize(m_len,1);
    EpzA.resize(m_len,1);
    EpxA.setZero();
    EpyA.setZero();
    EpzA.setZero();
    linshi.resize(datasetSize,1);
    polarization=1;//Represents the A-polarization source
    Exx1.resize(fmodel.NX);
    for(auto& row : Exx1) {
        row.resize(fmodel.NX); 
    }

    Eyy1.resize(fmodel.NX); 
    for(auto& row : Eyy1) {
        row.resize(fmodel.NX); 
    }

    Ezz1.resize(fmodel.NX); 
    for(auto& row : Ezz1) {
        row.resize(fmodel.NX); 
    }

    Hxx1.resize(fmodel.NX); 
    for(auto& row : Hxx1) {
        row.resize(fmodel.NX); 
    }

    Hyy1.resize(fmodel.NX); 
    for(auto& row : Hyy1) {
        row.resize(fmodel.NX); 
    }

    Hzz1.resize(fmodel.NX);
    for(auto& row : Hzz1) {
        row.resize(fmodel.NX); 
    }
}
v1Assembly::~v1Assembly()
{
    
}

void v1Assembly::creatv1(v0Assembly v0A,Fmodel fmodel,MPI_Comm curComm)
{
    //
    PetscLogDouble start,end1,end2,end3;  
    PetscTime(&start);
    int curRank, curSize;
    MPI_Comm_size(curComm, &curSize);
    MPI_Comm_rank(curComm, &curRank);

    int start1 =MPI_Wtime();
    int NX=fmodel.NX;
    int NY=fmodel.NX;
    int NZ=fmodel.NZ;
    int NL=fmodel.NL;
    complex<double> myzero(0.0,0.0);
    complex<double> myone(1.0,0.0);
    Mat v1;
    MatCreate(curComm, &v1);
    MatSetSizes(v1, PETSC_DECIDE, PETSC_DECIDE,NL, NL);
    MatSetFromOptions(v1);
    MatSetUp(v1);
    MatDuplicate(v0A.v0, MAT_COPY_VALUES, &v1);

    Mat copyK3;
    MatCreate(curComm, &copyK3);
    MatSetSizes(copyK3, PETSC_DECIDE, PETSC_DECIDE,NL, NL);
    MatSetFromOptions(copyK3);
    MatSetUp(copyK3);
    MatDuplicate(v0A.K3, MAT_COPY_VALUES, &copyK3);
    Vec P;
    VecCreate(curComm, &P);
    VecSetSizes(P, PETSC_DECIDE,NL);
    VecSetFromOptions(P);
    VecSetUp(P);
    Vec Ep;
    VecCreate(curComm, &Ep);
    VecSetSizes(Ep, PETSC_DECIDE,NL);
    VecSetFromOptions(Ep);
    VecSetUp(Ep);
    Vec Xe;
    VecCreate(PETSC_COMM_WORLD, &Xe);
    VecSetSizes(Xe, PETSC_DECIDE, NL);
    VecSetFromOptions(Xe);
    VecSetUp(Xe);

    //引入源项Ep
    int ExpAstart=NPML_X;
    int ExpAend=N_MODEL_X+NPML_X+1;
    for(int i=ExpAstart;i<ExpAend;i++)
    {
        EpxA(i)=linshi(i-ExpAstart)/linshi(4);//normalization processing
    }
    for(int i=0;i<ExpAstart;i++)
    {
        EpxA(i)=EpxA(ExpAstart);
    }
    for(int i=ExpAend;i<m_len;i++)
    {
        EpxA(i)=EpxA(ExpAend-1);
    }

    vector<complex<double>> Bp(NL+1,myzero);
    for(int i=1;i<=NZ+1;i++) 
    {
        for(int j=1;j<=NY+1;j++) 
        {
            for(int k=1;k<=NX;k++) 
            {
                int kk=(i-1)*((NX+1)*(NY+1)+(NX+1)*NY+(NY+1)*NX)+(j-1)*(NX+NX+1)+k;
                int h=floor((kk-1)/(NX*(NY+1)+(NX+1)*NY+(NX+1)*(NY+1)))+1;
                Bp[kk-1]=EpxA(h-1);
            
            }
        }    
    }   
    for(int i=1;i<=NZ+1;i++)
    {
        for(int j=1;j<=NY;j++)
        {
            for(int k=1;k<=NX+1;k++)
            {
                int kk=(i-1)*((NX+1)*(NY+1)+(NX+1)*NY+(NY+1)*NX)+(j-1)*(NX+NX+1)+NX+k;
                int h=floor((kk-1)/(NX*(NY+1)+(NX+1)*NY+(NX+1)*(NY+1)))+1;
                Bp[kk-1]=EpyA(h-1);
            }    
        }    
    }
    for(int i=1;i<=NZ;i++)
    {
        for(int j=1;j<=NY+1;j++) 
        {
            for(int k=1;k<=NX+1;k++)
            {
                int kk=(i-1)*((NX+1)*(NY+1)+(NX+1)*NY+(NY+1)*NX)+(NX+1)*NY+(NY+1)*NX+(j-1)*(NX+1)+k;
                int h=floor((kk-1)/(NX*(NY+1)+(NX+1)*NY+(NX+1)*(NY+1)))+1;
                Bp[kk-1]=EpzA(h-1);
            }    
        }    
    }

    for(int i=1;i<=NL;i++)
    {
        VecSetValue(Ep, i-1, Bp[i-1], INSERT_VALUES);   
    }   

    VecAssemblyBegin(Ep); 
    VecAssemblyEnd(Ep);
    
    MatAXPY(copyK3, 1, v0A.K4, UNKNOWN_NONZERO_PATTERN);
    MatMult(copyK3, Ep, P);//The right end of the equation [(K3+K4)*Ep]
    
    //Calculate the index values of the edges on the six outer boundary surfaces
    int idn, idn1, idn2, idn3, idn4, idn5, idn6, idn7;
    int i,j,k;
    int BDN = ((NX+1)*NY+(NY+1)*NY)*2+NZ*(NX+1+NX+1+NY-1+NY-1+NX+NX+NY+NY)-NX-NX-NY-NY;
    Eigen::RowVectorXi Bid(BDN+1);
    Bid.setZero(); 

    idn = 0;
    for(j=1;j<=NY+1;j++)
        for(i=1;i<=NX;i++){
            idn+=1;
            Bid(idn) = (j-1)*(NX+NX+1)+i;
        }

    idn1 = idn;
    PetscPrintf(PETSC_COMM_WORLD, "A Bid(%d)\n", idn1);
    for(k=1;k<=NZ-1;k++){
        for(i=1;i<=NX;i++){
            idn+=1;
            Bid(idn) = k*(NX*(NY+1)+(NX+1)*NY+(NX+1)*(NY+1))+i;
        }
        for(i=1;i<=NX;i++){
            idn+=1;
            Bid(idn)= k*(NX*(NY+1)+(NX+1)*NY+(NX+1)*(NY+1))+NY*(NX+1)+NX*NY+i;
        }
    }

    idn2 = idn;
    PetscPrintf(PETSC_COMM_WORLD, "A Bid(%d)\n", idn2);
    for(j=1;j<=NY+1;j++)
        for(i=1;i<=NX;i++){
            idn+=1;
            Bid(idn) = NZ*(NX*(NY+1)+(NX+1)*NY+(NX+1)*(NY+1))+(j-1)*(NX+NX+1)+i;
        }

    idn3 = idn;
    PetscPrintf(PETSC_COMM_WORLD, "A Bid(%d)\n", idn3);
    for(j=1;j<=NY;j++)
        for(i=1;i<=NX+1;i++){
            idn+=1;
            Bid(idn) = (j-1)*(NX+NX+1)+NX+i;
        }

    idn4 = idn;
    PetscPrintf(PETSC_COMM_WORLD, "A Bid(%d)\n", idn4);
    for(k=1;k<=NZ-1;k++)
        for(j=1;j<=NY;j++){
            idn+=1;
            Bid(idn) = k*(NX*(NY+1)+(NX+1)*NY+(NX+1)*(NY+1))+(j-1)*(NX+NX+1)+NX+1;

            idn+=1;
            Bid(idn) = k*(NX*(NY+1)+(NX+1)*NY+(NX+1)*(NY+1))+(j-1)*(NX+NX+1)+NX+NX+1;
        }
    
    idn5 = idn;
    PetscPrintf(PETSC_COMM_WORLD, "A Bid(%d)\n", idn5);
    for(j=1;j<=NY;j++)
        for(i=1;i<=NX+1;i++){
            idn+=1;
            Bid(idn) = NZ*(NX*(NY+1)+(NX+1)*NY+(NX+1)*(NY+1))+(j-1)*(NX+NX+1)+NX+i;
        }

    idn6 = idn;
    PetscPrintf(PETSC_COMM_WORLD, "A Bid(%d)\n", idn6);
    for(k=1;k<=NZ;k++){
        for(j=1;j<=1;j++){
            for(i=1;i<=NX+1;i++){
                idn+=1;
                Bid(idn) = (k-1)*((NY+1)*NX+(NX+1)*NY+(NX+1)*(NY+1))+(NY+1)*NX+(NX+1)*NY+i;
            }
        }

        for(j=2;j<=NY;j++){
            idn+=1;
            Bid(idn) = (k-1)*((NY+1)*NX+(NX+1)*NY+(NX+1)*(NY+1))+(NY+1)*NX+(NX+1)*NY+(j-1)*(NX+1)+1;

            idn+=1;
            Bid(idn) = (k-1)*((NY+1)*NX+(NX+1)*NY+(NX+1)*(NY+1))+(NY+1)*NX+(NX+1)*NY+j*(NX+1);
        }

        for(j=NY+1;j<=NY+1;j++){
            for(i=1;i<=NX+1;i++){
                idn+=1;
                Bid(idn) = (k-1)*((NY+1)*NX+(NX+1)*NY+(NX+1)*(NY+1))+(NY+1)*NX+(NX+1)*NY+(NX+1)*NY+i;
            }
        }
    }
        idn7 = idn;
        PetscPrintf(PETSC_COMM_WORLD, "A Bid(idn7=%d)\n", idn7);
        PetscPrintf(PETSC_COMM_WORLD, "A 2\n");
        //PEC boundary conditions

        int slide, st_idx, ed_idx;
        slide = floor((double)idn7 / curSize);
        st_idx = curRank*slide+1;
        if(curRank != curSize-1)
            ed_idx = (curRank+1)*slide;
        else
            ed_idx = idn7;
        for(i=st_idx;i<=ed_idx;i++)
        {   
            VecSetValue(P, Bid(i)-1, myzero, INSERT_VALUES);
        }
        VecAssemblyBegin(P);
        VecAssemblyEnd(P);
        VecDuplicate(P, &Xe); 
        
        int count_row=0;
        PetscInt row[ed_idx-st_idx+1];
        PetscPrintf(curComm, "ed_idx-st_idx:%d\n",ed_idx-st_idx+1);
        for(int i=st_idx;i<=ed_idx;i++)
        {   

            PetscInt base_i = Bid(i) - 1;
            row[count_row]=base_i;
            count_row++;
        }     
        MatZeroRowsColumns(v1, count_row, row, 1.0, NULL, NULL);
    //Direct solver solves the quadratic field
    PetscTime(&end1);
    KSP ksp_bicg;
    PC pc_bicg;
    KSPCreate(curComm, &ksp_bicg);
    KSPSetType(ksp_bicg, KSPPREONLY);
    KSPSetOperators(ksp_bicg, v1, v1);
    KSPGetPC(ksp_bicg, &pc_bicg);
    PCSetType(pc_bicg, PCLU);
    PCFactorSetMatSolverType(pc_bicg, MATSOLVERSUPERLU_DIST);
    PetscPrintf(curComm, "=======kspSolve start======\n");
    KSPSolve(ksp_bicg, P, Xe);
    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp_bicg, &reason);
    if (reason == KSP_DIVERGED_INDEFINITE_PC) {
        PetscPrintf(PETSC_COMM_WORLD, "\nDivergence because of indefinite preconditioner;\n");
        PetscPrintf(PETSC_COMM_WORLD, "Run the executable again but with -pc_factor_shift_positive_definite option.\n");
    } else if (reason < 0) {
        PetscPrintf(PETSC_COMM_WORLD, "\nOther kind of divergence: this should not happen.-->%d\n",reason);
    } else{
        PetscPrintf(PETSC_COMM_WORLD, "\nOther  not happen.\n");
    }
    PetscTime(&end2);
    PetscPrintf(curComm, "A solve time(MPI):%f \n",end2-end1);
    //Calculate the total area
    VecAXPY(Xe,1,Ep);//E=Ep+Es;
    
    PetscScalar *bufA;
    PetscInt Xe_s;
    VecGetArray(Xe, &bufA);
    VecGetLocalSize(Xe,&Xe_s);
    int recv_curSize[curSize]={0};
    complex<double> bufVal[NL]={myzero};
    int displs[curSize]={0};

    MPI_Gather(&Xe_s,1,MPIU_INT,recv_curSize,1,MPIU_INT,0,curComm);
    if(!curRank)
    {
        displs[0]=0;
        for(int i=1;i<curSize;i++)
        {
             displs[i]=displs[i-1]+recv_curSize[i-1];
        }
    }
     //
    MPI_Gatherv(bufA,Xe_s,MPIU_SCALAR,bufVal,recv_curSize,displs,MPIU_SCALAR,0,curComm);
   
    VecRestoreArray(Xe, &bufA);
    if(curRank==0)
    {
        int countbuf=0;
        {
            {
            }
                
        }
        cout<<"countbuf="<<countbuf<<endl;
        vector<vector<vector<complex<double>>>> Ex1(NX,vector<vector<complex<double>>>(NY+1,vector<complex<double>>(NZ+1)));
        vector<vector<vector<complex<double>>>> Ey1(NX+1,vector<vector<complex<double>>>(NY,vector<complex<double>>(NZ+1)));
        vector<vector<vector<complex<double>>>> Ez1(NX+1,vector<vector<complex<double>>>(NY+1,vector<complex<double>>(NZ)));
        //Convert the node values into edge values
        // %%%%%%%%% Ex %%%%%%%%%%%%%
        for(int k=1;k<=NZ+1;k++) 
        {
            for(int j=1;j<=NY+1;j++)
            {
                for(int i=1;i<=NX;i++) 
                {
                    int MEij=(k-1)*(NX*(NY+1)+NY*(NX+1)+(NX+1)*(NY+1))+(j-1)*(2*NX+1)+i;
                    Ex1[i-1][j-1][k-1]=bufVal[MEij-1];
                }
            }      
        }     
        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
                Exx1[i-1][j-1]=((Ex1[i-1][j][Nair]+Ex1[i-1][j-1][Nair])/2.0+(Ex1[i-1][j][Nair+1]+Ex1[i-1][j-1][Nair+1])/2.0)/2.0;
            }
        }
        // %%%%%%%%% Ey %%%%%%%%%%%%%
        for(int k=1;k<=NZ+1;k++) 
        {
            for(int i=1;i<=NX+1;i++)
            {
                for(int j=1;j<=NY;j++) 
                {
                    int MEij=(k-1)*(NX*(NY+1)+NY*(NX+1)+(NX+1)*(NY+1))+(j-1)*(2*NX+1)+NX+i;
                    Ey1[i-1][j-1][k-1]=bufVal[MEij-1];
                }
            }      
        }
        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
                Eyy1[i-1][j-1]=((Ey1[i-1][j-1][Nair]+Ey1[i][j-1][Nair])/2.0+(Ey1[i-1][j-1][Nair+1]+Ey1[i][j-1][Nair+1])/2.0)/2.0;
            }
        }
        //%%%%%%%%% Ez %%%%%%%%%%%%%
        for(int k=1;k<=NZ;k++) 
        {
            for(int i=1;i<=NX+1;i++)
            {
                for(int j=1;j<=NY+1;j++) 
                {
                    int MEij=(k-1)*(NX*(NY+1)+NY*(NX+1)+(NX+1)*(NY+1))+NX*(NY+1)+NY*(NX+1)+(j-1)*(NX+1)+i;
                    Ez1[i-1][j-1][k-1]=bufVal[MEij-1];
                }
            }      
        }
        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
                Ezz1[i-1][j-1]=((Ez1[i-1][j-1][Nair]+Ez1[i-1][j][Nair])/2.0+(Ez1[i][j-1][Nair]+Ez1[i][j][Nair])/2.0)/2.0;
            }
        }

        //Start calculating the partial derivatives in order to perform the subsequent calculations of the magnetic field.
        vector<vector<complex<double>>> Ezx1(NX, vector<complex<double>>(NY));
        vector<vector<complex<double>>> Exz1(NX, vector<complex<double>>(NY));

        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
                Ezx1[i-1][j-1]=(Ez1[i][j-1][Nair]-Ez1[i-1][j-1][Nair]+Ez1[i][j][Nair]-Ez1[i-1][j][Nair])/2.0/fmodel.A_X[i-1];
                Exz1[i-1][j-1]=(Ex1[i-1][j-1][Nair+1]-Ex1[i-1][j-1][Nair]+Ex1[i-1][j][Nair+1]-Ex1[i-1][j][Nair])/2.0/fmodel.C_Z[Nair];
            }
        }
        
        vector<vector<complex<double>>> Eyx1(NX, vector<complex<double>>(NY));
        vector<vector<complex<double>>> Exy1(NX, vector<complex<double>>(NY));

        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
                Eyx1[i-1][j-1]=(Ey1[i][j-1][Nair]-Ey1[i-1][j-1][Nair]+Ey1[i][j-1][Nair+1]-Ey1[i-1][j-1][Nair+1])/2.0/fmodel.A_X[i-1];
                Exy1[i-1][j-1]=(Ex1[i-1][j][Nair]-Ex1[i-1][j-1][Nair]+Ex1[i-1][j][Nair+1]-Ex1[i-1][j-1][Nair+1])/2.0/fmodel.B_Y[j-1];
            }
        }

        vector<vector<complex<double>>> Ezy1(NX, vector<complex<double>>(NY));
        vector<vector<complex<double>>> Eyz1(NX, vector<complex<double>>(NY));

        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
                Ezy1[i-1][j-1]=(Ez1[i-1][j][Nair]-Ez1[i-1][j-1][Nair]+Ez1[i][j][Nair]-Ez1[i][j-1][Nair])/2.0/fmodel.B_Y[j-1];
                Eyz1[i-1][j-1]=(Ey1[i][j-1][Nair+1]-Ey1[i][j-1][Nair]+Ey1[i-1][j-1][Nair+1]-Ey1[i-1][j-1][Nair])/2.0/fmodel.C_Z[Nair];
            }
        }
        //Calculate the magnetic field
        complex<double> mysqrt(0,1);
        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
                Hxx1[i-1][j-1]=1/(mysqrt*w*mu0)*(Ezy1[i-1][j-1]-Eyz1[i-1][j-1]);
                Hyy1[i-1][j-1]=1/(mysqrt*w*mu0)*(Exz1[i-1][j-1]-Ezx1[i-1][j-1]);
                Hzz1[i-1][j-1]=1/(mysqrt*w*mu0)*(Eyx1[i-1][j-1]-Exy1[i-1][j-1]);
            }
        }
    }
    VecDestroy(&P);
    VecDestroy(&Xe);
    KSPDestroy(&ksp_bicg);
    MatDestroy(&v1);
    PetscTime(&end3);
    PetscPrintf(curComm, "A_Ep and boundary time(MPI):%f \n",end3-start-end2+end1);
}
void v1Assembly::fDirBdaries(int groupid)
{
    string s=to_string(groupid);
    string f="./data/dataset"+s;
        if(polarization==1)
        {
            importdata(f+"/P2REx.dat",f+"/P2IEx.dat",linshi);//Read the complex value of the plane wave source
        }
}

void v1Assembly::importdata(string f1,string f2,MatrixX &E)
{
    fstream finput1(f1);
    ifstream finput2(f2); 
    string s1;
    string s2;
    int k=0;

    if (!finput1.is_open()) {
        std::cerr << "Failed to open file: " << f1 << std::endl;
    }

    while(!finput1.eof())
    {
        complex<double> c1;
        if(finput1.fail()||k==datasetSize)
        {
            break;
        }
        finput1>>s1;
        finput2>>s2;
        c1.real(stod(s1));
        c1.imag(stod(s2));
        E(k)=c1;
        k++;
    }
    finput1.close();
    finput2.close();
}