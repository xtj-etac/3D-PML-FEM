#include "../include/main.h"
#include "../include/fmodel.h"
using namespace std;

v1Assembly::v1Assembly(Fmodel fmodel,double f)
{
    
    w=2*M_PI*f;
    mu0=(4e-7)*M_PI;
    //int m_len=fmodel.NZ+1;
    m_len=fmodel.NZ+1;
    
   // printf("m_len=%d\n",m_len);
    //闂佸憡甯楃换鍌烇綖閹版澘绀岄柡鍫熻寖+1=25
    //ExA.resize(m_len,1);
    EpxA.resize(m_len,1);
    EpyA.resize(m_len,1);
    EpzA.resize(m_len,1);
    EpxA.setZero();
    EpyA.setZero();
    EpzA.setZero();
    //EpxB.resize(m_len,1);
   // EpyB.resize(m_len,1);
   // EpzB.resize(m_len,1);
    linshi.resize(m_len,1);
    polarization=1;
    //
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
    int curRank, curSize;
    MPI_Comm_size(curComm, &curSize);
    MPI_Comm_rank(curComm, &curRank);

    int start1 =MPI_Wtime();
   // int NP=fmodel.NP;
    int NX=fmodel.NX;
    int NY=fmodel.NX;
    int NZ=fmodel.NZ;
    int NL=fmodel.NL;
    //int Nair=fmodel.Nair;
    
    complex<double> myzero(0.0,0.0);
    complex<double> myone(1.0,0.0);
    PetscLogDouble start,endt;  
    PetscTime(&start);
       // cout<<"source A ??"<<endl;
    //v1=sparse(NL,NL);
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
    //P=sparse(NL,1);
    Vec P;
    VecCreate(curComm, &P);
    VecSetSizes(P, PETSC_DECIDE,NL);
    VecSetFromOptions(P);
    VecSetUp(P);//VecSet(P,myzero);
    
    //Ep=sparse(NL,1);
    Vec Ep;
    VecCreate(curComm, &Ep);
    VecSetSizes(Ep, PETSC_DECIDE,NL);
    VecSetFromOptions(Ep);
    VecSetUp(Ep);
    //X1=zeros(1,NL);
    Vec Xe;
    VecCreate(PETSC_COMM_WORLD, &Xe);
    VecSetSizes(Xe, PETSC_DECIDE, NL);
    VecSetFromOptions(Xe);
    VecSetUp(Xe);

    //fDir
    for(int i=0;i<m_len;i++)//18-6+1=13 47=NZ+1
    {
        EpxA(i)=linshi(i)/linshi(Nair);//ExA(7:19,1)=linshi(:,1)/linshi(3); 
       // printf("EpxA(i)=%e,%e linshi(i-6)=(%e,%e) linshi(2)=(%e,%e)\n",real(EpxA(i)),imag(EpxA(i)),real(linshi(i-6)),imag(linshi(i-6)),real(linshi(2)),imag(linshi(2)));
       // EpyB(i)=linshi(i-6)/linshi(2);
    }
    //for(int i=0;i<47;i++)
    {
       // printf("EpxA[i]=%e,%e\n",real(EpxA(i)),imag(EpxA(i)));
    }

    vector<complex<double>> Bp(NL+1,myzero);
    //printf("start Bp\n");
    //============EP===============  % Exp
    for(int i=1;i<=NZ+1;i++) 
    {
        for(int j=1;j<=NY+1;j++) // j=1:NY+1
        {
            for(int k=1;k<=NX;k++) //k=1:NX
            {
                int kk=(i-1)*((NX+1)*(NY+1)+(NX+1)*NY+(NY+1)*NX)+(j-1)*(NX+NX+1)+k;
                int h=floor((kk-1)/(NX*(NY+1)+(NX+1)*NY+(NX+1)*(NY+1)))+1;
                Bp[kk-1]=EpxA(h-1);//Ep(kk,1)=EpxA(h,fnn);
            
            }
        }    
    }   
    for(int i=1;i<=NZ+1;i++)// i=1:NZ+1  % Eyp
    {
        for(int j=1;j<=NY;j++)// j=1:NY
        {
            for(int k=1;k<=NX+1;k++)// k=1:NX+1
            {
                int kk=(i-1)*((NX+1)*(NY+1)+(NX+1)*NY+(NY+1)*NX)+(j-1)*(NX+NX+1)+NX+k;
                int h=floor((kk-1)/(NX*(NY+1)+(NX+1)*NY+(NX+1)*(NY+1)))+1;
                Bp[kk-1]=EpyA(h-1);//Ep(kk,1)=EpyA(h,fnn);
            }    
        }    
    }
    for(int i=1;i<=NZ;i++)//% Ezp
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
    //int countEp=0;
    for(int i=1;i<=NL;i++)
    {
        VecSetValue(Ep, i-1, Bp[i-1], INSERT_VALUES);
        //if(curRank==0)
        //if(Bp[i-1].real()!=0.0||Bp[i-1].imag()!=0.0)
        //{
           // printf("Bp[%d]=%e,%e\n",i-1,Bp[i-1].real(),Bp[i-1].imag());
          //  countEp++;
       // }
        
    }   
    //cout<<"countEp="<<countEp<<endl;
    VecAssemblyBegin(Ep); 
    VecAssemblyEnd(Ep);
    //闁哄鐗婇幐鎼佸吹閻犲埦濠碘槅鍋€閸嬫捇鏌＄仦璇插姕妞ゆ帗绮庡☉鐢割敊鐟欙絽浜鹃柨鐕傛嫹
   /* PetscScalar *bufA1;
    PetscInt sizeEp;
    VecGetLocalSize(Ep,&sizeEp);
    VecGetArray(Ep, &bufA1);
    int countEp=0;
    //if(rank>=70&&rank<=80)
    {
        for(int i=0;i<sizeEp;i++)
        {
            if(bufA1[i].real()!=0.0||bufA1[i].imag()!=0.0)
            {
                printf( "Ep,i=%d=(%e,%e)\n",i,bufA1[i].real(),bufA1[i].imag());
                countEp++;
            }
                
        }
            
    }
    cout<<"coutEp="<<countEp<<endl;
    VecRestoreArray(Ep, &bufA1);*/
    MatAXPY(copyK3, 1, v0A.K4, UNKNOWN_NONZERO_PATTERN); //闂佹椿鍘搁弲鐘参熼埀顒勬煕閺冨倸校缁绢厾顪�3=K3+K4
    MatMult(copyK3, Ep, P);// 闂佹椿鍘搁弲鐘参熼埀顒€鈽夐弮鈧敋闁诡喗绮撻弻宀勬晸閿燂拷//P=(K3+K4)*Ep;  
   /* PetscScalar *bufA1;
    PetscInt sizeP;
    VecGetLocalSize(P,&sizeP);
    VecGetArray(P, &bufA1);
    int countP=0;
    //if(rank>=70&&rank<=80)
    {
        for(int i=0;i<sizeP;i++)
        {
            if(bufA1[i].real()!=0.0||bufA1[i].imag()!=0.0)
            {
                printf( "P,i=%d=(%e,%e)\n",i,bufA1[i].real(),bufA1[i].imag());
                countP++;
            }
                
        }
            
    }
    cout<<"coutEp="<<countP<<endl;
    VecRestoreArray(P, &bufA1);*/
    
    //PetscPrintf(PETSC_COMM_WORLD, "P=(K3+K4)*Ep");
    
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
        

        PetscLogDouble start_matv1,end_matv1;
        PetscTime(&start_matv1);
        int count_row=0;
        PetscInt row[ed_idx-st_idx+1];
        PetscPrintf(curComm, "ed_idx-st_idx:%d\n",ed_idx-st_idx+1);
        for(int i=st_idx;i<=ed_idx;i++)
        {   

            PetscInt base_i = Bid(i) - 1;
            row[count_row]=base_i;
            count_row++;
            //PetscPrintf(curComm,"i =%d ,base_i=%d\n",i,base_i );
        }     
        MatZeroRowsColumns(v1, count_row, row, 1.0, NULL, NULL);//闂佺儵鍋撻崝瀣姳椤掑倷鐒婇柛鐐存▓w闂佹眹鍔岀€氼垶銆侀幋鐐碘枖濠电姴瀚悾杈ㄦ叏閿濆懐绠ｉ梺鍙夌矒瀹曟劙鎳濋柇锕€鏅紓鍌氬枤閸犳洜鎷归敓锟�0
        PetscTime(&end_matv1);
       // PetscPrintf(curComm, "MatZeroRowsColumns_V1 time(MPI):%f\n",end_matv1-start_matv1);
        
    
   //闁哄鏅╅崢娲船椤掑倹鍠嗛柨婵嗗閻熸繈鏌￠埀顒勬嚋閸偅鐣�
   /* KSP ksp_bicg;
    PC pc_bicg;
    KSPCreate(curComm, &ksp_bicg);
    KSPSetType(ksp_bicg, KSPBCGS);
    KSPSetOperators(ksp_bicg, v1, v1);
    PetscInt    its;
    KSPGetPC(ksp_bicg, &pc_bicg);
    PCSetType(pc_bicg, PCJACOBI); 
    //PCFactorSetMatSolverType(pc_bicg, MATSOLVERSUPERLU_DIST);//
    KSPSetTolerances(ksp_bicg,1e-14, PETSC_DEFAULT, PETSC_DEFAULT, 10000);
    KSPSetFromOptions(ksp_bicg);

    PetscPrintf(curComm, "=======kspSolve start======\n");
    KSPSolve(ksp_bicg, P, Xe);
    KSPGetIterationNumber(ksp_bicg, &its);*/

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
    PetscTime(&endt);
    //闂佸吋鍎抽崲鑼躲亹閸ヮ剙纾绘繝濠傚閸撶钉SP闁哄鏅╅崢娲船椤掑嫭鍎嶉柛鏇ㄥ亜閺傃囨煕閵壯冃為柍褜鍓ㄩ幏锟�
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
    PetscTime(&endt);
    PetscPrintf(curComm, "V1+P solve time(MPI):%f \n",endt-start);
    VecAXPY(Xe,1,Ep);//X1(:)=X1(:)+Ep;
    
    PetscScalar *bufA;
    PetscInt Xe_s;
    VecGetArray(Xe, &bufA);
    VecGetLocalSize(Xe,&Xe_s);
    //int *recv_curSize=new int[curSize];
    int recv_curSize[curSize]={0};
    //complex<double> *bufVal=new complex<double>[NL];
    //vector<complex<double>> bufVal(NL);
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
   //delete[] recv_curSize;
    //if(curRank==0)
        //PetscPrintf(PETSC_COMM_SELF, "curRank=%d Xe=(%e,%e)  curSize=%d\n",curRank,bufVal[1].real(),bufVal[1].imag(),Xe_s);
    VecRestoreArray(Xe, &bufA);
    if(curRank==0)
    {
        int countbuf=0;
        //for(int i=0; i<NL; i++)
        {
            //if(bufVal[i].real()!=0.0||bufVal[i].imag()!=0)
            {
               // printf("bufVal[%d]=(%e,%e)\n",i,bufVal[i].real(),bufVal[i].imag());
               // countbuf++;
            }
                
        }
       // cout<<"countbuf="<<countbuf<<endl;
        //%%%%%%%%%%%%%%  Ex  %%%%%%%%%%%%%%% Ex1=zeros(NX,NY+1,NZ+1);
        vector<vector<vector<complex<double>>>> Ex1(NX,vector<vector<complex<double>>>(NY+1,vector<complex<double>>(NZ+1)));
        //Ey1=zeros(NX+1,NY,NZ+1);
        vector<vector<vector<complex<double>>>> Ey1(NX+1,vector<vector<complex<double>>>(NY,vector<complex<double>>(NZ+1)));
        //Ez1=zeros(NX+1,NY+1,NZ);
        vector<vector<vector<complex<double>>>> Ez1(NX+1,vector<vector<complex<double>>>(NY+1,vector<complex<double>>(NZ)));
        //printf("vector<vector<vector<complex<double>>>> Ez1\n");
        for(int k=1;k<=NZ+1;k++) 
        {
            for(int j=1;j<=NY+1;j++)
            {
                for(int i=1;i<=NX;i++) 
                {
                    int MEij=(k-1)*(NX*(NY+1)+NY*(NX+1)+(NX+1)*(NY+1))+(j-1)*(2*NX+1)+i;
                    Ex1[i-1][j-1][k-1]=bufVal[MEij-1];
                    //printf("MEij=%d \n",MEij);
                }
            }      
        }     
       // printf("Ex1\n");
        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
              //Exx1[i-1][j-1]=((Ex1(i,j+1,Nair+1)+Ex1(i,j,Nair+1))/2.d0+(Ex1(i,j+1,Nair+2)+Ex1(i,j,Nair+2))/2.d0)/2.d0;
                Exx1[i-1][j-1]=((Ex1[i-1][j][Nair]+Ex1[i-1][j-1][Nair])/2.0+(Ex1[i-1][j][Nair+1]+Ex1[i-1][j-1][Nair+1])/2.0)/2.0;

               // printf("Exx1(%d,%d)=%e,%e\n",i,j,Exx1[i-1][j-1].real(),Exx1[i-1][j-1].imag());
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
        //printf("Eyy1\n");
        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
              //Eyy1(i,j)=((Ey1(i,j,Nair+1)+Ey1(i+1,j,Nair+1))/2.d0+(Ey1(i,j,Nair+2)+Ey1(i+1,j,Nair+2))/2.d0)/2.d0;
                Eyy1[i-1][j-1]=((Ey1[i-1][j-1][Nair]+Ey1[i][j-1][Nair])/2.0+(Ey1[i-1][j-1][Nair+1]+Ey1[i][j-1][Nair+1])/2.0)/2.0;
                //printf("Eyy1(%d,%d)=%e,%e\n",i,j,Eyy1[i-1][j-1].real(),Eyy1[i-1][j-1].imag());
            }
        }
        //printf("Ez\n");
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
       // printf("Ezz1\n");
        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
              //Ezz1(i,j)=((Ez1(i,j,Nair+1)+Ez1(i,j+1,Nair+1))/2.d0+(Ez1(i+1,j,Nair+1)+Ez1(i+1,j+1,Nair+1))/2.d0)/2.d0;
                Ezz1[i-1][j-1]=((Ez1[i-1][j-1][Nair]+Ez1[i-1][j][Nair])/2.0+(Ez1[i][j-1][Nair]+Ez1[i][j][Nair])/2.0)/2.0;
                //printf("Ezz1(%d,%d)=%e,%e\n",i,j,Ezz1[i-1][j-1].real(),Ezz1[i-1][j-1].imag());
            }
        }
        //printf("Exz Ezx 婵炴垶鎸撮崑鎾绘⒑鐠恒劌鏋戦柡瀣煼楠炲繘鎮滈懞銉︽闂佸搫鍊堕崐鏍偓姘秺閺屻劑鎮㈤崨濠勪紕闂佸綊顥撻崗姗€寮幘璇叉闁靛牆妫楅鍫曟⒑鐠恒劌鏋戦柡瀣煼楠炲繘鏌ㄧ€ｅ灚楠勯梻浣芥硶閸犳挾绮婄缓锟�");
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Exz Ezx 婵炴垶鎸撮崑鎾绘⒑鐠恒劌鏋戦柡瀣煼楠炲繘鎮滈懞銉︽闂佸搫鍊堕崐鏍偓姘秺閺屻劑鎮㈤崨濠勪紕闂佸綊顥撻崗姗€寮幘璇叉闁靛牆妫楅鍫曟⒑鐠恒劌鏋戦柡瀣煼楠炲繘鏌ㄧ€ｅ灚楠勯梻浣芥硶閸ｏ箓骞忛敓锟�?)  %%%%%%%%%%%%%%%%%%%%%%%%%%
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
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Eyx Exy 婵炴垶鎸撮崑鎾绘⒑鐠恒劌鏋戦柡瀣煼楠炲繘鎮滈懞銉︽闂佸搫鍊堕崐鏍偓姘秺閺屻劑鎮㈤崨濠勪紕闂佸綊顥撻崗姗€寮幘璇叉闁靛牆妫楅鍫曟⒑鐠恒劌鏋戦柡瀣煼楠炲繘鏌ㄧ€ｅ灚楠勯梻浣芥硶閸ｏ箓骞忛敓锟�?)  %%%%%%%%%%%%%%%%%%%%%%%%%%
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
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Ezy Eyz 婵炴垶鎸撮崑鎾绘⒑鐠恒劌鏋戦柡瀣煼楠炲繘鎮滈懞銉︽闂佸搫鍊堕崐鏍偓姘秺閺屻劑鎮㈤崨濠勪紕闂佸綊顥撻崗姗€寮幘璇叉闁靛牆妫楅鍫曟⒑鐠恒劌鏋戦柡瀣煼楠炲繘鏌ㄧ€ｅ灚楠勯梻浣芥硶閸ｏ箓骞忛敓锟�?)  %%%%%%%%%%%%%%%%%%%%%%%%%%
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
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Res And Phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        complex<double> mysqrt(0,1);
        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
                Hxx1[i-1][j-1]=1/(mysqrt*w*mu0)*(Ezy1[i-1][j-1]-Eyz1[i-1][j-1]);
                Hyy1[i-1][j-1]=1/(mysqrt*w*mu0)*(Exz1[i-1][j-1]-Ezx1[i-1][j-1]);
                Hzz1[i-1][j-1]=1/(mysqrt*w*mu0)*(Eyx1[i-1][j-1]-Exy1[i-1][j-1]);
                //printf("Hxx1=%e,%e \n",Hxx1[i-1][j-1].real(),Hxx1[i-1][j-1].imag());
                //printf("Hyy1=%e,%e \n",Hyy1[i-1][j-1].real(),Hyy1[i-1][j-1].imag());
               // printf("Hzz1=%e,%e\n",Hzz1[i-1][j-1].real(),Hzz1[i-1][j-1].imag());
            }
        }
    }
    VecDestroy(&P);
    VecDestroy(&Xe);
    KSPDestroy(&ksp_bicg);
    //delete []bufVal;
    //delete  []Bp;
    MatDestroy(&v1);
}
void v1Assembly::fDirBdaries(int groupid)
{
    string s=to_string(groupid);//0/1/2/3
    string f="./data/dataset"+s;
    //{
        if(polarization==1)
        {
            
            importdata(f+"/P2REx.dat",f+"/P2IEx.dat",linshi);
            //importdata(f+"/P2REy.dat",f+"/P2IEy.dat",EyA,m_len);
            //importdata(f+"/P2REz.dat",f+"/P2IEz.dat",EzA,m_len);
            //importdata(f+"/P2RHx.dat",f+"/P2IHx.dat",HxA,m_len);
            //importdata(f+"/P2RHy.dat",f+"/P2IHy.dat",HyA,m_len);
            //Eigen::VectorXd v1=Eigen::VectorXd::Zero(13);
           // HzA.col(0)<<v1; 
            //printf("HzA\n");
        }
   // }

}
//
void v1Assembly::importdata(string f1,string f2,MatrixX &E)
{
    fstream finput1(f1);
    ifstream finput2(f2); 
    string s1;
    string s2;
    int k=0;
    //cout<<f1<<endl;
    //cout<<f2<<endl;
    if (!finput1.is_open()) {
        std::cerr << "Failed to open file: " << f1 << std::endl;
       
    }

    while(!finput1.eof())
    {
       // printf("闁诲孩顔栭崰鎺楀磻閹炬枼鏀芥い鏃傗拡閸庢挻銇勯弴姘祮鐎殿噮鍓熷畷鍫曟晜缁涘浠篭n");
        complex<double> c1;
        if(finput1.fail()||k==m_len)
        {
           // printf("闂佽崵濮村ú鈺咁敋瑜戦妵鎰板炊閵娿儳绐為柡澶婄墱閸嬪顤傞梻浣烘嚀閻°劑鎮ч悙鍝勭劦妞ゆ巻鍋撶€规洩鎷�==13\n");
            break;
        }
        finput1>>s1;
        finput2>>s2;
        c1.real(stod(s1));
        c1.imag(stod(s2));
        E(k)=c1;
        //printf("c1(%d)=(%e,%e)\n",k,real(c1),imag(c1));
        k++;
    }
    finput1.close();
    finput2.close();
}