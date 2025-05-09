#include "../include/main.h"
#include "../include/fmodel.h"
using namespace std;

v2Assembly::v2Assembly(Fmodel fmodel,double f)
{
    
    w=2*M_PI*f;
    mu0=(4e-7)*M_PI;
    m_len=fmodel.NZ+1;
    datasetSize=N_MODEL_Z+1;
    //printf("m_len=%d\n",m_len);
    //闂傚倷绀侀幉锛勬暜濡ゅ啯宕查柛宀€鍎戠紞鏍煙閻楀牊绶茬紒鈧畝鍕厸闁割偆鍠曠€碉拷+1=25
    //ExA.resize(m_len,1);
    EpxB.resize(m_len,1);
    EpyB.resize(m_len,1);
    EpzB.resize(m_len,1);
    EpxB.setZero();
    EpyB.setZero();
    EpzB.setZero();
    //EpxB.resize(m_len,1);
   // EpyB.resize(m_len,1);
   // EpzB.resize(m_len,1);
    linshi.resize(datasetSize,1);
    polarization=1;
    //
    Exx2.resize(fmodel.NX);
    for(auto& row : Exx2) {
        row.resize(fmodel.NX); 
    }

    Eyy2.resize(fmodel.NX); 
    for(auto& row : Eyy2) {
        row.resize(fmodel.NX); 
    }

    Ezz2.resize(fmodel.NX); 
    for(auto& row : Ezz2) {
        row.resize(fmodel.NX); 
    }

    Hxx2.resize(fmodel.NX); 
    for(auto& row : Hxx2) {
        row.resize(fmodel.NX); 
    }

    Hyy2.resize(fmodel.NX); 
    for(auto& row : Hyy2) {
        row.resize(fmodel.NX); 
    }

    Hzz2.resize(fmodel.NX);
    for(auto& row : Hzz2) {
        row.resize(fmodel.NX); 
    }
}
v2Assembly::~v2Assembly()
{
    
}

void v2Assembly::creatv2(v0Assembly v0A,Fmodel fmodel,MPI_Comm curComm)
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
    //v2=sparse(NL,NL);
    Mat v2;
    MatCreate(curComm, &v2);
    MatSetSizes(v2, PETSC_DECIDE, PETSC_DECIDE,NL, NL);
    MatSetFromOptions(v2);
    MatSetUp(v2);
    MatDuplicate(v0A.v0, MAT_COPY_VALUES, &v2);

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
    int EpyBstart=NPML_X;
   // cout<<"EpyBstart="<<EpyBstart<<endl;      
    int EpyBend=N_MODEL_X+NPML_X+1;
    for(int i=EpyBstart;i<EpyBend;i++)//18-6+1=13
    {
        EpyB(i)=linshi(i-EpyBstart)/linshi(4);//ExA(7:19,1)=linshi(:,1)/linshi(3); 
       // printf("EpxB(i)=%e,%e linshi(i-6)=(%e,%e) linshi(2)=(%e,%e)\n",real(EpxB(i)),imag(EpxB(i)),real(linshi(i-6)),imag(linshi(i-6)),real(linshi(2)),imag(linshi(2)));
       // EpyB(i)=linshi(i-6)/linshi(2);
    }
    for(int i=0;i<EpyBstart;i++)
    {
        EpyB(i)=EpyB(EpyBstart);//ExA(1:6,1)=ExA(7,1);
        //EpyB(i)=EpyB(6);
    }
    for(int i=EpyBend;i<m_len;i++)
    {
        EpyB(i)=EpyB(EpyBend-1);//ExA(20:25,1)=ExA(19,1);
        //EpyB(i)=EpyB(18);
    }
   // for(int i=0;i<25;i++)
    {
        //printf("EpyB[i]=%e,%e\n",real(EpyB(i)),imag(EpyB(i)));
    }

    vector<complex<double>> Bp(NL+1,myzero);
   // printf("start Bp\n");
    //============EP===============  % Exp
    for(int i=1;i<=NZ+1;i++) 
    {
        for(int j=1;j<=NY+1;j++) // j=1:NY+1
        {
            for(int k=1;k<=NX;k++) //k=1:NX
            {
                int kk=(i-1)*((NX+1)*(NY+1)+(NX+1)*NY+(NY+1)*NX)+(j-1)*(NX+NX+1)+k;
                int h=floor((kk-1)/(NX*(NY+1)+(NX+1)*NY+(NX+1)*(NY+1)))+1;
                Bp[kk-1]=EpxB(h-1);//Ep(kk,1)=EpxB(h,fnn);
            
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
                Bp[kk-1]=EpyB(h-1);//Ep(kk,1)=EpyB(h,fnn);
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
                Bp[kk-1]=EpzB(h-1);
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
    //闂備礁鎼ˇ顖炴偋婵犲洤绠伴柟闂寸閸氬綊鏌ｉ悩鎻掔參濠电姷顣藉Σ鍛村磻閳ь剟鏌涚€ｎ偅宕岄柡宀嬬磿娴狅妇鎷犻幓鎺戭潥婵＄偑鍊栫敮妤冨垝鎼粹檧妲堥柣銏犲閺佸﹪鎮峰▎娆戝埌濞存粓绠栭弻銊╂偆閸屾稑顏�
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
    MatAXPY(copyK3, 1, v0A.K4, UNKNOWN_NONZERO_PATTERN); //闂傚倷鐒﹀鍧楀储閹间礁鐤鹃柣妯哄棘閻旂厧鐒垫い鎺戝閻撴洟鏌￠崘銊モ偓鍛婄墡缂備胶铏庨崢楣冣€旈敓锟�3=K3+K4
    MatMult(copyK3, Ep, P);// 闂傚倷鐒﹀鍧楀储閹间礁鐤鹃柣妯哄棘閻旂厧鐒垫い鎺嗗亾闁宠棄顦靛顕€鍩€椤掑嫭鏅梻浣筋嚃閸犳鍒掗幘璇茶摕鐎光偓閸曨剚娅㈤梺璺ㄥ櫐閹凤拷//P=(K3+K4)*Ep;  
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
    
    //PetscPrintf(PETSC_COMM_WORLD, "P=(K3+K4)*Ep闂備浇顕уù鐑藉箠閹捐瀚夋い鎺戝閸ㄥ倹鎱ㄥ┑鍥ㄢ枎");
    
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
        

        PetscLogDouble start_matv2,end_matv2;
        PetscTime(&start_matv2);
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
        MatZeroRowsColumns(v2, count_row, row, 1.0, NULL, NULL);//闂傚倷鑳堕崕鐢稿磻閹捐绀夐悗锝庡墲婵櫕銇勯幒鎴濃偓鐑芥倿婵犲洦鐓曢柣鎰摠閳绘悩闂傚倷鐒﹂惇褰掑礉瀹€鈧埀顒佸嚬閸ㄥ爼濡存笟鈧獮瀣倷绾版ɑ鐏冨┑鐘垫暩婵鈧凹鍙冮幃鐐綇閵婏箑寮块梺鎸庣箓閹虫劗绮婚敐澶嬧拺闁告瑥顦遍惌鎺斺偓瑙勬礃閸旀瑩骞婂┑瀣劷闁挎洍鍋撻柡鍜佸幘缁辨捇宕掑顒佺亾闂佸摜濮靛ú婊堝箯瑜版帗鏅搁柨鐕傛嫹0
        PetscTime(&end_matv2);
        //PetscPrintf(curComm, "MatZeroRowsColumns_V1 time(MPI):%f\n",end_matv2-start_matv2);
        
    //PetscPrintf(curComm, "==========\n");
   //闂備礁鎼ˇ顐﹀疾閳哄懎鍌ㄦ繛宸簻閼歌銇勯幒鎴濃偓褰掑窗閸℃稒鐓ユ繝闈涙椤忣偊鏌ｉ悢鍝ョ疄闁哄矉绻濋崺鈧い鎺戝閸ゅ鏌涢…鎴濅簼闁伙綇鎷�
   /* KSP ksp_bicg;
    PC pc_bicg;
    KSPCreate(curComm, &ksp_bicg);
    KSPSetType(ksp_bicg, KSPBCGS);
    KSPSetOperators(ksp_bicg, v2, v2);
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
    KSPSetOperators(ksp_bicg, v2, v2);
    KSPGetPC(ksp_bicg, &pc_bicg);
    PCSetType(pc_bicg, PCLU);
    PCFactorSetMatSolverType(pc_bicg, MATSOLVERSUPERLU_DIST);
    PetscPrintf(curComm, "=======kspSolve start======\n");
    KSPSolve(ksp_bicg, P, Xe);
    PetscTime(&endt);
    //
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
    PetscPrintf(curComm, "V2+P solve time(MPI):%f \n",endt-start);
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
        //{
            //if(bufVal[i].real()!=0.0||bufVal[i].imag()!=0)
         //   {
               // printf("bufVal[%d]=(%e,%e)\n",i,bufVal[i].real(),bufVal[i].imag());
               // countbuf++;
           // }
                
       // }
        //cout<<"countbuf="<<countbuf<<endl;
        //%%%%%%%%%%%%%%  Ex  %%%%%%%%%%%%%%% Ex2=zeros(NX,NY+1,NZ+1);
        vector<vector<vector<complex<double>>>> Ex2(NX,vector<vector<complex<double>>>(NY+1,vector<complex<double>>(NZ+1)));
        //Ey2=zeros(NX+1,NY,NZ+1);
        vector<vector<vector<complex<double>>>> Ey2(NX+1,vector<vector<complex<double>>>(NY,vector<complex<double>>(NZ+1)));
        //Ez2=zeros(NX+1,NY+1,NZ);
        vector<vector<vector<complex<double>>>> Ez2(NX+1,vector<vector<complex<double>>>(NY+1,vector<complex<double>>(NZ)));
        //printf("vector<vector<vector<complex<double>>>> Ez2\n");
        for(int k=1;k<=NZ+1;k++) 
        {
            for(int j=1;j<=NY+1;j++)
            {
                for(int i=1;i<=NX;i++) 
                {
                    int MEij=(k-1)*(NX*(NY+1)+NY*(NX+1)+(NX+1)*(NY+1))+(j-1)*(2*NX+1)+i;
                    Ex2[i-1][j-1][k-1]=bufVal[MEij-1];
                    //printf("MEij=%d \n",MEij);
                }
            }      
        }     
       // printf("Ex2\n");
        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
              //Exx2[i-1][j-1]=((Ex2(i,j+1,Nair+1)+Ex2(i,j,Nair+1))/2.d0+(Ex2(i,j+1,Nair+2)+Ex2(i,j,Nair+2))/2.d0)/2.d0;
                Exx2[i-1][j-1]=((Ex2[i-1][j][Nair]+Ex2[i-1][j-1][Nair])/2.0+(Ex2[i-1][j][Nair+1]+Ex2[i-1][j-1][Nair+1])/2.0)/2.0;

               // printf("Exx2(%d,%d)=%e,%e\n",i,j,Exx2[i-1][j-1].real(),Exx2[i-1][j-1].imag());
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
                    Ey2[i-1][j-1][k-1]=bufVal[MEij-1];
                }
            }      
        }
        //printf("Eyy2\n");
        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
              //Eyy2(i,j)=((Ey2(i,j,Nair+1)+Ey2(i+1,j,Nair+1))/2.d0+(Ey2(i,j,Nair+2)+Ey2(i+1,j,Nair+2))/2.d0)/2.d0;
                Eyy2[i-1][j-1]=((Ey2[i-1][j-1][Nair]+Ey2[i][j-1][Nair])/2.0+(Ey2[i-1][j-1][Nair+1]+Ey2[i][j-1][Nair+1])/2.0)/2.0;
                //printf("Eyy2(%d,%d)=%e,%e\n",i,j,Eyy2[i-1][j-1].real(),Eyy2[i-1][j-1].imag());
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
                    Ez2[i-1][j-1][k-1]=bufVal[MEij-1];
                }
            }      
        }
       // printf("Ezz2\n");
        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
              //Ezz2(i,j)=((Ez2(i,j,Nair+1)+Ez2(i,j+1,Nair+1))/2.d0+(Ez2(i+1,j,Nair+1)+Ez2(i+1,j+1,Nair+1))/2.d0)/2.d0;
                Ezz2[i-1][j-1]=((Ez2[i-1][j-1][Nair]+Ez2[i-1][j][Nair])/2.0+(Ez2[i][j-1][Nair]+Ez2[i][j][Nair])/2.0)/2.0;
                //printf("Ezz2(%d,%d)=%e,%e\n",i,j,Ezz2[i-1][j-1].real(),Ezz2[i-1][j-1].imag());
            }
        }
        //printf("Exz Ezx 婵犵數鍋為崹鍫曞箰閹绢喖纾婚柟鍓х帛閳锋垿鎮归幁鎺戝闁哄鍨块弻锛勨偓锝庝憾閻撳吋顨ラ悙鑼闁诡喗绮撻幊鐐哄Ψ閿旂瓔浠ч梻鍌欑閹碱偊宕愰崼鏇炵９闁哄稁鍋€閸嬫挸顫濋鍌溞ㄩ梺鍝勮閸旀垿骞冮姀銈呭窛濠电姴瀚槐鏇㈡⒒娴ｅ摜绉烘い銉︽崌瀹曟顫滈埀顒€顕ｉ锕€绠婚悹鍥у级椤ユ繈姊洪棃娑氬婵☆偅顨婇、鏃堝醇閺囩啿鎷洪柣鐘充航閸斿矂寮搁幋锔界厸閻庯綆浜堕悡鍏碱殽閻愯尙绠婚柡灞诲妿閳ь剨绲介悘姘殽閸曨垱鈷掑ù锝堝Г绾爼鏌涢悩铏暗缂侇喖锕︾紓鎾绘晸閿燂拷");
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Exz Ezx 婵犵數鍋為崹鍫曞箰閹绢喖纾婚柟鍓х帛閳锋垿鎮归幁鎺戝闁哄鍨块弻锛勨偓锝庝憾閻撳吋顨ラ悙鑼闁诡喗绮撻幊鐐哄Ψ閿旂瓔浠ч梻鍌欑閹碱偊宕愰崼鏇炵９闁哄稁鍋€閸嬫挸顫濋鍌溞ㄩ梺鍝勮閸旀垿骞冮姀銈呭窛濠电姴瀚槐鏇㈡⒒娴ｅ摜绉烘い銉︽崌瀹曟顫滈埀顒€顕ｉ锕€绠婚悹鍥у级椤ユ繈姊洪棃娑氬婵☆偅顨婇、鏃堝醇閺囩啿鎷洪柣鐘充航閸斿矂寮搁幋锔界厸閻庯綆浜堕悡鍏碱殽閻愯尙绠婚柡灞诲妿閳ь剨绲介悘姘殽閸曨垱鈷掑ù锝堝Г绾爼鏌涢敐蹇曠暤妤犵偛绻橀弫鎾绘晸閿燂拷?)  %%%%%%%%%%%%%%%%%%%%%%%%%%
        vector<vector<complex<double>>> Ezx2(NX, vector<complex<double>>(NY));
        vector<vector<complex<double>>> Exz2(NX, vector<complex<double>>(NY));

        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
                Ezx2[i-1][j-1]=(Ez2[i][j-1][Nair]-Ez2[i-1][j-1][Nair]+Ez2[i][j][Nair]-Ez2[i-1][j][Nair])/2.0/fmodel.A_X[i-1];
                Exz2[i-1][j-1]=(Ex2[i-1][j-1][Nair+1]-Ex2[i-1][j-1][Nair]+Ex2[i-1][j][Nair+1]-Ex2[i-1][j][Nair])/2.0/fmodel.C_Z[Nair];
            }
        }
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Eyx Exy 婵犵數鍋為崹鍫曞箰閹绢喖纾婚柟鍓х帛閳锋垿鎮归幁鎺戝闁哄鍨块弻锛勨偓锝庝憾閻撳吋顨ラ悙鑼闁诡喗绮撻幊鐐哄Ψ閿旂瓔浠ч梻鍌欑閹碱偊宕愰崼鏇炵９闁哄稁鍋€閸嬫挸顫濋鍌溞ㄩ梺鍝勮閸旀垿骞冮姀銈呭窛濠电姴瀚槐鏇㈡⒒娴ｅ摜绉烘い銉︽崌瀹曟顫滈埀顒€顕ｉ锕€绠婚悹鍥у级椤ユ繈姊洪棃娑氬婵☆偅顨婇、鏃堝醇閺囩啿鎷洪柣鐘充航閸斿矂寮搁幋锔界厸閻庯綆浜堕悡鍏碱殽閻愯尙绠婚柡灞诲妿閳ь剨绲介悘姘殽閸曨垱鈷掑ù锝堝Г绾爼鏌涢敐蹇曠暤妤犵偛绻橀弫鎾绘晸閿燂拷?)  %%%%%%%%%%%%%%%%%%%%%%%%%%
        vector<vector<complex<double>>> Eyx2(NX, vector<complex<double>>(NY));
        vector<vector<complex<double>>> Exy2(NX, vector<complex<double>>(NY));

        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
                Eyx2[i-1][j-1]=(Ey2[i][j-1][Nair]-Ey2[i-1][j-1][Nair]+Ey2[i][j-1][Nair+1]-Ey2[i-1][j-1][Nair+1])/2.0/fmodel.A_X[i-1];
                Exy2[i-1][j-1]=(Ex2[i-1][j][Nair]-Ex2[i-1][j-1][Nair]+Ex2[i-1][j][Nair+1]-Ex2[i-1][j-1][Nair+1])/2.0/fmodel.B_Y[j-1];
            }
        }
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Ezy Eyz 婵犵數鍋為崹鍫曞箰閹绢喖纾婚柟鍓х帛閳锋垿鎮归幁鎺戝闁哄鍨块弻锛勨偓锝庝憾閻撳吋顨ラ悙鑼闁诡喗绮撻幊鐐哄Ψ閿旂瓔浠ч梻鍌欑閹碱偊宕愰崼鏇炵９闁哄稁鍋€閸嬫挸顫濋鍌溞ㄩ梺鍝勮閸旀垿骞冮姀銈呭窛濠电姴瀚槐鏇㈡⒒娴ｅ摜绉烘い銉︽崌瀹曟顫滈埀顒€顕ｉ锕€绠婚悹鍥у级椤ユ繈姊洪棃娑氬婵☆偅顨婇、鏃堝醇閺囩啿鎷洪柣鐘充航閸斿矂寮搁幋锔界厸閻庯綆浜堕悡鍏碱殽閻愯尙绠婚柡灞诲妿閳ь剨绲介悘姘殽閸曨垱鈷掑ù锝堝Г绾爼鏌涢敐蹇曠暤妤犵偛绻橀弫鎾绘晸閿燂拷?)  %%%%%%%%%%%%%%%%%%%%%%%%%%
        vector<vector<complex<double>>> Ezy2(NX, vector<complex<double>>(NY));
        vector<vector<complex<double>>> Eyz2(NX, vector<complex<double>>(NY));

        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
                Ezy2[i-1][j-1]=(Ez2[i-1][j][Nair]-Ez2[i-1][j-1][Nair]+Ez2[i][j][Nair]-Ez2[i][j-1][Nair])/2.0/fmodel.B_Y[j-1];
                Eyz2[i-1][j-1]=(Ey2[i][j-1][Nair+1]-Ey2[i][j-1][Nair]+Ey2[i-1][j-1][Nair+1]-Ey2[i-1][j-1][Nair])/2.0/fmodel.C_Z[Nair];
            }
        }
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Res And Phase %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        complex<double> mysqrt(0,1);
        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
                Hxx2[i-1][j-1]=1/(mysqrt*w*mu0)*(Ezy2[i-1][j-1]-Eyz2[i-1][j-1]);
                Hyy2[i-1][j-1]=1/(mysqrt*w*mu0)*(Exz2[i-1][j-1]-Ezx2[i-1][j-1]);
                Hzz2[i-1][j-1]=1/(mysqrt*w*mu0)*(Eyx2[i-1][j-1]-Exy2[i-1][j-1]);
                //printf("Hxx2=%e,%e \n",Hxx2[i-1][j-1].real(),Hxx2[i-1][j-1].imag());
                //printf("Hyy2=%e,%e \n",Hyy2[i-1][j-1].real(),Hyy2[i-1][j-1].imag());
               // printf("Hzz2=%e,%e\n",Hzz2[i-1][j-1].real(),Hzz2[i-1][j-1].imag());
            }
        }
    }
    VecDestroy(&P);
    VecDestroy(&Xe);
    KSPDestroy(&ksp_bicg);
    //delete []bufVal;
    //delete  []Bp;
    MatDestroy(&v2);
}
void v2Assembly::fDirBdaries(int groupid)
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
            //Eigen::VectorXd v2=Eigen::VectorXd::Zero(13);
           // HzA.col(0)<<v2; 
            //printf("HzA\n");
        }
   // }

}
//
void v2Assembly::importdata(string f1,string f2,MatrixX &E)
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
       // printf("闂備浇顕х€涒晠顢欓弽顓炵獥闁圭儤顨呯壕濠氭煙閻愵剚鐏遍柡鈧懞銉ｄ簻闁哄啫鍊甸幏锟犳煕鎼淬垺灏伴柕鍥у瀵潙顫濋鈧粊顕€鎮楀▓鍨珮闁告挾鍠庨悾鐑藉醇閺囩喐娅嗙紓浣圭☉椤戝懏绂掔弧鐠�");
        complex<double> c1;
        if(finput1.fail()||k==datasetSize)
        {
           // printf("闂傚倷娴囧畷鍨叏閺夋嚚娲煛閸滀焦鏅悷婊勫灴婵＄敻骞囬弶璺ㄥ€為梺闈浨归崕宕囩玻閻愮儤鐓涘璺猴功婢ч亶鏌涚€ｎ亶妲兼い銈呭€垮缁樻媴閻戞ê娈岄梺瀹︽澘濮傞柟顔ㄥ洦鍊烽柛婵嗗閸旓箑顪冮妶鍡楀潑闁稿鎸鹃埀顒冾潐濞测晠骞忛敓锟�==13\n");
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