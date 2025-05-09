#include "../include/main.h"
#include "../include/fmodel.h"
using namespace std;

v2Assembly::v2Assembly(Fmodel fmodel,double f)
{
    
    w=2*M_PI*f;
    mu0=(4e-7)*M_PI;
    m_len=fmodel.NZ+1;
    //printf("m_len=%d\n",m_len);
    //闂備礁鎲＄敮妤冩崲閸岀儑缍栭柟鐗堟緲缁€宀勬煛閸喕瀵�+1=25
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
    linshi.resize(m_len,1);
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
    for(int i=0;i<m_len;i++)//18-6+1=13
    {
        EpyB(i)=linshi(i)/linshi(Nair);//ExA(7:19,1)=linshi(:,1)/linshi(3); 
       // printf("EpxA(i)=%e,%e linshi(i-6)=(%e,%e) linshi(2)=(%e,%e)\n",real(EpxA(i)),imag(EpxA(i)),real(linshi(i-6)),imag(linshi(i-6)),real(linshi(2)),imag(linshi(2)));
       // EpyB(i)=linshi(i-6)/linshi(2);
    }
   /* for(int i=0;i<47;i++)
    {
        printf("EpyB[i]=%e,%e\n",real(EpyB(i)),imag(EpyB(i)));
    }*/
   /* for(int i=6;i<19;i++)//18-6+1=13
    {
        EpyB(i)=linshi(i-6)/linshi(2);//ExA(7:19,1)=linshi(:,1)/linshi(3); 
       // printf("EpxB(i)=%e,%e linshi(i-6)=(%e,%e) linshi(2)=(%e,%e)\n",real(EpxB(i)),imag(EpxB(i)),real(linshi(i-6)),imag(linshi(i-6)),real(linshi(2)),imag(linshi(2)));
       // EpyB(i)=linshi(i-6)/linshi(2);
    }
    for(int i=0;i<6;i++)
    {
        EpyB(i)=EpyB(6);//ExA(1:6,1)=ExA(7,1);
        //EpyB(i)=EpyB(6);
    }
    for(int i=19;i<25;i++)
    {
        EpyB(i)=EpyB(18);//ExA(20:25,1)=ExA(19,1);
        //EpyB(i)=EpyB(18);
    }*/
    /*for(int i=0;i<25;i++)
    {
        printf("EpxB[i]=%e,%e\n",real(EpxB(i)),imag(EpxB(i)));
    }*/

    vector<complex<double>> Bp(NL+1,myzero);
    printf("start Bp\n");
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
    //闂佸搫顦悧濠囧箰閹间礁鍚归柣鐘插煢婵犵妲呴崑鈧柛瀣崌閺岋紕浠︾拠鎻掑濡炪倖甯楃划搴♀槈閻㈠壊鏁婇悷娆欑到娴滈箖鏌ㄩ悤鍌涘
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
    MatAXPY(copyK3, 1, v0A.K4, UNKNOWN_NONZERO_PATTERN); //闂備焦妞块崢鎼佸疾閻樺弬鐔煎焵椤掑嫭鐓曢柡鍐ㄥ€告牎缂佺虎鍘鹃—锟�3=K3+K4
    MatMult(copyK3, Ep, P);// 闂備焦妞块崢鎼佸疾閻樺弬鐔煎焵椤掆偓閳藉寮埀顒勬晪闂佽鍠楃划鎾诲蓟瀹€鍕櫢闁跨噦鎷�//P=(K3+K4)*Ep;  
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
    
    //PetscPrintf(PETSC_COMM_WORLD, "P=(K3+K4)*Ep闂佽娴烽幊鎾诲嫉椤掑嫬鍨傛慨婵囨▎");
    
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
        MatZeroRowsColumns(v2, count_row, row, 1.0, NULL, NULL);//闂備胶鍎甸崑鎾诲礉鐎ｎ剝濮虫い鎺戝€烽悞濠囨煕閻愬瓨鈻搘闂備焦鐪归崝宀€鈧凹鍨堕妴渚€骞嬮悙纰樻灃婵犵數濮寸€氼參鎮炬潏銊﹀弿闁挎繂鎳愮粻锝夋⒑閸欏鐭掔€规洘鍔欓幊婵嬫焽閿曗偓閺咁厾绱撻崒姘灓闁哥姵娲滈幏褰掓晸閿燂拷0
        PetscTime(&end_matv2);
        //PetscPrintf(curComm, "MatZeroRowsColumns_V1 time(MPI):%f\n",end_matv2-start_matv2);
        
    //PetscPrintf(curComm, "====闁诲孩顔栭崰鎺楀磻閹炬枼鏀芥い鏃傗拡閸庢劙鏌ｉ鍛棆闁瑰弶鎸抽弫鎾绘晸閿燂拷======\n");
   //闂佸搫顦弲鈺呭储濞差亜鑸规い鎺戝€归崰鍡涙煥濠靛棗顏柣鐔哥箞閺岋繝鍩€椤掑嫭鍤嬮柛顭戝亝閻ｏ拷
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
    //闂備礁鍚嬮崕鎶藉床閼艰翰浜归柛銉墮绾剧粯绻濇繝鍌氼伌闁告挾閽塖P闂佸搫顦弲鈺呭储濞差亜鑸规い鎺戝閸庡秹鏌涢弴銊ヤ簻闁哄們鍥ㄧ厱闁靛／鍐冪偤鏌嶈閸撱劑骞忛敓锟�
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
        //printf("Exz Ezx 濠电偞鍨堕幐鎾磻閹剧粯鈷戦悹鎭掑妼閺嬫垿鏌＄€ｎ亶鐓兼鐐茬箻閹粓鎳為妷锔筋仧闂備礁鎼崐鍫曞磹閺嶎偀鍋撳顒傜Ш闁哄被鍔戦幃銏ゅ川婵犲嫪绱曢梻浣哥秺椤ユ捇宕楀鈧顐﹀箻鐠囧弶顥濋梺闈涚墕濡顢旈崼鏇熲拺閻犳亽鍔岄弸鎴︽煛鐎ｎ亶鐓兼鐐茬箻閺屻劎鈧絽鐏氭鍕⒒娴ｈ姤纭堕柛鐘虫尵缁﹦缂撻敓锟�");
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Exz Ezx 濠电偞鍨堕幐鎾磻閹剧粯鈷戦悹鎭掑妼閺嬫垿鏌＄€ｎ亶鐓兼鐐茬箻閹粓鎳為妷锔筋仧闂備礁鎼崐鍫曞磹閺嶎偀鍋撳顒傜Ш闁哄被鍔戦幃銏ゅ川婵犲嫪绱曢梻浣哥秺椤ユ捇宕楀鈧顐﹀箻鐠囧弶顥濋梺闈涚墕濡顢旈崼鏇熲拺閻犳亽鍔岄弸鎴︽煛鐎ｎ亶鐓兼鐐茬箻閺屻劎鈧絽鐏氭鍕⒒娴ｈ姤纭堕柛锝忕畵楠炲繘鏁撻敓锟�?)  %%%%%%%%%%%%%%%%%%%%%%%%%%
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
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Eyx Exy 濠电偞鍨堕幐鎾磻閹剧粯鈷戦悹鎭掑妼閺嬫垿鏌＄€ｎ亶鐓兼鐐茬箻閹粓鎳為妷锔筋仧闂備礁鎼崐鍫曞磹閺嶎偀鍋撳顒傜Ш闁哄被鍔戦幃銏ゅ川婵犲嫪绱曢梻浣哥秺椤ユ捇宕楀鈧顐﹀箻鐠囧弶顥濋梺闈涚墕濡顢旈崼鏇熲拺閻犳亽鍔岄弸鎴︽煛鐎ｎ亶鐓兼鐐茬箻閺屻劎鈧絽鐏氭鍕⒒娴ｈ姤纭堕柛锝忕畵楠炲繘鏁撻敓锟�?)  %%%%%%%%%%%%%%%%%%%%%%%%%%
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
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Ezy Eyz 濠电偞鍨堕幐鎾磻閹剧粯鈷戦悹鎭掑妼閺嬫垿鏌＄€ｎ亶鐓兼鐐茬箻閹粓鎳為妷锔筋仧闂備礁鎼崐鍫曞磹閺嶎偀鍋撳顒傜Ш闁哄被鍔戦幃銏ゅ川婵犲嫪绱曢梻浣哥秺椤ユ捇宕楀鈧顐﹀箻鐠囧弶顥濋梺闈涚墕濡顢旈崼鏇熲拺閻犳亽鍔岄弸鎴︽煛鐎ｎ亶鐓兼鐐茬箻閺屻劎鈧絽鐏氭鍕⒒娴ｈ姤纭堕柛锝忕畵楠炲繘鏁撻敓锟�?)  %%%%%%%%%%%%%%%%%%%%%%%%%%
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
   // cout<<f2<<endl;
    if (!finput1.is_open()) {
        std::cerr << "Failed to open file: " << f1 << std::endl;
       
    }

    while(!finput1.eof())
    {
       // printf("闂佽瀛╅鏍窗閹烘纾婚柟鐐灱閺€鑺ャ亜閺冨倵鎷￠柛搴㈡尰閵囧嫰寮村顒€绁悗娈垮櫘閸撶喎鐣烽崼鏇熸櫆缂佹稑顑呮禒绡璶");
        complex<double> c1;
        if(finput1.fail()||k==m_len)
        {
           // printf("闂備浇宕垫慨鏉懨洪埡鍜佹晪鐟滄垿濡甸幇鏉跨倞闁靛ǹ鍎崇粣鐐烘煛婢跺﹦澧遍柛瀣槼椤ゅ倿姊绘担鐑樺殌闁宦板姂閹囨倷閸濆嫮鍔﹀銈嗗坊閸嬫挾鈧娲╅幏锟�==13\n");
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