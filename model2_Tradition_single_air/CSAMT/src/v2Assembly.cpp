#include "../include/main.h"
#include "../include/fmodel.h"
using namespace std;

v2Assembly::v2Assembly(Fmodel fmodel,double f)
{
    
    w=2*M_PI*f;
    mu0=(4e-7)*M_PI;
    m_len=fmodel.NZ+1;
    //printf("m_len=%d\n",m_len);
    //闂傚倸鍊风粈渚€骞夐敍鍕殰婵°倕鍟畷鏌ユ煕瀹€鈧崕鎴犵礊閺嶎厽鐓欓柣妤€鐗婄欢鑼磼閳ь剙鐣濋崟顒傚幐闂佸壊鍋嗛崰鏇犫偓纰夋嫹+1=25
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
    //闂傚倷绀侀幖顐λ囬鐐村亱濠电姴娲ょ粻浼存煙闂傚顦﹂柛姘秺閺岋綁鎮╅幓鎺斿弮婵犵數濮烽。钘壩ｉ崨鏉戠；闁逞屽墴閺屾稓鈧綆鍋呭畷宀勬煛瀹€瀣？濞寸媴濡囬幏鐘诲箵閹烘埈娼ュ┑锛勫亼閸婃牜鏁Δ鍐ㄥ灊閹肩补妾уΣ鍫ユ煟閵忕姴顥忛柡浣革躬閹嘲鈻庡▎鎴濆煂婵炲瓨绮撶粻鏍蓟閵娾晜鍋嗛柛灞剧☉椤忥拷
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
    MatAXPY(copyK3, 1, v0A.K4, UNKNOWN_NONZERO_PATTERN); //闂傚倸鍊烽悞锕€顪冮崸妤€鍌ㄩ柟闂寸閻ら箖鏌ｅΟ鍝勬闁绘梻鍘ч悞鍨亜閹烘垵顏柣鎾存礋閺岋繝宕橀妸銉㈠亾閸涘﹦澧＄紓鍌欒兌閾忓酣宕㈡ィ鍐ｂ偓鏃堟晸閿燂拷3=K3+K4
    MatMult(copyK3, Ep, P);// 闂傚倸鍊烽悞锕€顪冮崸妤€鍌ㄩ柟闂寸閻ら箖鏌ｅΟ鍝勬闁绘梻鍘ч悞鍨亜閹哄棗浜鹃梺瀹犳椤﹂潧顕ｉ鈧崺鈧い鎺戝閺咁亪姊绘担绛嬪殐闁哥姵顨堥崚鎺楀箻鐠囪尪鎽曢悗鍏夊亾闁告洦鍓氬▍銏ゆ⒑鐠恒劌娅愰柟鍑ゆ嫹//P=(K3+K4)*Ep;  
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
    
    //PetscPrintf(PETSC_COMM_WORLD, "P=(K3+K4)*Ep闂傚倷娴囬褍霉閻戣棄绠犻柟鎹愵嚙鐎氬銇勯幒鎴濐仼闁搞劌鍊归幈銊モ攽閸ャ劉鏋�");
    
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
        MatZeroRowsColumns(v2, count_row, row, 1.0, NULL, NULL);//闂傚倸鍊烽懗鍫曞磿閻㈢ǹ纾婚柟鎹愵嚙缁€澶愭倵閿濆骸澧插┑顔挎珪閵囧嫰骞掗幋婵冨亾閻戣姤鍊垮┑鐘叉处閻撴洟鏌ｉ幇顒傛憼闁崇粯鎮╅梻鍌氬€烽悞锕傛儑瑜版帒绀夌€光偓閳ь剟鍩€椤掍礁鍤柛銊ョ埣婵″瓨绗熼埀顒勭嵁鐎ｎ喗鍊风痪鐗埳戦悘鍐ㄢ攽閻樺灚鏆╁┑顔碱嚟閳ь剚鍑归崣鍐箖閻愵剚缍囬柕濠忕畱瀵潡姊洪幐搴ｇ畵闁硅櫕鍔楃划濠氭晲婢跺鎷洪梺鍛婄懃椤﹂亶鎯岄幒鏂哄亾鐟欏嫭绀冮柛鏃€鐟╅獮濠傗攽鐎ｎ偆鍔烽梺鎸庢磵閸嬫捇鏌￠崪浣稿箻缂佽鲸鎹囧畷鎺戭潩椤掍胶浜鹃梻浣告憸婵潧煤濠婂牆绠憸鐗堝笚閺呮悂鏌ㄩ悤鍌涘0
        PetscTime(&end_matv2);
       // PetscPrintf(curComm, "MatZeroRowsColumns_V1 time(MPI):%f\n",end_matv2-start_matv2);
        
    //PetscPrintf(curComm, "====闂備浇顕х€涒晠顢欓弽顓炵獥闁圭儤顨呯壕濠氭煙閻愵剚鐏遍柡鈧懞銉ｄ簻闁哄啫鍊甸幏锟犳煕鎼淬垹濮嶉柡宀嬬秮椤㈡棃宕ㄩ褎顥嬮梻浣烘嚀瀵爼骞愰幎钘夋瀬闁瑰墽绮弲鎼佹煥閻曞倹瀚�======\n");
   //闂傚倷绀侀幖顐λ囬锕€鐤鹃柍鍝勬噹閸屻劍绻涘顔荤盎闁兼瓕顫夐妵鍕箳閹存績鍋撹ぐ鎺戠獥闁糕剝绋掗悡銉︾節闂堟稒顥滄い蹇ｅ亰閺岋綁鎮㈤崫銉х杽闂佸搫鐭夌换婵嬪春閳ь剚銇勯幒鎴濐仾闁搞倕顑夐弻娑⑩€﹂幋婵呯凹闂佷紮缍囬幏锟�
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
    //闂傚倸鍊风粈渚€宕ョ€ｎ喖纾块柟鎯版鎼村﹪鏌ら懝鎵牚濞存粌缍婇弻娑㈠Ψ椤旂厧顫╃紒鎯у⒔缁垳鎹㈠┑鍥╃瘈闁稿本鍑规导宀勬⒑閸涘﹥灏甸梺钘夘敍P闂傚倷绀侀幖顐λ囬锕€鐤鹃柍鍝勬噹閸屻劍绻涘顔荤盎闁兼瓕顫夐妵鍕箳閹存繍浠鹃梺绋块缁夊綊寮诲☉銏犲嵆闁靛ǹ鍎扮花濠氭⒑閸濆嫬鈧垿宕堕妸褍骞堥梻渚€娼ч敍蹇涘礃閸愵亜浠遍柡灞界Х椤т線鏌涢幘鍗炲妤犵偛绻橀弫鎾绘晸閿燂拷
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
        //printf("Exz Ezx 濠电姷鏁搁崑鐐哄垂閸洖绠伴柟缁㈠枛绾惧鏌熼崜褏甯涢柍閿嬪灴閹綊骞侀幒鎴濐瀳闂佸搫顑嗛崹鍧楀蓟閿涘嫧鍋撻敐搴濇喚闁绘挸鍚嬮〃銉╂倷閼碱剛顔掗梺璇″枟缁捇骞婇悙鍝勎ㄩ柨鏃傜摂娴犙囨⒒閸屾瑧顦﹂柟纰卞亰瀹曟劙宕奸弴鐐碉紮闂佸搫绋侀崑鈧柛瀣尭椤繈顢楅崒婧炪劑姊洪崫鍕潶闁告梹鍨块獮鍐閵堝懎绐涙繝鐢靛Т鐎氼亞妲愰弴銏♀拻濞达絽鎽滅粔鐑樸亜閵夛附宕岀€规洘顨呴～婊堝焵椤掆偓椤曪綁顢曢敃鈧粻濠氭偣閸パ冪骇妞ゃ儲绻堝娲濞戞艾顣哄┑鈽嗗亝椤ㄥ﹪銆侀弮鍫濋唶闁哄洨鍟块幏娲煟閻樺厖鑸柛鏂跨焸瀵悂骞嬮敂鐣屽幐闁诲函缍嗘禍鍫曟偂閸忕⒈娈介柣鎰皺缁犲鏌＄仦璇插闁逞屽墾缁蹭粙鎮樺顓ф闁告洦鍨遍埛鎺懨归敐鍫澬撶痪顓炵埣閺屾盯鎮╅搹顐㈡殫缂備緡鍠栭敃锔剧磽閹剧粯鏅搁柨鐕傛嫹");
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Exz Ezx 濠电姷鏁搁崑鐐哄垂閸洖绠伴柟缁㈠枛绾惧鏌熼崜褏甯涢柍閿嬪灴閹綊骞侀幒鎴濐瀳闂佸搫顑嗛崹鍧楀蓟閿涘嫧鍋撻敐搴濇喚闁绘挸鍚嬮〃銉╂倷閼碱剛顔掗梺璇″枟缁捇骞婇悙鍝勎ㄩ柨鏃傜摂娴犙囨⒒閸屾瑧顦﹂柟纰卞亰瀹曟劙宕奸弴鐐碉紮闂佸搫绋侀崑鈧柛瀣尭椤繈顢楅崒婧炪劑姊洪崫鍕潶闁告梹鍨块獮鍐閵堝懎绐涙繝鐢靛Т鐎氼亞妲愰弴銏♀拻濞达絽鎽滅粔鐑樸亜閵夛附宕岀€规洘顨呴～婊堝焵椤掆偓椤曪綁顢曢敃鈧粻濠氭偣閸パ冪骇妞ゃ儲绻堝娲濞戞艾顣哄┑鈽嗗亝椤ㄥ﹪銆侀弮鍫濋唶闁哄洨鍟块幏娲煟閻樺厖鑸柛鏂跨焸瀵悂骞嬮敂鐣屽幐闁诲函缍嗘禍鍫曟偂閸忕⒈娈介柣鎰皺缁犲鏌＄仦璇插闁逞屽墾缁蹭粙鎮樺顓ф闁告洦鍨遍埛鎺懨归敐鍫澬撶痪顓炵埣閺屾盯鏁愯箛鏇犳殼濡ょ姷鍋涚换姗€寮幘缁樻櫢闁跨噦鎷�?)  %%%%%%%%%%%%%%%%%%%%%%%%%%
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
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Eyx Exy 濠电姷鏁搁崑鐐哄垂閸洖绠伴柟缁㈠枛绾惧鏌熼崜褏甯涢柍閿嬪灴閹綊骞侀幒鎴濐瀳闂佸搫顑嗛崹鍧楀蓟閿涘嫧鍋撻敐搴濇喚闁绘挸鍚嬮〃銉╂倷閼碱剛顔掗梺璇″枟缁捇骞婇悙鍝勎ㄩ柨鏃傜摂娴犙囨⒒閸屾瑧顦﹂柟纰卞亰瀹曟劙宕奸弴鐐碉紮闂佸搫绋侀崑鈧柛瀣尭椤繈顢楅崒婧炪劑姊洪崫鍕潶闁告梹鍨块獮鍐閵堝懎绐涙繝鐢靛Т鐎氼亞妲愰弴銏♀拻濞达絽鎽滅粔鐑樸亜閵夛附宕岀€规洘顨呴～婊堝焵椤掆偓椤曪綁顢曢敃鈧粻濠氭偣閸パ冪骇妞ゃ儲绻堝娲濞戞艾顣哄┑鈽嗗亝椤ㄥ﹪銆侀弮鍫濋唶闁哄洨鍟块幏娲煟閻樺厖鑸柛鏂跨焸瀵悂骞嬮敂鐣屽幐闁诲函缍嗘禍鍫曟偂閸忕⒈娈介柣鎰皺缁犲鏌＄仦璇插闁逞屽墾缁蹭粙鎮樺顓ф闁告洦鍨遍埛鎺懨归敐鍫澬撶痪顓炵埣閺屾盯鏁愯箛鏇犳殼濡ょ姷鍋涚换姗€寮幘缁樻櫢闁跨噦鎷�?)  %%%%%%%%%%%%%%%%%%%%%%%%%%
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
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Ezy Eyz 濠电姷鏁搁崑鐐哄垂閸洖绠伴柟缁㈠枛绾惧鏌熼崜褏甯涢柍閿嬪灴閹綊骞侀幒鎴濐瀳闂佸搫顑嗛崹鍧楀蓟閿涘嫧鍋撻敐搴濇喚闁绘挸鍚嬮〃銉╂倷閼碱剛顔掗梺璇″枟缁捇骞婇悙鍝勎ㄩ柨鏃傜摂娴犙囨⒒閸屾瑧顦﹂柟纰卞亰瀹曟劙宕奸弴鐐碉紮闂佸搫绋侀崑鈧柛瀣尭椤繈顢楅崒婧炪劑姊洪崫鍕潶闁告梹鍨块獮鍐閵堝懎绐涙繝鐢靛Т鐎氼亞妲愰弴銏♀拻濞达絽鎽滅粔鐑樸亜閵夛附宕岀€规洘顨呴～婊堝焵椤掆偓椤曪綁顢曢敃鈧粻濠氭偣閸パ冪骇妞ゃ儲绻堝娲濞戞艾顣哄┑鈽嗗亝椤ㄥ﹪銆侀弮鍫濋唶闁哄洨鍟块幏娲煟閻樺厖鑸柛鏂跨焸瀵悂骞嬮敂鐣屽幐闁诲函缍嗘禍鍫曟偂閸忕⒈娈介柣鎰皺缁犲鏌＄仦璇插闁逞屽墾缁蹭粙鎮樺顓ф闁告洦鍨遍埛鎺懨归敐鍫澬撶痪顓炵埣閺屾盯鏁愯箛鏇犳殼濡ょ姷鍋涚换姗€寮幘缁樻櫢闁跨噦鎷�?)  %%%%%%%%%%%%%%%%%%%%%%%%%%
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
       // printf("闂傚倷娴囬褏鈧稈鏅犻、娆撳冀椤撶偟鐛ラ梺鍦劋椤ㄥ懐澹曟繝姘厵闁绘劦鍓氶悘閬嶆煛閳ь剟鎳為妷锝勭盎闂佸搫鍟崐鐢稿箯閿熺姵鐓曢幖娣灪鐏忎即鏌曢崶褍顏€殿喗娼欓～婵嬵敆閳ь剛绮婇鈧幃妤€鈻撻崹顔界彯闂佸憡鎸鹃崰搴ㄦ偩閻戣棄閱囬柡鍥╁枑濞呭棛绱撴担鍦槈妞ゆ垵鎳忕粋鎺斿姬閻狅拷");
        complex<double> c1;
        if(finput1.fail()||k==m_len)
        {
           // printf("闂傚倸鍊峰ù鍥х暦閸偅鍙忛柡澶嬪殮濞差亜鐓涢柛婊€鐒﹂弲顏堟偡濠婂嫬鐏村┑锛勬暬楠炲洭寮剁捄銊モ偓鐐烘⒑闂堟胆褰掑磿瀹曞洨鐜婚柣鎰劋閻撴稑顭跨捄鐚村姛濠⒀囦憾閺屾稓鈧綆浜跺Σ鍏笺亜閵堝懎鈧灝顫忕紒妯诲闁绘垶锚濞堝矂姊虹€癸附婢樻慨鍌炴煙椤斻劌娲﹂崐鐑芥煕濠靛棗顏柛鏃撶畱椤啴濡堕崱妤€娼戦梺绋款儐閹搁箖鍩€椤掑喚娼愭繛娴嬫櫊楠炲繘鏁撻敓锟�==13\n");
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