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
    //闂傚倸鍊搁崐鎼佸磹妞嬪海鐭嗗〒姘ｅ亾妤犵偛顦甸弫宥夊礋椤掍焦顔囨繝寰锋澘鈧洟宕姘辨殾闁哄被鍎查悡鏇犫偓鍏夊亾闁逞屽墴瀹曟洟骞嬮悩鐢殿槸闂佸搫绋侀崢浠嬫偂濞嗘挻鐓熸俊銈傚亾闁绘锕﹀▎銏ゆ嚑椤掑倻锛滈梺閫炲苯澧柣锝嗙箞瀹曠喖顢楅崒姘闂傚倷绀佹竟濠囧磻閸℃稑绐楅柡鍥╁Л閸嬫挾鍠婃径瀣伓+1=25
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
    int NN=fmodel.NN;
    
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
        EpyB(i)=linshi(i-EpyBstart)/linshi(Nair-NPML_Z);//ExA(7:19,1)=linshi(:,1)/linshi(3); 
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
    //闂傚倸鍊搁崐椋庣矆娓氣偓楠炴牠顢曚綅閸ヮ剦鏁嶉柣鎰綑娴滆鲸绻濋悽闈浶㈡繛灞傚€楃划缁樺鐎涙鍘甸梻鍌氬€搁顓⑺囬敃鍌涚厱婵﹩鍓涚粔娲煛瀹€瀣М闁诡喒鏅犻獮鎾诲箳閺傚灝绫嶆繝鐢靛Х閺佸憡鎱ㄩ悜濮愨偓鍌炴寠婢光晪缍佸畷銊╁级閹寸媴绱遍梻渚€鈧偛鑻晶鎾煛鐏炲墽銆掗柍褜鍓ㄧ紞鍡涘磻閸涱厾鏆︾€光偓閸曨剛鍘搁悗鍏夊亾閻庯綆鍓涢敍鐔哥箾鐎电ǹ顎撳┑鈥虫喘楠炲繘鎮╃拠鑼唽闂佸湱鍎ら崺鍫濐焽閵夈儮鏀介柨娑樺娴滃ジ鏌涙繝鍐⒈闁轰緡鍠栬灃闁告劑鍔岄悘濠囨煙閼测晞藟婵＄嫏鍏撅綁宕奸妷锔惧幗闂侀潧绻掓慨鎾夎箛娑欑厸濞达綁娼婚煬顒勬煙椤旂厧妲婚柍璇查閳诲酣骞嬪┑鍡欏帓婵犵數鍋涢悺銊у垝閹惧墎涓嶉柡宥庡幖閽冪喖鏌曟繛鐐珕闁稿妫濋弻娑氫沪閸撗€妲堝銈呴獜閹凤拷
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
    MatAXPY(copyK3, 1, v0A.K4, UNKNOWN_NONZERO_PATTERN); //闂傚倸鍊搁崐鎼佸磹閻戣姤鍊块柨鏇楀亾妞ゎ亜鍟村畷绋课旈埀顒勫磼閵娾晜鐓欓梻鍌氼嚟椤︼箓鏌ｉ妶澶岀暫闁哄矉绲借灒闁告繂瀚ˉ婵嬫⒑缂佹ɑ鈷掗柛妯犲洦鍊块柛顭戝亖娴滄粓鏌熼悜妯虹仴妞ゅ浚浜弻锝夊箻鐎涙顦伴梺鍝勭焿缁绘繂鐣峰鈧俊鎼佸Ψ閵忕姳澹曢梺鍛婄☉閿曪附鏅堕敍鍕＝闁稿本鐟ㄩ崗宀勬煣韫囨捇鍙勭€规洏鍨洪妶锝夊礃閿濆倸浜鹃柡鍐ㄧ墛閺呮悂鏌ㄩ悤鍌涘3=K3+K4
    MatMult(copyK3, Ep, P);// 闂傚倸鍊搁崐鎼佸磹閻戣姤鍊块柨鏇楀亾妞ゎ亜鍟村畷绋课旈埀顒勫磼閵娾晜鐓欓梻鍌氼嚟椤︼箓鏌ｉ妶澶岀暫闁哄矉绲借灒闁告繂瀚ˉ婵嬫⒑缂佹ɑ鈷掗柛妯犲洦鍊块柛顭戝亖娴滄粓鏌熼崫鍕棞濞存粓绠栧铏光偓鍦У椤ュ銇勯敃鍌涙锭妞ゆ洩缍侀、姘跺焵椤掑嫬鏄ラ柍褜鍓氶妵鍕箳閹存繍浠鹃梺鍝勬媼娴滎亜顫忕紒妯诲缂佹稑顑呭▓鎰版⒑閸濄儱校妞ゃ劌鐗撳畷姘跺箳濡も偓缁犲鎮归崶顏勮敿闁硅姤娲熼幃妤呭礂婢跺﹣澹曢梻浣告啞濞诧箓宕滃顑芥瀺闁靛繈鍊栭埛鎴︽偣閹帒濡兼繛鍛姍閺岀喖宕欓妶鍡楊伓//P=(K3+K4)*Ep;  
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
    
    //PetscPrintf(PETSC_COMM_WORLD, "P=(K3+K4)*Ep闂傚倸鍊搁崐宄懊归崶顒夋晪鐟滃秹婀侀梺缁樺灱濡嫮绮婚悩缁樼厵闁硅鍔﹂崵娆撴倵濮橆剦妲归柕鍥у楠炴帡骞嬪┑鎰偅闂備焦鎮堕崝宀勫磹瑜版帒绠為柕濞垮劗閺€浠嬫煕閵夛絽濡烽柡瀣舵嫹");
    
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
        MatZeroRowsColumns(v2, count_row, row, 1.0, NULL, NULL);//闂傚倸鍊搁崐鎼佸磹閻戣姤鍤勯柛顐ｆ礀绾惧潡鏌ｉ姀鈶跺湱澹曟繝姘厵闁硅鍔﹂崵娆戠磼閳ь剚寰勯幇顓炩偓鐢告煥濠靛棝顎楀褎褰冮埞鎴︻敊閹稿海褰ч梺闈涙搐鐎氫即鐛幒妤€绠ｆ繝闈涘暙娴滈箖鏌ｉ幋锝呅撻柛濠傜仢閳规垿鎮╅崣澶婎槱闂佺粯鎸诲ú鐔煎蓟閿濆绠涙い鎺戝€归幉濂告⒑瀹曞洨甯涢柟顔规櫊濮婂宕掑顑藉亾閻戣姤鍊块柨鏇炲€归崕鎴犳喐閻楀牆绗掔紒鈧径灞稿亾閸忓浜鹃梺閫炲苯澧撮柛鈹惧亾濡炪倖甯婄粈渚€宕甸鍕厱闁靛ǹ鍎抽崺锝呪攽閳ュ磭鎽犵紒妤冨枛閸┾偓妞ゆ帒瀚畵渚€鎮楅敐搴℃灍闁稿﹪顥撻惀顏堟偋閸╄櫕鍨块幃姗€宕橀妸銏℃杸闂佺粯锚閻忔岸寮抽埡浣叉斀妞ゆ梻鈷堥崵鐔兼煃瑜滈崜姘跺礄瑜版帒鐭楅柛鎰靛枛缁犳牠鏌ｉ幇闈涘缂傚秴娲弻鏇熺節韫囨洜鏆犻悗娈垮枛濞尖€愁潖濞差亜绠伴幖杈剧悼閻ｇ敻姊虹涵鍛彧闁告梹顨堥崚鎺撶節濮橆厽娅滄繝銏ｆ硾椤戝倿骞忓ú顏呪拺闁告稑锕﹂幊鍐┿亜閿曞倷鎲鹃柟顖氱焸楠炴帡寮崫鍕闁荤喐鐟ョ€氼厾绮堥崘顔界厱闁哄啠鍋撻柣鐔叉櫊閻涱喗绻濋崒妤佹杸闁诲函缍嗛崑鍡涘礉閻戣姤鈷戦柟绋垮绾剧敻鏌涚€ｎ偅宕岄柡宀嬬節瀹曨亝鎷呯粙璺啋缂傚倷娴囨ご鎼佸箲閸パ呮殾闁圭儤鍩堝鈺傘亜閹哄秷鍏屽ù婊堢畺濮婄粯鎷呴崨濠冨枑婵犳鍣ｅ褏鍙呭┑鐘诧工閻楀棛绮婚婊勫枑闁绘鐗嗙粭姘舵煛閸涱喗鍊愰柡灞诲姂閹倝宕掑☉姗嗕紦0
        PetscTime(&end_matv2);
        //PetscPrintf(curComm, "MatZeroRowsColumns_V1 time(MPI):%f\n",end_matv2-start_matv2);
        
    //PetscPrintf(curComm, "==========\n");
   //闂傚倸鍊搁崐椋庣矆娓氣偓楠炴牠顢曚綅閸ヮ剦鏁冮柨鏇楀亾闁汇倝绠栭弻宥夊传閸曨剙娅ｉ梺绋胯閸斿秶鎹㈠☉姗嗗晠妞ゆ棁宕甸惄搴ㄦ⒑閸忚偐鎽冩い顐㈩樀婵＄敻宕熼姘鳖唺闂佺懓鐡ㄧ缓楣冨磻閹惧箍浜归柟鐑樺灩閻涖儵姊虹化鏇炲⒉缂佸甯￠幃锟犲Ψ閿斿墽鐦堥梻鍌氱墛缁嬫帡藟濠婂嫨浜滈煫鍥风到娴滀即鏌″畝瀣М闁诡喓鍨藉畷顐﹀Ψ瑜忛弶浠嬫⒒娴ｅ憡鎯堥柣顓烆槺閹广垹鈹戠€ｎ亝妲梺閫炲苯澧柕鍥у楠炴帡骞嬪┑鎰棯闂備焦鎮堕崐鏇㈩敄婢舵劕钃熸繛鎴旀噰閳ь剨绠撻獮瀣攽閸涱垰鍤梻鍌欐缁鳖喚绱為崶顒€绠柨鐕傛嫹
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
              //Exx2[i-1][j-1]=((Ex2(i,j+1,Nair+NN+1)+Ex2(i,j,Nair+NN+1))/2.d0+(Ex2(i,j+1,Nair+NN+2)+Ex2(i,j,Nair+NN+2))/2.d0)/2.d0;
                Exx2[i-1][j-1]=((Ex2[i-1][j][Nair+NN]+Ex2[i-1][j-1][Nair+NN])/2.0+(Ex2[i-1][j][Nair+NN+1]+Ex2[i-1][j-1][Nair+NN+1])/2.0)/2.0;

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
              //Eyy2(i,j)=((Ey2(i,j,Nair+NN+1)+Ey2(i+1,j,Nair+NN+1))/2.d0+(Ey2(i,j,Nair+NN+2)+Ey2(i+1,j,Nair+NN+2))/2.d0)/2.d0;
                Eyy2[i-1][j-1]=((Ey2[i-1][j-1][Nair+NN]+Ey2[i][j-1][Nair+NN])/2.0+(Ey2[i-1][j-1][Nair+NN+1]+Ey2[i][j-1][Nair+NN+1])/2.0)/2.0;
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
              //Ezz2(i,j)=((Ez2(i,j,Nair+NN+1)+Ez2(i,j+1,Nair+NN+1))/2.d0+(Ez2(i+1,j,Nair+NN+1)+Ez2(i+1,j+1,Nair+NN+1))/2.d0)/2.d0;
                Ezz2[i-1][j-1]=((Ez2[i-1][j-1][Nair+NN]+Ez2[i-1][j][Nair+NN])/2.0+(Ez2[i][j-1][Nair+NN]+Ez2[i][j][Nair+NN])/2.0)/2.0;
                //printf("Ezz2(%d,%d)=%e,%e\n",i,j,Ezz2[i-1][j-1].real(),Ezz2[i-1][j-1].imag());
            }
        }
        //printf("Exz Ezx 濠电姷鏁告慨鐑藉极閹间礁纾婚柣鎰惈閸ㄥ倿鏌涢锝嗙缂佺姳鍗抽弻鐔虹磼閵忕姵鐏堢紒鐐劤椤兘寮婚悢鐓庣鐟滃繒鏁☉銏＄厓闂佸灝顑呴悘鎾煙椤旇偐绉烘鐐扮窔楠炴帡骞嬪┑鎰偓鎶芥⒒娴ｅ憡鎯堟い鎴濇瀹曞綊宕稿Δ鈧拑鐔兼煥濞戞ê顏ら柛瀣崌閺佹劖鎯斿┑鍥ф灓闂備胶绮幐鎼佸触鐎ｎ喓鈧啴濡烽埡鍌氣偓鐑芥煠绾板崬澧版い鏃€甯″铏规嫚閳ヨ櫕鐏嶇紓渚囧枟閹瑰洭鐛繝鍥ㄥ€烽柛婵嗗珋閵娾晜鐓ラ柡鍐ㄥ€婚幗鍌毭归悩娆忔处閳锋帡鏌涚仦鍓ф噮妞わ讣绠撻弻鐔哄枈閸楃偘鍠婇悗瑙勬礃閸旀瑥鐣锋總绋垮嵆闁绘劗顣槐顕€姊绘担鍛婃儓缂佸绶氬畷鎴﹀焵椤掑嫭鐓曢悗锝庡亝鐏忣厽銇勯锝囩疄妞ゃ垺顨婂畷鎺戔攦閻愵亜濮傛慨濠冩そ瀹曨偊宕熼鍛晧闂備礁鎲″褰掑垂閸ф宓侀柛鎰靛枛椤懘鏌曢崼婵囧櫣缂佹劖绋掔换婵嬫偨闂堟刀銏ゆ倵濮橀棿绨芥俊鍙夊姍瀵挳濮€閳锯偓閹风粯绻涙潏鍓у埌闁硅姤绮庣划鏃堟倻濡晲绨婚梺闈涱檧闂勫嫬鐣峰畝鈧埀顒冾潐濞叉﹢銆冮崨杈剧稏婵犲﹤鐗嗛悞鍨亜閹哄棗浜惧銈嗘穿缂嶄線銆侀弴銏℃櫇闁逞屽墰缁粯绻濆顓炰化闂佹悶鍎滈崘顏堢崜婵＄偑鍊戦崕鑼崲閸繍娼栨繛宸簼椤ュ牊绻涢幋鐐跺妞わ絽鎼埞鎴﹀煡閸℃ぞ绨煎銈冨妼閿曨亪濡存笟鈧顕€宕煎┑瀣暪闂備礁鎼ú銊╁疮閸ф绠繛宸簼閻撶喖鏌ｅΟ鍝勫笭闁煎壊浜弻娑㈠棘鐠恒劎鍔悗娈垮枟閹倿鐛€ｎ喗鏅滈柣锝呰嫰楠炴劙姊虹拠鎻掑毐缂傚秴妫欑粋宥夊醇閺囩喎浜楅梺绋跨箳閳峰牆鈻撴禒瀣厽闁归偊鍨伴惃铏圭磼閻樺樊鐓奸柡宀嬬磿娴狅妇鎷犻幓鎺濈€烽梻渚€鈧偛鑻晶鍓х磼闊厾鐭欓柟顔矫～婵嬵敇瑜庨鏃堟⒑閸涘﹥澶勯柛銊╀憾閸╂盯骞掗幊銊ョ秺閺佹劙宕煎顒佹尵閻ヮ亪顢橀悙闈涚厽闂佸搫鏈惄顖炲箖閳哄懏鎯炴い鎰╁灪濞堫偆绱撻崒娆戣窗闁哥姵鐗犻弫鍐晜閸撗咃紳闂佺懓澧界划顖炲疾閹间焦鐓ラ柣鏇炲€圭€氾拷");
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Exz Ezx 濠电姷鏁告慨鐑藉极閹间礁纾婚柣鎰惈閸ㄥ倿鏌涢锝嗙缂佺姳鍗抽弻鐔虹磼閵忕姵鐏堢紒鐐劤椤兘寮婚悢鐓庣鐟滃繒鏁☉銏＄厓闂佸灝顑呴悘鎾煙椤旇偐绉烘鐐扮窔楠炴帡骞嬪┑鎰偓鎶芥⒒娴ｅ憡鎯堟い鎴濇瀹曞綊宕稿Δ鈧拑鐔兼煥濞戞ê顏ら柛瀣崌閺佹劖鎯斿┑鍥ф灓闂備胶绮幐鎼佸触鐎ｎ喓鈧啴濡烽埡鍌氣偓鐑芥煠绾板崬澧版い鏃€甯″铏规嫚閳ヨ櫕鐏嶇紓渚囧枟閹瑰洭鐛繝鍥ㄥ€烽柛婵嗗珋閵娾晜鐓ラ柡鍐ㄥ€婚幗鍌毭归悩娆忔处閳锋帡鏌涚仦鍓ф噮妞わ讣绠撻弻鐔哄枈閸楃偘鍠婇悗瑙勬礃閸旀瑥鐣锋總绋垮嵆闁绘劗顣槐顕€姊绘担鍛婃儓缂佸绶氬畷鎴﹀焵椤掑嫭鐓曢悗锝庡亝鐏忣厽銇勯锝囩疄妞ゃ垺顨婂畷鎺戔攦閻愵亜濮傛慨濠冩そ瀹曨偊宕熼鍛晧闂備礁鎲″褰掑垂閸ф宓侀柛鎰靛枛椤懘鏌曢崼婵囧櫣缂佹劖绋掔换婵嬫偨闂堟刀銏ゆ倵濮橀棿绨芥俊鍙夊姍瀵挳濮€閳锯偓閹风粯绻涙潏鍓у埌闁硅姤绮庣划鏃堟倻濡晲绨婚梺闈涱檧闂勫嫬鐣峰畝鈧埀顒冾潐濞叉﹢銆冮崨杈剧稏婵犲﹤鐗嗛悞鍨亜閹哄棗浜惧銈嗘穿缂嶄線銆侀弴銏℃櫇闁逞屽墰缁粯绻濆顓炰化闂佹悶鍎滈崘顏堢崜婵＄偑鍊戦崕鑼崲閸繍娼栨繛宸簼椤ュ牊绻涢幋鐐跺妞わ絽鎼埞鎴﹀煡閸℃ぞ绨煎銈冨妼閿曨亪濡存笟鈧顕€宕煎┑瀣暪闂備礁鎼ú銊╁疮閸ф绠繛宸簼閻撶喖鏌ｅΟ鍝勫笭闁煎壊浜弻娑㈠棘鐠恒劎鍔悗娈垮枟閹倿鐛€ｎ喗鏅滈柣锝呰嫰楠炴劙姊虹拠鎻掑毐缂傚秴妫欑粋宥夊醇閺囩喎浜楅梺绋跨箳閳峰牆鈻撴禒瀣厽闁归偊鍨伴惃铏圭磼閻樺樊鐓奸柡宀嬬磿娴狅妇鎷犻幓鎺濈€烽梻渚€鈧偛鑻晶鍓х磼闊厾鐭欓柟顔矫～婵嬵敇瑜庨鏃堟⒑閸涘﹥澶勯柛銊╀憾閸╂盯骞掗幊銊ョ秺閺佹劙宕煎顒佹尵閻ヮ亪顢橀悙闈涚厽闂佸搫鏈惄顖炲极閹邦垳鐤€闁哄洨濮靛▓鍏间繆閵堝洤啸闁稿绋撻幑銏狀潨閳ь剙顕ｉ锕€绠荤紓浣股戝▍銏ゆ⒑鐠恒劌娅愰柟鍑ゆ嫹?)  %%%%%%%%%%%%%%%%%%%%%%%%%%
        vector<vector<complex<double>>> Ezx2(NX, vector<complex<double>>(NY));
        vector<vector<complex<double>>> Exz2(NX, vector<complex<double>>(NY));

        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
                Ezx2[i-1][j-1]=(Ez2[i][j-1][Nair+NN]-Ez2[i-1][j-1][Nair+NN]+Ez2[i][j][Nair+NN]-Ez2[i-1][j][Nair+NN])/2.0/fmodel.A_X[i-1];
                Exz2[i-1][j-1]=(Ex2[i-1][j-1][Nair+NN+1]-Ex2[i-1][j-1][Nair+NN]+Ex2[i-1][j][Nair+NN+1]-Ex2[i-1][j][Nair+NN])/2.0/fmodel.C_Z[Nair+NN];
            }
        }
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Eyx Exy 濠电姷鏁告慨鐑藉极閹间礁纾婚柣鎰惈閸ㄥ倿鏌涢锝嗙缂佺姳鍗抽弻鐔虹磼閵忕姵鐏堢紒鐐劤椤兘寮婚悢鐓庣鐟滃繒鏁☉銏＄厓闂佸灝顑呴悘鎾煙椤旇偐绉烘鐐扮窔楠炴帡骞嬪┑鎰偓鎶芥⒒娴ｅ憡鎯堟い鎴濇瀹曞綊宕稿Δ鈧拑鐔兼煥濞戞ê顏ら柛瀣崌閺佹劖鎯斿┑鍥ф灓闂備胶绮幐鎼佸触鐎ｎ喓鈧啴濡烽埡鍌氣偓鐑芥煠绾板崬澧版い鏃€甯″铏规嫚閳ヨ櫕鐏嶇紓渚囧枟閹瑰洭鐛繝鍥ㄥ€烽柛婵嗗珋閵娾晜鐓ラ柡鍐ㄥ€婚幗鍌毭归悩娆忔处閳锋帡鏌涚仦鍓ф噮妞わ讣绠撻弻鐔哄枈閸楃偘鍠婇悗瑙勬礃閸旀瑥鐣锋總绋垮嵆闁绘劗顣槐顕€姊绘担鍛婃儓缂佸绶氬畷鎴﹀焵椤掑嫭鐓曢悗锝庡亝鐏忣厽銇勯锝囩疄妞ゃ垺顨婂畷鎺戔攦閻愵亜濮傛慨濠冩そ瀹曨偊宕熼鍛晧闂備礁鎲″褰掑垂閸ф宓侀柛鎰靛枛椤懘鏌曢崼婵囧櫣缂佹劖绋掔换婵嬫偨闂堟刀銏ゆ倵濮橀棿绨芥俊鍙夊姍瀵挳濮€閳锯偓閹风粯绻涙潏鍓у埌闁硅姤绮庣划鏃堟倻濡晲绨婚梺闈涱檧闂勫嫬鐣峰畝鈧埀顒冾潐濞叉﹢銆冮崨杈剧稏婵犲﹤鐗嗛悞鍨亜閹哄棗浜惧銈嗘穿缂嶄線銆侀弴銏℃櫇闁逞屽墰缁粯绻濆顓炰化闂佹悶鍎滈崘顏堢崜婵＄偑鍊戦崕鑼崲閸繍娼栨繛宸簼椤ュ牊绻涢幋鐐跺妞わ絽鎼埞鎴﹀煡閸℃ぞ绨煎銈冨妼閿曨亪濡存笟鈧顕€宕煎┑瀣暪闂備礁鎼ú銊╁疮閸ф绠繛宸簼閻撶喖鏌ｅΟ鍝勫笭闁煎壊浜弻娑㈠棘鐠恒劎鍔悗娈垮枟閹倿鐛€ｎ喗鏅滈柣锝呰嫰楠炴劙姊虹拠鎻掑毐缂傚秴妫欑粋宥夊醇閺囩喎浜楅梺绋跨箳閳峰牆鈻撴禒瀣厽闁归偊鍨伴惃铏圭磼閻樺樊鐓奸柡宀嬬磿娴狅妇鎷犻幓鎺濈€烽梻渚€鈧偛鑻晶鍓х磼闊厾鐭欓柟顔矫～婵嬵敇瑜庨鏃堟⒑閸涘﹥澶勯柛銊╀憾閸╂盯骞掗幊銊ョ秺閺佹劙宕煎顒佹尵閻ヮ亪顢橀悙闈涚厽闂佸搫鏈惄顖炲极閹邦垳鐤€闁哄洨濮靛▓鍏间繆閵堝洤啸闁稿绋撻幑銏狀潨閳ь剙顕ｉ锕€绠荤紓浣股戝▍銏ゆ⒑鐠恒劌娅愰柟鍑ゆ嫹?)  %%%%%%%%%%%%%%%%%%%%%%%%%%
        vector<vector<complex<double>>> Eyx2(NX, vector<complex<double>>(NY));
        vector<vector<complex<double>>> Exy2(NX, vector<complex<double>>(NY));

        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
                Eyx2[i-1][j-1]=(Ey2[i][j-1][Nair+NN]-Ey2[i-1][j-1][Nair+NN]+Ey2[i][j-1][Nair+NN+1]-Ey2[i-1][j-1][Nair+NN+1])/2.0/fmodel.A_X[i-1];
                Exy2[i-1][j-1]=(Ex2[i-1][j][Nair+NN]-Ex2[i-1][j-1][Nair+NN]+Ex2[i-1][j][Nair+NN+1]-Ex2[i-1][j-1][Nair+NN+1])/2.0/fmodel.B_Y[j-1];
            }
        }
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%  (Ezy Eyz 濠电姷鏁告慨鐑藉极閹间礁纾婚柣鎰惈閸ㄥ倿鏌涢锝嗙缂佺姳鍗抽弻鐔虹磼閵忕姵鐏堢紒鐐劤椤兘寮婚悢鐓庣鐟滃繒鏁☉銏＄厓闂佸灝顑呴悘鎾煙椤旇偐绉烘鐐扮窔楠炴帡骞嬪┑鎰偓鎶芥⒒娴ｅ憡鎯堟い鎴濇瀹曞綊宕稿Δ鈧拑鐔兼煥濞戞ê顏ら柛瀣崌閺佹劖鎯斿┑鍥ф灓闂備胶绮幐鎼佸触鐎ｎ喓鈧啴濡烽埡鍌氣偓鐑芥煠绾板崬澧版い鏃€甯″铏规嫚閳ヨ櫕鐏嶇紓渚囧枟閹瑰洭鐛繝鍥ㄥ€烽柛婵嗗珋閵娾晜鐓ラ柡鍐ㄥ€婚幗鍌毭归悩娆忔处閳锋帡鏌涚仦鍓ф噮妞わ讣绠撻弻鐔哄枈閸楃偘鍠婇悗瑙勬礃閸旀瑥鐣锋總绋垮嵆闁绘劗顣槐顕€姊绘担鍛婃儓缂佸绶氬畷鎴﹀焵椤掑嫭鐓曢悗锝庡亝鐏忣厽銇勯锝囩疄妞ゃ垺顨婂畷鎺戔攦閻愵亜濮傛慨濠冩そ瀹曨偊宕熼鍛晧闂備礁鎲″褰掑垂閸ф宓侀柛鎰靛枛椤懘鏌曢崼婵囧櫣缂佹劖绋掔换婵嬫偨闂堟刀銏ゆ倵濮橀棿绨芥俊鍙夊姍瀵挳濮€閳锯偓閹风粯绻涙潏鍓у埌闁硅姤绮庣划鏃堟倻濡晲绨婚梺闈涱檧闂勫嫬鐣峰畝鈧埀顒冾潐濞叉﹢銆冮崨杈剧稏婵犲﹤鐗嗛悞鍨亜閹哄棗浜惧銈嗘穿缂嶄線銆侀弴銏℃櫇闁逞屽墰缁粯绻濆顓炰化闂佹悶鍎滈崘顏堢崜婵＄偑鍊戦崕鑼崲閸繍娼栨繛宸簼椤ュ牊绻涢幋鐐跺妞わ絽鎼埞鎴﹀煡閸℃ぞ绨煎銈冨妼閿曨亪濡存笟鈧顕€宕煎┑瀣暪闂備礁鎼ú銊╁疮閸ф绠繛宸簼閻撶喖鏌ｅΟ鍝勫笭闁煎壊浜弻娑㈠棘鐠恒劎鍔悗娈垮枟閹倿鐛€ｎ喗鏅滈柣锝呰嫰楠炴劙姊虹拠鎻掑毐缂傚秴妫欑粋宥夊醇閺囩喎浜楅梺绋跨箳閳峰牆鈻撴禒瀣厽闁归偊鍨伴惃铏圭磼閻樺樊鐓奸柡宀嬬磿娴狅妇鎷犻幓鎺濈€烽梻渚€鈧偛鑻晶鍓х磼闊厾鐭欓柟顔矫～婵嬵敇瑜庨鏃堟⒑閸涘﹥澶勯柛銊╀憾閸╂盯骞掗幊銊ョ秺閺佹劙宕煎顒佹尵閻ヮ亪顢橀悙闈涚厽闂佸搫鏈惄顖炲极閹邦垳鐤€闁哄洨濮靛▓鍏间繆閵堝洤啸闁稿绋撻幑銏狀潨閳ь剙顕ｉ锕€绠荤紓浣股戝▍銏ゆ⒑鐠恒劌娅愰柟鍑ゆ嫹?)  %%%%%%%%%%%%%%%%%%%%%%%%%%
        vector<vector<complex<double>>> Ezy2(NX, vector<complex<double>>(NY));
        vector<vector<complex<double>>> Eyz2(NX, vector<complex<double>>(NY));

        for(int i=1;i<=NX;i++) 
        {
            for(int j=1;j<=NY;j++)
            {
                Ezy2[i-1][j-1]=(Ez2[i-1][j][Nair+NN]-Ez2[i-1][j-1][Nair+NN]+Ez2[i][j][Nair+NN]-Ez2[i][j-1][Nair+NN])/2.0/fmodel.B_Y[j-1];
                Eyz2[i-1][j-1]=(Ey2[i][j-1][Nair+NN+1]-Ey2[i][j-1][Nair+NN]+Ey2[i-1][j-1][Nair+NN+1]-Ey2[i-1][j-1][Nair+NN])/2.0/fmodel.C_Z[Nair+NN];
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
       // printf("闂傚倸鍊搁崐宄懊归崶顒夋晪鐟滃繘鍩€椤掍胶鈻撻柡鍛█閵嗕礁鈻庨幘鍐插敤濡炪倖鎸鹃崑鐔兼偘閵夆晜鈷戦柛锔诲幖閸斿銇勯妸銉﹀櫧濠㈣娲樼换婵嗩潩椤撶姴骞嶉梻浣虹帛閸旓箓宕滃鑸靛€堕梺顒€绉甸悡娑㈡煃瑜滈崜鐔煎箠閻愬搫唯闁挎繂瀚惄搴ㄦ⒒娴ｅ憡鎯堥柛鐕佸亰瀹曟劙鎮界粙璺唵闂佽法鍠撴慨鐢告偂閺囥垹绠规繛锝庡墮閻忣亪鎮樿箛搴″祮闁哄本娲熷畷鎯邦槻妞ゅ浚鍘鹃埀顒侇問閸犳顭囧▎鎿勭稏婵犻潧顑愰弫鍡涙煃瑜滈崜娑氬垝婵犲浂鏁嬮柍褜鍓熼獮鍐ㄎ旈埀顒勫煝閹捐鍨傛い鏃傛櫕瑜邦垶姊绘担鍛婂暈闁规悂绠栧畷鐗堟償閵婏箑浠奸梺缁樺灱濡嫰鏌嬮崶顒佺厸闁搞儮鏅涢弸鎴炵箾閸涱厽顥犵紒杈ㄦ尰閹峰懘宕烽娑欘潔婵＄偑鍊栭崹鐢稿箠韫囨洜鐭夐柟鐑樻煥婵剟鏌ｉ悪鍛");
        complex<double> c1;
        if(finput1.fail()||k==datasetSize)
        {
           // printf("闂傚倸鍊搁崐鎼佸磹瀹勬噴褰掑炊瑜忛弳锕傛煕椤垵浜濋柛娆忕箻閺屸剝寰勭€ｎ亝顔呭┑鐐叉▕娴滄粓鎮″☉銏＄厱婵犲ň鍋撻柣鎺炵畵瀵煡顢旈崼鐔蜂画濠电姴锕ょ€氼剟鎮橀弶鎴旀斀闁挎稑瀚弳顒侇殽閻愬弶鍠樼€殿喖澧庨幑鍕Ω閵夈垹浜鹃柣鎰劋閳锋垿姊婚崼鐔诲剱鐟滅増甯掔壕璺ㄢ偓瑙勬礀濞层劑鎮虫繝姘厽闁归偊鍠栭崝瀣煟閹惧鈽夋い顓℃硶閹瑰嫰鎮弶鎴濐潬濠电姭鎷冮崶锔藉枤闂佸搫鏈粙鎾诲焵椤掑﹦绉靛ù婊嗘硾鍗遍柛蹇曨儠娴滄粓鏌曢崼婵囧櫣闁逞屽墮閻忔繈锝炶箛鏇犵＜婵☆垵顕ч鎾绘⒑缂佹ê鐏﹂柨姘箾閸繄鐒告慨濠呮閳ь剛娅㈤梽鍕熆濡粯鍙忛柛宀€鍋為悡娆愩亜閺傝濡兼繛璇х畵瀹曟劙鎮滈懞銉у幈濠电娀娼уΛ妤咁敂椤愶附鐓曢柡鍐╂尵閻ｈ鲸銇勯鍕殻濠碘€崇埣瀹曞崬螖閳ь剙顭囬幋锔解拺缂佸顑欓崕鎰版煙閹间胶鐣洪柛鈹惧亾濡炪倖甯掗崰姘焽閹邦厾绠炬繛鏉戭儐濞呭﹥顨ラ悙鑼闁轰焦鎹囬弫鎾绘晸閿燂拷==13\n");
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