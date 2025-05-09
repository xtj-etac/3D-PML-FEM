/*
 * @Author: Shutian Shutian@653.com
 * @Date: 2023-11-10 22:46:51
 * @LastEditors: Shutian Shutian@653.com
 * @LastEditTime: 2023-12-13 17:32:44
 * @FilePath: /st/project2/matlabC_v3/HYMAX-APP/codeMAX/rhoAndpha.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#ifndef RHOANDPHA_H
#define RHOANDPHA_H
using namespace std;

class rhoAndpha
{
//private:
    /* data */
    
    vector<vector<complex<double>>> Zxx;
    vector<vector<complex<double>>> Zxy;
    vector<vector<complex<double>>> Zyx;
    vector<vector<complex<double>>> Zyy;

    vector<vector<double>> PhaseXX;
    vector<vector<double>> PhaseXY;
    vector<vector<double>> PhaseYX;
    vector<vector<double>> PhaseYY;

    vector<vector<double>> ResXX;
    vector<vector<double>> ResXY;
    vector<vector<double>> ResYX;
    vector<vector<double>> ResYY;

public:
    rhoAndpha(int NX);
    ~rhoAndpha();

    void rho(v1Assembly v1,v2Assembly v2,Fmodel fmodel,int group_id);
    void printResXY();
};

#endif