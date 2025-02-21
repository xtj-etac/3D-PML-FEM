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