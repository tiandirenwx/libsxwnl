#pragma once
#include <string>
#include <array>
#include "eph.h"

struct _SECXY
{
    long double mCJ;
    long double mCW;
    long double mR;
    long double mCJ2;
    long double mCW2;
    long double mR2;

    long double sCJ;
    long double sCW;
    long double sR;
    long double sCJ2;
    long double sCW2;
    long double sR2;

    long double mr;
    long double sr;
    long double x;
    long double y;
    long double t;
};

struct _ZB
{
    Vector3 S;
    Vector3 M;
    long double sr;
    long double mr;
    long double x;
    long double y;
    long double g;
};

struct _GJW
{
    long double J;
    long double W;
    std::string c;
};

typedef struct SolarEclipseMaxData
{
    long double sf;
    long double sf2; // 食分(日出食分)
    long double sf3; // 食分(日没食分)

    long double b1;
    long double dur;
    long double sun_s;
    long double sun_j;
    long double P1;
    long double V1;
    long double P2;
    long double V2;

    bool nasa_r;                   // 为1表示采用NASA的视径比
    std::string sflx;              // 食分类型
    std::array<long double, 5> sT; // 地方日食时间表
    std::string LX;
}SE_MAX_DATA;

typedef struct SolarEclipseNbjData
{
    Vector3 A;
    Vector3 B;  // 本半影锥顶点坐标
    _ZB P;                         // t1时刻的日月坐标,g为恒星时
    _ZB Q;                         // t2时刻的日月坐标
    std::array<long double, 10> V; // 食界表
    std::string Vc;
    std::string Vb; // 食中心类型,本影南北距离
}SE_NBJ_DATA;

class EphRspl
{
public:
     EphRspl();
     ~EphRspl();
    
     void secMax(long double jd, long double L, long double fa, long double high);
     void nbj(long double jd);
     
     long double lineT(_SECXY G, long double v, long double u, long double r, bool n);
     void zbXY(_ZB &p, long double L, long double fa);
     void zb0(long double jd);
     void p2p(long double L, long double fa, _GJW &re, bool fAB, int f);
     void pp0(_GJW &re);
     void secXY(long double jd, long double L, long double fa, long double high, _SECXY &re);

     SE_MAX_DATA &getMaxData();
     SE_NBJ_DATA &getNbjData();
private:
    SE_MAX_DATA *pstSeMax;
    SE_NBJ_DATA *pstNbj;
};