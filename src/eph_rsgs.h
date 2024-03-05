#pragma once
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include "eph.h"

struct _VXY
{
    long double vx;
    long double vy;
    long double Vx;
    long double Vy;
    long double V;
};

struct _RSM
{
    long double r1;
    long double r2;
    long double ar2;
    long double sf;
};

struct _FEATURE
{
    long double jdSuo;
    long double dT;
    long double ds;
    long double vx;
    long double vy;
    long double ax;
    long double ay;
    long double v;
    long double k;

    long double t0;
    long double jd;
    long double xc;
    long double yc;
    long double zc;
    long double D;
    long double d;
    Vector3 I;
    Vector3 gk1;
    Vector3 gk2;
    Vector3 gk3;
    Vector3 gk4;
    Vector3 gk5;
    std::string lx;

    long double zxJ;
    long double zxW;
    long double dw;
    long double sf;
    long double tt;
    Vector3 Sdp;

    std::vector<long double> p1;
    std::vector<long double> p2;
    std::vector<long double> p3;
    std::vector<long double> p4;
    std::vector<long double> q1;
    std::vector<long double> q2;
    std::vector<long double> q3;
    std::vector<long double> q4;
    std::vector<long double> L0;
    std::vector<long double> L1;
    std::vector<long double> L2;
    std::vector<long double> L3;
    std::vector<long double> L4;
    std::vector<long double> L5;
    std::vector<long double> L6;
};

struct _JIEX2
{
    std::vector<long double> p1;
    std::vector<long double> p2;
    std::vector<long double> p3;
};

struct _FLAG
{
    int f;
    int f2;
};

class RSGS_PARA
{
public:   
    void setZjd(long double jd)
    {
        Zjd = jd;
    }
    long double getZjd() const
    {
        return Zjd;
    }
    void setZdt(long double dt)
    {
        Zdt = dt;
    }

    long double getZdt() const
    {
        return Zdt;
    }
    void setdT(long double dt)
    {
        dT = dt;
    }
    long double getdT() const
    {
        return dT;
    }

    void setTanf1(long double f)
    {
        tanf1 = f;
    }
    long double getTanf1() const
    {
        return tanf1;
    }

    void setTanf2(long double f)
    {
        tanf2 = f;
    }
    long double getTanf2() const
    {
        return tanf2;
    }

    void setSrad(long double d)
    {
        srad = d;
    }
    long double getSrad() const
    {
        return srad;
    }

    void setBba(long double a)
    {
        bba = a;
    }
    long double getBba() const
    {
        return bba;
    }

    void setBhc(long double b)
    {
        bhc = b;
    }
    long double getBhc() const
    {
        return bhc;
    }

    void setDyj(long double d)
    {
        dyj = d;
    }
    long double getDyj() const
    {
        return dyj;
    }

private:
    long double Zjd;
    long double Zdt;
    long double dT;
    long double tanf1;
    long double tanf2;
    long double srad;
    long double bba;
    long double bhc;
    long double dyj;
};

class EphRsgs
{
public:
    static EphRsgs &getInstance()
    {
        static EphRsgs instance;
        return instance;
    }

public:
    void init(long double jd, int n);
    _FEATURE feature(long double jd);
    // static _FEATURE __rsGS::jieX(long double jd);
    // static _JIEX2 __rsGS::jieX2(long double jd);
    std::string jieX3(long double jd);
    RSGS_PARA &getRsgsParaData();
    inline Vector3 sun(long double jd) { return chazhi(jd, 0); } // 传回值可能超过360度
    inline Vector3 moon(long double jd) { return chazhi(jd, 1); }
    inline Vector3 bse(long double jd) { return chazhi(jd, 2); }

private:
    EphRsgs();
    ~EphRsgs();
    EphRsgs(const EphRsgs &) = delete;
    EphRsgs &operator=(const EphRsgs &) = delete;

private:
    std::vector<long double> Zs; // 日月赤道坐标插值表
    RSGS_PARA *pstPara;

    inline bool approximatelyEqual(long double a, long double b,long double epsilon = std::numeric_limits<long double>::epsilon() )
    {
        return std::abs(a - b) < epsilon;
    }

    inline bool notEqual(long double a, long double b, long double epsilon = std::numeric_limits<long double>::epsilon()) 
    {
        return std::abs(a - b) > epsilon;
    }

    Vector3 chazhi(long double jd, int xt);
    Vector3 cd2bse(Vector3 z, Vector3 I);
    Vector3 bse2cd(Vector3 z, Vector3 I);
    Vector3 bse2db(Vector3 z, Vector3 I, bool f);
    Vector3 bseXY2db(long double x, long double y, Vector3 I, bool f);
    Vector3 bseM(long double jd);
    _VXY Vxy(long double x, long double y, long double s, long double vx, long double vy);
    _RSM rSM(long double mR);
    Vector3 qrd(long double jd, long double dx, long double dy, bool fs);
    void push(Vector3 z, std::vector<long double> &p);
    Vector4 nanbei(Vector3 M, long double vx0, long double vy0, long double h, long double r, Vector3 I);
    bool mDian(Vector3 M, long double vx0, long double vy0, bool AB, long double r, Vector3 I, std::vector<long double> &A);
    // static void __rsGS::elmCpy(std::vector<long double> &a,int n,std::vector<long double> b,int m);
    // static void __rsGS::mQie(Vector3 M,long double vx0,long double vy0,long double h, long double r,Vector3 I, std::vector<long double> &A,_FLAG &FLAG);
};