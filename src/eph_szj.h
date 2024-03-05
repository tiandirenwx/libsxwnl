#pragma once
#include <string>
#include <vector>
#include "const.h"

/****************************************
以下是天文计算部分,包含有：
物件 SZJ   : 用来计算日月的升起、中天、降落
注意，上述函数或物件是纯天文学的，根据实际需要组合使用可以得到所需要的各种日月坐标，计算精度及计算速度也是可以根据需要有效控制的。
*****************************************/
struct SJ
{
    long double z;
    long double x;
    long double s;
    long double j;
    long double c;
    long double h;
    long double c2;
    long double h2;
    long double c3;
    long double h3;
    long double H0;
    long double H;
    long double H1;
    long double H2;
    long double H3;
    long double H4;
    std::string sm;
};

struct SJ_S
{
    std::string s;
    std::string z;
    std::string j;
    std::string c;
    std::string h;
    std::string ch;
    std::string sj;
    std::string Ms;
    std::string Mz;
    std::string Mj;
};

// 日月的升中天降,不考虑气温和气压的影响
class SunMoonRiseSet
{
public:
    static SunMoonRiseSet &getInstance()
    {
        static SunMoonRiseSet instance;
        return instance;
    }

public:
    const long double Epsilon = 1e-12;
    long double getH(long double h, long double w);
    void Mcoord(long double jd, long double H0, SJ &r);
    void Scoord(long double jd, int xm, SJ &r);
    SJ Mt(long double jd);
    // SJ Qt(long double jd);
    SJ St(long double jd);
    void calcRTS(long double jd, int n, long double Jdl, long double Wdl, long double sq);

    inline long double mod2(long double a, long double b)
    { // 临界余数(a与最近的整倍数b相差的距离)

        long double c = a / b;
        c -= int64(c);
        if (c > 0.5)
        {
            c -= 1;
        }
        return c * b;
    }

    bool areAlmostEqual(long double a, long double b)
    {
        return std::abs(a - b) < Epsilon;
    }

    void setE(long double x);
    long double getE() const;
    void setDt(long double x);
    long double getDt() const;
    void setL(long double x);
    long double getL() const;
    void setFa(long double x);
    long double getFa() const;

private:
    SunMoonRiseSet();
    ~SunMoonRiseSet();
    SunMoonRiseSet(const SunMoonRiseSet &) = delete;
    SunMoonRiseSet &operator=(const SunMoonRiseSet &) = delete;

private:
    long double E;         // 黄赤交角
    long double dt;        // TD-UT
    std::vector<SJ_S> rts; // 多天的升中降
    long double L;         // 站点地理经度,向东测量为正
    long double fa;        // 站点地理纬度
};