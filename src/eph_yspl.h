#pragma once
#include <array>
#include <string>

struct RE0
{
    long double e_mRad;
    long double eShadow;
    long double eShadow2;
    long double x;
    long double y;
    long double mr;
    long double er;
    long double Er;
    long double t;
};

// 月食快速计算器
class EphYspl
{
public:
    EphYspl();
    ~EphYspl();

public:
    void lecXY(long double jd, RE0 &re);
    long double lineT(RE0 G, long double v, long double u, long double r, bool n);
    void lecMax(long double jd);
    std::string getLX() const;
    long double getSF() const;
    std::array<long double,7> getLT() const;

private:
    std::array<long double, 7> lT;
    std::string LX;
    long double sf;
};