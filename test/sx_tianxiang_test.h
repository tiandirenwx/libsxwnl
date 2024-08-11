#pragma once
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include "const.h"
#include "JD.h"
#include "eph.h"
#include "sx_lang_zh.h"

std::string txFormatT(long double t)
{ // 天象时间格式化输出
    long double t1 = t * 36525 + J2000;
    long double t2 = t1 - dt_T(t1 - J2000) + 8.0 / 24;
    return JD::JD2str(t1) + " TD " + JD::JD2str(t2).substr(9, 11) + " UT ";
}

void tianXiang(int xm, int xm2, Time &dat, int n = 10)
{
    static const char *xxName[] = {"地球", "水星", "金星", "火星", "木星", "土星", "天王星", "海王星", "冥王星"};
    double jd = JD::toJD(dat) - J2000; // 取屏幕时间
    std::string s = "";
    int i;
    long double re0;
    Vector2 re;
    Vector4 re2;
    jd /= 36525.0;
    if (xm == 1 || xm == 2)
    { // 求月亮近远点
        if (xm == 1)
        {
            s += "---------------------月亮近点-----------------------\n";
        }

        if (xm == 2)
        {
            s += "---------------------月亮远点-----------------------\n";
        }
        
        
        for (i = 0; i < n; i++, jd = re[0] + 27.555 / 36525.0)
        {
            if (xm == 1)
            {
                re = XL::moonMinR(jd, true); // 求近点
            }
            if (xm == 2)
            {
                re = XL::moonMinR(jd, false); // 求远点
            }
            s += txFormatT(re[0]) + toFixed(re[1], 2) + "千米\n";
        }
    }
    if (xm == 3 || xm == 4)
    { // 求月亮升降交点
        if (xm == 3)
        {
            s += "---------------------月亮升交-----------------------\n";
        }

        if (xm == 4)
        {
            s += "---------------------月亮降交-----------------------\n";
        }
        
        for (i = 0; i < n; i++, jd = re[0] + 27.555 / 36525.0)
        {
            if (xm == 3)
            {
                re = XL::moonNode(jd, true); // 求升
            }
            if (xm == 4)
            {
                re = XL::moonNode(jd, false); // 求降
            }
            s += txFormatT(re[0]) + rad2str(rad2mrad(re[1]), false) + "\n";
        }
    }
    if (xm == 5 || xm == 6)
    { // 求地球近远点
        if (xm == 5)
        {
            s += "---------------------地球近日-----------------------\n";
        }

        if (xm == 6)
        {
            s += "---------------------地球远日-----------------------\n";
        }
        
        for (i = 0; i < n; i++, jd = re[0] + 365.259636 / 36525.0)
        {
            if (xm == 5)
            {
                re = XL::earthMinR(jd, true); // 求近点
            }
            if (xm == 6)
            {
                re = XL::earthMinR(jd, false); // 求远点
            }
            s += txFormatT(re[0]) + toFixed(re[1], 8) + " AU\n";
        }
    }
    if (xm == 7 || xm == 8)
    { // 大距计算
        if (xm == 7)
        {
            s += "---------------------水东大距-----------------------\n";
        }

        if (xm == 8)
        {
            s += "---------------------水西大距-----------------------\n";
        }

        for (i = 0; i < n; i++, jd = re[0] + 115.8774777586 / 36525.0)
        {
            if (xm == 7)
            {
                re = XL::daJu(1, jd, true); // 求水星东大距
            }
            if (xm == 8)
            {
                re = XL::daJu(1, jd, false); // 求水星东西距
            }
            s += txFormatT(re[0]) + toFixed((re[1] / M_PI * 180), 5) + "度\n";
        }
    }
    if (xm == 9 || xm == 10)
    { // 大距计算
        if (xm == 9)
        {
            s += "---------------------金东大距-----------------------\n";
        }

        if (xm == 10)
        {
            s += "---------------------金西大距-----------------------\n";
        }
        for (i = 0; i < n; i++, jd = re[0] + 583.9213708245 / 36525.0)
        {
            if (xm == 9)
            {
                re = XL::daJu(2, jd, true); // 求金星东大距
            }
            if (xm == 10)
            {
                re = XL::daJu(2, jd, false); // 求金星东西距
            }
            s += txFormatT(re[0]) + toFixed((re[1] / M_PI * 180), 5) + "度\n";
        }
    }
    if (xm == 11)
    { // 合月计算
        s += std::string(xxName[xm2]) + "赤经合月\n";
        s += "合月时间(TD UT) 星月赤纬差(小于1度可能月掩星,由视差决定)\n";
        for (i = 0; i < n; i++, jd = re[0] + 28.0 / 36525.0)
        {
            re = XL::xingHY(xm2, jd);
            s += txFormatT(re[0]) + toFixed((-re[1] / M_PI * 180), 5) + "度\n";
        }
    }
    if (xm == 12 || xm == 13)
    {
        if (xm == 12)
        {
            s = std::string(xxName[xm2]) + "合日(地内行星上合)\n";
        }
        if (xm == 13)
        {
            s = std::string(xxName[xm2]) + "冲日(地内行星下合)\n";
        }
        s += "黄经合/冲日时间(TD UT) 星日赤纬差\n";
        for (i = 0; i < n; i++, jd = re[0] + cs_xxHH[xm2 - 1] / 36525.0)
        {
            if (xm == 12)
            {
                re = XL::xingHR(xm2, jd, false);
            }
            if (xm == 13)
            {
                re = XL::xingHR(xm2, jd, true);
            }
            s += txFormatT(re[0]) + toFixed((-re[1] / M_PI * 180), 5) + "度\n";
        }
    }

    if (xm == 14 || xm == 15)
    { // 顺留
        if (xm == 14)
        {
            s = std::string(xxName[xm2]) + "顺留\n";
        }
        if (xm == 15)
        {
            s = std::string(xxName[xm2]) + "逆留\n";
        }
        s += "留时间(TD UT)\n";
        for (i = 0; i < n; i++, jd = re0 + cs_xxHH[xm2 - 1] / 36525.0)
        {
            if (xm == 14)
            {
                re0 = XL::xingLiu(xm2, jd, true);
            }
            if (xm == 15)
            {
                re0 = XL::xingLiu(xm2, jd, false);
            }
            s += txFormatT(re0) + "\n";
        }
    }
    std::cout << s << std::endl;
}
