#pragma once
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include "const.h"
#include "JD.h"
#include "eph.h"
#include "sx_lang_zh.h"
#include "geo.h"
void pCalc(int xt, Time &dat, int n = 10, int dt = 1, bool Cd_ut = true)
{                                      // 行星星历计算
    double jd = JD::toJD(dat) - J2000; // 取屏幕时间
    if (Cd_ut)
    {
        jd += -8.0 / 24 + dt_T(jd); // 转为力学时
    }
    GeoPostion &gep = GeoPostion::getInstance();
	auto jw = gep.getDefaultGeoPos();
    double L = jw.J / 180 * M_PI; // 地标
    double fa = jw.W / 180 * M_PI;
    if (n > 1000)
    {
        std::cout << "个数太多了" << std::endl;
        return;
    }
    std::string s = "";
    int i;
    // 求星历
    s += "**************************" + std::string(EphemerisList[xt]) + "星历**************************\n";

    for (i = 0; i < n; i++, jd += dt)
    {
        double jd2 = jd + J2000;
        std::ostringstream oss;
        oss<< std::fixed<< std::setprecision(7)<< jd2;
        std::string result = oss.str();
        s += JD::JD2str(jd2) + "TD, JED = " + result + " " + "\n";
        s += XL::xingX(xt, jd, L, fa) + "\n";
    }
    std::cout << s << std::endl;
}
