#pragma once
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <chrono>
#include "eph.h"
#include "JD.h"
#include "sx_lang_zh.h"

void dingQi_v()
{ // 定气计算速度测试
    int i;
    auto d1 = std::chrono::system_clock::now();
    for (i = 0; i < 10000; i++)
    {
        XL::S_aLon_t(0);
    }
    auto d2 = std::chrono::system_clock::now();
    for (i = 0; i < 10000; i++)
    {
        XL::S_aLon_t2(0);
    }

    auto d3 = std::chrono::system_clock::now();
    std::cout << "定气测试\n\033[31;1m高精度:" << std::chrono::duration<double>(d2 - d1).count() * 1000 << "毫秒/10千个\n"
              << "低精度:" << std::chrono::duration<double>(d3 - d2).count() * 1000 << "毫秒/10千个\n\033[0m";
}

void dingSuo_v()
{ // 定朔计算速度测试
    int i;
    auto d1 = std::chrono::system_clock::now();
    for (i = 0; i < 10000; i++)
    {
        XL::MS_aLon_t(0);
    }

    auto d2 = std::chrono::system_clock::now();
    for (i = 0; i < 10000; i++)
    {
        XL::MS_aLon_t2(0);
    }

    auto d3 = std::chrono::system_clock::now();
    std::cout << "\033[31;1m高精度:" << std::chrono::duration<double>(d2 - d1).count() * 1000 << "毫秒/10千个\n"
              << "低精度:" << std::chrono::duration<double>(d3 - d2).count() * 1000 << "毫秒/10千个\n\033[0m";
}

void suoCalc(int y, int n = 24, int jiao = 0)
{ // 定朔测试函数
    y -= 2000;
    if (jiao == -1)
    {
        std::cout << "请输入角度(0朔,90上弦,180望,270下弦,或其它):" << std::endl;
        std::cin >> jiao;
    }
    int i;
    double r, T;
    std::string s = "月-日黄经差" + std::to_string(jiao) + "\n", s2 = "";
    int n0 = int2(y * (365.2422 / 29.53058886)); // 截止当年首经历朔望的个数
    for (i = 0; i < n; i++)
    {
        T = XL::MS_aLon_t((n0 + i + jiao / 360.0) * 2 * M_PI);                                                    // 精确时间计算,入口参数是当年各朔望黄经
        r = XL1_calc(2, T, -1); 
        std::ostringstream oss;
        oss<< std::fixed<< std::setprecision(3)<< r;
        std::string result = oss.str();                                                                                  // 计算月亮
        s2 += JD::JD2str(T * 36525 + J2000 + 8.0 / 24 - dt_T(T * 36525)) + " " + result + "千米\n"; // 日期转为字串
        if (i % 50 == 0)
        {
            s += s2, s2 = "";
        }
    }

    std::cout << s + s2 << std::endl;
}

void qiCalc(int y, int n = 24)
{ // 定气测试函数
    y -= 2000;
    int i;
    double T;
    std::string s = "", s2 = "";

    for (i = 0; i < n; i++)
    {
        T = XL::S_aLon_t((y + i * 15 / 360.0 + 1) * 2 * M_PI);                                 // 精确节气时间计算
        s2 += JD::JD2str(T * 36525 + J2000 + 8 / 24.0 - dt_T(T * 36525)) + Jqmc[(i + 6) % 24]; // 日期转为字串
        if (i % 2 == 1)
        {
            s2 += " 视黄经" + std::to_string(i * 15) + "\n";
        }
        else
        {
            s2 += " ";
        }
        if (i % 50 == 0)
        {
            s += s2, s2 = "";
        }
    }
    std::cout << s + s2 << std::endl;
}

void houCalc(int y, int n = 24)
{ // 定候测试函数
    y -= 2000;
    int i;
    double T;
    std::string s = "初候　　　　　　　　　　　　二候　　　　　　　　　三候", s2 = "";
    for (i = 0; i < n * 3; i++)
    {
        T = XL::S_aLon_t((y + i * 5 / 360.0 + 1) * 2 * M_PI); // 精确节气时间计算
        if (i % 3 == 0)
        {
            s2 = s2 + "\n" + Jqmc[(i / 3 + 6) % 24];
        }
        else
        {
            s2 += " ";
        }
        s2 += JD::JD2str(T * 36525 + J2000 + 8.0 / 24.0 - dt_T(T * 36525)); // 日期转为字串
        if (i % 50 == 0)
        {
            s += s2, s2 = "";
        }
    }
    std::cout << s + s2 << std::endl;
}

void dingQi_cmp(int y = 2000, int N = 10)
{ // 定气误差测试
    y -= 2000;
    int i;
    double W, T, maxT = 0;
    for (i = 0; i < N; i++)
    {
        W = (y + i / 24) * 2 * M_PI;
        T = XL::S_aLon_t2(W) - XL::S_aLon_t(W); // 节气粗算与精算的差异
        T = int2(abs(T * 36525 * 86400));
        if (T > maxT)
        {
            maxT = T;
        }
    }
    //输出一个设置终端文本颜色为红色并且粗体的空字符串
    std::cout << "\033[31;1m" << std::to_string(2000 + y) << "年之后" << std::to_string(N) << "个朔日粗算与精算的最大差异:" << std::to_string(maxT) << "秒。\033[0m" << std::endl;
}

void dingSuo_cmp(int y = 2000, int N = 10)
{ // 定朔测试函数
    y -= 2000;
    int i;
    double T, maxT = 0, W;
    int n = int2(y * (365.2422 / 29.53058886)); // 截止当年首经历朔望的个数
    for (i = 0; i < N; i++)
    {
        W = (n + i / 24.0) * 2 * M_PI;
        T = XL::MS_aLon_t2(W) - XL::MS_aLon_t(W); // 合塑粗算与精算的差异
        T = int2(abs(T * 36525 * 86400));
        if (T > maxT)
        {
            maxT = T;
        }
    }
    std::cout << "\033[31;1m" << std::to_string(2000 + y) << "年之后" << std::to_string(N) << "个朔日粗算与精算的最大差异:" << std::to_string(maxT) << "秒。\033[0m" << std::endl;
}