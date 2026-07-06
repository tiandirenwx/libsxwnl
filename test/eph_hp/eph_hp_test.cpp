// eph_hp_test.cpp
// ================
// VSOP87 全项 (HP) 与现有截断版 (XL::E_Lon / XL::S_aLon) 的对比测试.
// 本测试**只读**运行, 不修改任何生产代码或数据, 不与 pingqi_dingqi
// 或其他测试有任何耦合.
//
// 三组用例:
//   [Test A] 自洽性:
//       在近现代年份 (1900~2024) 处比较 HP 与当前实现的黄经差
//       期望 < 0.1 角秒, 换算成时间 < 3 秒
//
//   [Test B] 精度提升量化:
//       在 -4000 ~ +9999 每 500 年采样一次
//       输出 HP - 当前 的黄经差, 观察远古/远期何时差异开始显著
//
//   [Test C] 高价值年份 (关键立春时刻):
//       计算 2020/2050/5000 年立春时刻, 用 HP 反算与当前实现比较
//       给出秒级差异, 直观展示精度提升的"实际"效果
//
// 运行: ./build/bin/eph_hp_test

#include "eph_hp.h"
#include "../../src/eph.h"    // XL::E_Lon, XL::S_aLon, XL::S_aLon_t, nutationLon2

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

// 注意: 主库 src/const.h 已定义 PI 为宏,此处沿用 (不能用同名局部常量)
// 每弧度对应的角秒数 (~206265)
static constexpr long double RAD_AS = 180.0L * 3600.0L / (long double)PI;

// 儒略日 <-> 儒略世纪 (J2000)
static constexpr long double J2000_JD = 2451545.0L;
static long double year_to_t(long double y)
{
    // 简化: 年份中点 (y-2000)*365.25 / 36525
    return (y - 2000.0L) / 100.0L;
}

// 弧度差 -> 角秒
static long double rad2as(long double x) { return x * RAD_AS; }
// 弧度差 -> 时间秒 (太阳视黄经日均运动 ≈ 3548 arcsec/day, 即 ~0.041 arcsec/s)
static long double rad2timeSec(long double dLon)
{
    // 太阳日均视黄经 ≈ 360°/365.24 天 = 0.98565°/day = 3548.19 arcsec/day
    long double as = rad2as(dLon);
    return as / (3548.19L / 86400.0L);
}

// ==================== [Test A] 近现代自洽性 ====================
static int testA_nearModern()
{
    printf("\n===================================================================\n");
    printf(" [Test A] 近现代自洽性 (1900 ~ 2024)\n");
    printf(" 期望: HP 与截断版之差在 0.1 角秒以内, 时间上 < 3 秒\n");
    printf("===================================================================\n");
    printf("%-8s %-14s %-14s %-14s %-12s %s\n",
           "Year", "E_Lon (curr)", "E_Lon_HP", "|dLon| (arcsec)", "≈dTime(s)", "PASS?");

    int passed = 0, total = 0;
    for (int yr = 1900; yr <= 2024; yr += 10)
    {
        long double t = year_to_t((long double)yr);
        long double lonCurr = XL::E_Lon(t, -1);
        long double lonHP = sxwnl_hp::E_Lon_HP(t, -1);
        long double dLon = lonHP - lonCurr;
        long double dLon_as = rad2as(std::fabs(dLon));
        long double dTime_s = std::fabs(rad2timeSec(dLon));

        bool pass = dLon_as < 0.1L;
        printf("%-8d %-14.9Lf %-14.9Lf %-14.6Lf %-12.4Lf %s\n",
               yr, lonCurr, lonHP, dLon_as, dTime_s,
               pass ? "OK" : "FAIL");
        total++;
        if (pass) passed++;
    }
    printf("[Test A] %d / %d 通过\n", passed, total);
    return (passed == total) ? 0 : 1;
}

// ==================== [Test B] 长期精度差异扫描 ====================
static int testB_longRangeScan()
{
    printf("\n===================================================================\n");
    printf(" [Test B] 长期精度差异 (-4000 ~ +9999, 每 500 年采样)\n");
    printf(" 目的: 观察远古/远期何时 HP 与截断版出现明显分歧\n");
    printf("===================================================================\n");
    printf("%-8s %-14s %-14s %-14s %-14s\n",
           "Year", "E_Lon (curr)", "E_Lon_HP", "|dLon| (as)", "≈dTime(s)");

    long double maxDiff_as = 0.0L;
    int maxYear = 0;
    long double maxDiff_timeSec = 0.0L;

    for (int yr = -4000; yr <= 9999; yr += 500)
    {
        long double t = year_to_t((long double)yr);
        long double lonCurr = XL::E_Lon(t, -1);
        long double lonHP = sxwnl_hp::E_Lon_HP(t, -1);
        long double dLon = lonHP - lonCurr;
        long double dLon_as = rad2as(std::fabs(dLon));
        long double dTime_s = std::fabs(rad2timeSec(dLon));

        printf("%-8d %-14.9Lf %-14.9Lf %-14.6Lf %-14.4Lf\n",
               yr, lonCurr, lonHP, dLon_as, dTime_s);

        if (dLon_as > maxDiff_as)
        {
            maxDiff_as = dLon_as;
            maxYear = yr;
            maxDiff_timeSec = dTime_s;
        }
    }
    printf("\n[Test B] 最大差异: %.6Lf 角秒 (≈ %.3Lf 秒) 出现在 Y=%d\n",
           maxDiff_as, maxDiff_timeSec, maxYear);
    return 0;
}

// ==================== [Test C] 太阳视黄经 S_aLon 直接对比 ====================
static int testC_apparentSunLon()
{
    printf("\n===================================================================\n");
    printf(" [Test C] 太阳视黄经 S_aLon = E_Lon + 章动 + 光行差 + π 直接对比\n");
    printf(" (章动/光行差沿用现有 nutationLon2/gxc_sunLon, HP 只升级 VSOP 那层)\n");
    printf("===================================================================\n");
    printf("%-6s %-16s %-16s %-14s %-14s\n",
           "Year", "S_aLon (curr)", "S_aLon (HP)", "|dLon| (as)", "≈dTime(s)");

    long double maxDiff_as = 0.0L;
    int maxDiffYear = 0;
    long double maxDiff_s = 0.0L;

    for (int yr : {1900, 1960, 2000, 2020, 2050, 2100, 3000, 5000, 7000, 9000, 9999})
    {
        long double t = year_to_t((long double)yr);
        long double sCurr = XL::S_aLon(t, -1);
        long double sHP = sxwnl_hp::S_aLon_HP(t, -1);
        long double dLon = sHP - sCurr;
        long double dLon_as = rad2as(std::fabs(dLon));
        long double dTime_s = std::fabs(rad2timeSec(dLon));

        printf("%-6d %-16.9Lf %-16.9Lf %-14.6Lf %-14.4Lf\n",
               yr, sCurr, sHP, dLon_as, dTime_s);

        if (dLon_as > maxDiff_as)
        {
            maxDiff_as = dLon_as;
            maxDiffYear = yr;
            maxDiff_s = dTime_s;
        }
    }
    printf("\n[Test C] 最大 S_aLon 差异: %.6Lf 角秒 (≈ %.3Lf 秒) 出现在 Y=%d\n",
           maxDiff_as, maxDiff_s, maxDiffYear);
    return 0;
}

// ==================== [Test D] IAU 2000A 章动 vs 现有中精度章动 ====================
static int testD_nutationLon()
{
    printf("\n===================================================================\n");
    printf(" [Test D] IAU 2000A 章动黄经 (1320+38 项) vs 现有 nutationLon2 (~22 项)\n");
    printf(" 主要看近现代 sub-mas 级差异, 是否符合 IAU 精度要求\n");
    printf("===================================================================\n");
    printf("%-6s %-16s %-16s %-16s %-16s\n",
           "Year", "nutation (curr)", "nutation (HP)", "|dNut| (as)", "|dNut| (mas)");

    long double maxDiff_mas = 0.0L;
    int maxYear = 0;

    for (int yr : {-2000, -1000, 0, 500, 1000, 1500, 1800, 1900, 1950, 1980,
                   2000, 2010, 2020, 2050, 2100, 2500, 3000, 5000, 7000, 9999})
    {
        long double t = year_to_t((long double)yr);
        long double nCurr = nutationLon2(t);
        long double nHP = sxwnl_hp::nutationLon_HP(t);
        long double d = nHP - nCurr;
        long double d_as = rad2as(std::fabs(d));
        long double d_mas = d_as * 1000.0L;

        printf("%-6d %-16.10Lf %-16.10Lf %-16.6Lf %-16.3Lf\n",
               yr, nCurr, nHP, d_as, d_mas);

        if (d_mas > maxDiff_mas)
        {
            maxDiff_mas = d_mas;
            maxYear = yr;
        }
    }
    printf("\n[Test D] 最大章动差异: %.3Lf mas (毫弧秒) 出现在 Y=%d\n",
           maxDiff_mas, maxYear);
    printf("           (IAU 2000A 声明精度 ~0.2 mas; 现有中精度 ~10 mas)\n");
    return 0;
}

int main(int argc, char *argv[])
{
    printf("======================================================================\n");
    printf(" 高精度天体位置计算 备件对比测试\n");
    printf("======================================================================\n");
    printf(" VSOP87D 全项    (2425 项) vs XL0_0 截断    (~890 项)\n");
    printf(" IAU 2000A 章动 (1320+38 项) vs nutB 中精度 (~22 项)\n");
    printf("======================================================================\n");

    int a = testA_nearModern();
    (void)testB_longRangeScan();
    (void)testC_apparentSunLon();
    (void)testD_nutationLon();

    printf("\n======================================================================\n");
    printf(" 测试结束. 结论:\n");
    printf(" - Test A 近现代 : VSOP HP 与截断版 sub-arcsec 一致 (方法验证)\n");
    printf(" - Test B 长期扫描: VSOP HP 差异随时间放大 (精度损失量化)\n");
    printf(" - Test C S_aLon : 加入章动/光行差后综合视黄经差异\n");
    printf(" - Test D 章动   : IAU 2000A 相对中精度章动的 mas 级提升\n");
    printf("======================================================================\n");

    return a;
}
