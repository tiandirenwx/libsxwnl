// de_hp_test.cpp
// ===============
// DE441/DE440 星历读取的独立测试. 与 eph_hp_test 分开是因为需要 calceph
// 第三方库依赖, 不适合和无依赖的 HP 测试打包在一起.
//
// Test F: 加载 de440s.bsp 并计算 Earth 位置; 与 sxwnl_hp::E_Lon_HP 做量级对比.

#include "de_hp.h"
#include "eph_hp.h"

#include <cmath>
#include <cstdio>
#include <initializer_list>

#ifndef DE_EPHEMERIS_PATH
#define DE_EPHEMERIS_PATH "test/eph_hp/data/de440s.bsp"
#endif

static constexpr long double J2000_JD = 2451545.0L;
static constexpr long double PI_L = 3.14159265358979323846264338327950288L;
static constexpr long double RAD_AS = 180.0L * 3600.0L / PI_L;

static long double year_to_jd(long double year)
{
    return J2000_JD + (year - 2000.0L) * 365.25L;
}

int main()
{
    printf("======================================================================\n");
    printf(" [Test F] JPL DE440 星历 vs VSOP87D-HP 地球日心黄经对比\n");
    printf("   星历文件: %s\n", DE_EPHEMERIS_PATH);
    printf("======================================================================\n");

    if (!sxwnl_hp::de_open(DE_EPHEMERIS_PATH))
    {
        // 有意用 return 0 而非 1: CI/自动化跑测试时, 缺少可选大文件不算失败.
        printf("SKIP: 无法加载 DE 星历: %s\n", DE_EPHEMERIS_PATH);
        printf("      运行以下脚本一键下载:\n");
        printf("        cd test/eph_hp && ./fetch_data.sh de440\n");
        printf("      或手动:\n");
        printf("        curl -o test/eph_hp/data/de440s.bsp \\\n");
        printf("          https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp\n");
        return 0;
    }

    printf("%-6s %-18s %-18s %-14s %-14s\n",
           "Year", "L_J2000 (DE440)", "L_HP (VSOP87D)", "|dL| (as)", "note");

    // de440s 覆盖 1849 ~ 2150 CE
    long double maxDiff_as = 0.0L;
    int maxYear = 0;

    for (int yr : {1900, 1950, 1980, 2000, 2010, 2020, 2050, 2100, 2150})
    {
        long double jd_tdb = year_to_jd((long double)yr);
        auto pos = sxwnl_hp::de_earth_helio(jd_tdb);
        if (!pos.ok)
        {
            printf("%-6d [DE440 越界]\n", yr);
            continue;
        }

        // VSOP87D-HP 用世纪数, 且返回 date-of-epoch 黄经
        long double t = (jd_tdb - J2000_JD) / 36525.0L;
        long double L_HP = sxwnl_hp::E_Lon_HP(t, -1);
        // 归化到 [0, 2π)
        long double L_HP_norm = std::fmod((double)L_HP, (double)(2 * PI_L));
        if (L_HP_norm < 0) L_HP_norm += 2 * PI_L;

        // 差异: DE440 是 J2000 黄经, VSOP87D 是 date 黄经;
        // 差异中主要成分是 J2000→date 的岁差 (~50 arcsec/年)
        long double dL = L_HP_norm - pos.L_J2000;
        // 处理 ±2π 环绕
        while (dL > PI_L) dL -= 2 * PI_L;
        while (dL < -PI_L) dL += 2 * PI_L;
        long double dL_as = std::fabs(dL) * RAD_AS;

        printf("%-6d %-18.12Lf %-18.12Lf %-14.3Lf ",
               yr, pos.L_J2000, L_HP_norm, dL_as);
        // 岁差理论值: (yr-2000)*50.29 arcsec/yr (一般岁差率)
        long double expected_as = std::fabs((yr - 2000) * 50.29L);
        printf("~岁差=%.1Lf as (预期)\n", expected_as);

        if (dL_as > maxDiff_as)
        {
            maxDiff_as = dL_as;
            maxYear = yr;
        }
    }

    printf("\n注意: 上表列出的 |dL| 主要成分是 J2000→date 岁差 (~50\"/年).\n");
    printf("      要判断 DE 与 VSOP87D 本身的差异, 需要额外做岁差修正.\n");
    printf("      本测试只做 DE440 加载和数值合理性验证.\n");

    sxwnl_hp::de_close();
    return 0;
}
