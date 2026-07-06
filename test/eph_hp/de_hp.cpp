// de_hp.cpp
// ==========
// 使用 IMCCE calceph 读取 JPL DE441/DE440 二进制星历, 提供地球日心位置.
//
// 只在 test/eph_hp/de_hp_test 目标里编译, 不进入 libsxwnl 主库.
//
// calceph 输出 (target=EARTH, center=SUN, unit=AU/AU per day) 是
// **ICRF (J2000 平赤道)** 下的三维直角坐标, 需要:
//   1. 旋转到 J2000 黄道系(即绕 X 轴逆转 ε_J2000 = 23°26'21.406" ≈ 84381.406 arcsec)
//   2. 转球坐标 (L, B, R)
//
// 若将来要与 XL::E_Lon (date-of-epoch 黄经) 直接比较, 还需要:
//   3. 从 J2000 岁差到 date epoch (IAU 2006 岁差角 P03)
// 这一步暂未做, 只输出 J2000 黄经供参考.

#include "de_hp.h"

extern "C" {
#include <calceph.h>
}

#include <cmath>
#include <cstddef>

namespace sxwnl_hp
{

// 单例句柄. 测试是单线程, 无并发问题.
static t_calcephbin *g_eph = nullptr;

// J2000.0 平黄道倾角 (arcsec, IAU 2006 数值 ε₀ = 84381.406")
static constexpr long double EPS0_AS = 84381.406L;
static constexpr long double DAS2R = 4.848136811095359935899141e-6L;  // π/(180·3600)
// 1 AU (km, IAU 2012 决议 B2 定义值)
// SPK 二进制文件默认输出 km, 需要手工换算到 AU 才能和 VSOP87D 的距离比较.
static constexpr long double KM_PER_AU = 149597870.7L;

bool de_open(const char *path)
{
    if (g_eph)
    {
        calceph_close(g_eph);
        g_eph = nullptr;
    }
    g_eph = calceph_open(path);
    if (!g_eph)
    {
        return false;
    }
    // 显式声明 AU, day (SPK 通常已经是 AU/day, 保险起见明示单位)
    calceph_prefetch(g_eph);
    return true;
}

void de_close()
{
    if (g_eph)
    {
        calceph_close(g_eph);
        g_eph = nullptr;
    }
}

DE_HelioPos de_earth_helio(long double jd_tdb)
{
    DE_HelioPos out{0, 0, 0, false};
    if (!g_eph)
    {
        return out;
    }

    // calceph_compute_unit: 使用 KM 长度 + 日 时间单位
    // (SPK 二进制不带 AU 常量, 用 km 最稳妥, 再手动 / KM_PER_AU 转 AU)
    // target=NAIFID_EARTH(399), center=NAIFID_SUN(10), CALCEPH_USE_NAIFID
    double xyz[6] = {0};
    double jd0 = (double)std::floor((double)jd_tdb);
    double jd1 = (double)jd_tdb - jd0;
    int rc = calceph_compute_unit(
        g_eph, jd0, jd1,
        NAIFID_EARTH, NAIFID_SUN,
        CALCEPH_UNIT_KM + CALCEPH_UNIT_DAY + CALCEPH_USE_NAIFID,
        xyz);
    if (!rc)
    {
        return out;
    }
    // xyz[0..2]: 地球位置 (km), xyz[3..5]: 速度 (km/day)
    // 目前是 ICRF (J2000 平赤道)
    long double X = (long double)xyz[0] / KM_PER_AU;
    long double Y = (long double)xyz[1] / KM_PER_AU;
    long double Z = (long double)xyz[2] / KM_PER_AU;

    // ICRF (J2000 平赤道) -> J2000 平黄道:
    //   绕 X 轴顺时针旋转 ε₀ (让赤道抬起到黄道)
    //   x' = x
    //   y' = y·cos ε + z·sin ε
    //   z' = -y·sin ε + z·cos ε
    long double eps = EPS0_AS * DAS2R;
    long double cos_e = std::cos((double)eps);
    long double sin_e = std::sin((double)eps);
    long double xe = X;
    long double ye = Y * cos_e + Z * sin_e;
    long double ze = -Y * sin_e + Z * cos_e;

    long double R = std::sqrt((double)(xe * xe + ye * ye + ze * ze));
    long double L = std::atan2((double)ye, (double)xe);
    if (L < 0) L += 2.0L * 3.14159265358979323846L;
    long double B = std::asin((double)(ze / R));

    out.L_J2000 = L;
    out.B_J2000 = B;
    out.R = R;
    out.ok = true;
    return out;
}

}  // namespace sxwnl_hp
