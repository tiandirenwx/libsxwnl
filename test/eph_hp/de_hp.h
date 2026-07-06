// de_hp.h
// ========
// JPL DE441 / DE440 星历包装 (使用 IMCCE calceph 库读取二进制文件).
//
// 位置: 只在 test/eph_hp/de_hp_test 目标里编译, 不进入 libsxwnl 主库.
//       主库使用者(Android/iOS 集成方)不需要 calceph 依赖.
//
// 覆盖时段:
//   de440s.bsp : 1849-2150 CE (32 MB, 快速获得可用星历)
//   de440.bsp  : 1550-2650 CE (~110 MB)
//   de441.bsp  : -13200 ~ 17191 CE (~3.2 GB, 完整长时段, 需 lfs 或按需下载)
//   任何 SPK/INPOP 格式二进制均可.

#pragma once

namespace sxwnl_hp
{

struct DE_HelioPos
{
    long double L_J2000;  // 地球日心 J2000 黄道黄经 (弧度)
    long double B_J2000;  // 地球日心 J2000 黄道黄纬 (弧度)
    long double R;        // 距离 (AU)
    bool ok;              // false 表示未初始化 / JD 越界
};

// 加载星历文件. 支持 SPK (.bsp), INPOP (.dat) 等 calceph 认识的格式.
// 返回 true 表示加载成功. 生命周期由 sxwnl_hp 内部维护(线程不安全,
// 单例用法即可, 目前是给测试用).
bool de_open(const char *ephemeris_path);
void de_close();

// 计算地球相对太阳的位置 (heliocentric).
// 输入: TDB (Barycentric Dynamical Time) 儒略日
// 输出: J2000 惯性系(ICRF) 下的 (L, B, R). 若未打开或越界返回 ok=false
DE_HelioPos de_earth_helio(long double jd_tdb);

}  // namespace sxwnl_hp
