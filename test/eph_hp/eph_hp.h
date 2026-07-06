// eph_hp.h
// ========
// 高精度(High-Precision)天体位置计算接口 —— "备件"版本
//
// 目的:
//   把 VSOP87 全项(2564 项,而非现有 XL0_0 截断版的 ~890 项)封装成一组
//   与 src/eph.h 中 XL::E_Lon() 接口对齐的函数,方便将来在需要更高
//   精度时直接切换,而不需改动生产代码.
//
// 使用范围:
//   目前仅提供 地球日心黄经/黄纬/距离 的高精度版本 (对应现有 XL0_0)
//   其他行星、月球、章动、岁差等仍复用 src/ 里的现有实现.
//   将来若接入 IAU 2000A 全项章动,可在本文件同名 namespace 内追加.
//
// 依赖:
//   - test/eph_hp/vsop87d_earth_hp.h    (由 tools/vsop87d_to_header.py 生成)
//     数据源: IMCCE VSOP87D.ear (Date-of-epoch, 与现有 XL0_0 同坐标系)
//   - test/eph_hp/iau2000a_nut_lon_hp.h (由 tools/iau2000a_to_header.py 生成)
//     数据源: IERS Conventions 2010 Chapter 5 Table 5.3a (章动黄经全项)
//   - src/eph.h 的 gxc_sunLon() 保留复用 (光行差, 与章动/VSOP 独立)
//
// 编译方式:
//   本文件不会被打包进主库 libsxwnl,仅在 pingqi_dingqi 之外的独立
//   test/eph_hp/eph_hp_test 目标里编译,防止影响生产代码.

#pragma once

namespace sxwnl_hp
{

// 计算地球日心黄道坐标某一分量(J2000 分点黄道, VSOP87B 完整项).
//   zn: 0 = 黄经 L, 1 = 黄纬 B, 2 = 距离 R
//   t : 儒略世纪数(自 J2000)
//   n : 项数选择规则. n<0 使用全部项;n>=0 走与 XL0_calc 相同的自适应
//       截断规则,便于和现有实现做同口径对比
long double XL0_calc_HP(int zn, long double t, int n);

// 地球日心黄经(高精度),对应 XL::E_Lon
long double E_Lon_HP(long double t, int n);

// IAU 2000A 全项章动黄经(单位: 弧度)
// 与 src/eph.cpp 中中精度 nutationLon2() 对应, 项数 1320+38 vs ~22
long double nutationLon_HP(long double t);

// 太阳视黄经(高精度),对应 XL::S_aLon
// 内部使用: E_Lon_HP + nutationLon_HP + gxc_sunLon + π
// 即 VSOP 数据源和章动模型都升级, 光行差仍复用中精度实现
long double S_aLon_HP(long double t, int n);

// 由已知太阳视黄经反求时间(高精度),对应 XL::S_aLon_t
// 用两步 Newton 迭代,与 src/eph.cpp 的 S_aLon_t 结构相同
long double S_aLon_t_HP(long double W);

}  // namespace sxwnl_hp
