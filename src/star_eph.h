#pragma once

#include <array>
#include "eph.h"

// ═══════════════════════════════════════════════════════════════
//  star_eph —— 恒星星历计算核心
//
//  对应上游 sxwnl/ephB.js, 移植以下函数:
//    evSSB(t)    : 地球 SSB 速度(直角坐标, AU/儒略世纪)
//    epSSB(t)    : 地球 SSB 位置(直角坐标, AU)
//    ylpz(z, a)  : 引力偏转修正
//    scGxc(z,v,f): 严格周年视差(f=1) 或 光行差(f=0) 改正
//    sun2000(t)  : 太阳 J2000 球面坐标
//
//  精度: 全部使用 long double, 与 libsxwnl 其它模块一致。
//  坐标说明:
//    - SSB(Solar System Barycenter, 太阳系质心) 直角坐标系, J2000 赤道
//    - 球面坐标 z[0]=赤经/黄经, z[1]=赤纬/黄纬, z[2]=向径
// ═══════════════════════════════════════════════════════════════

namespace star_eph
{

// 地球的 SSB 速度(直角坐标, AU/儒略世纪).
// t: 儒略世纪TD, 自 J2000 起算。平均误差 4×10⁻⁸ AU/日。
Vector3 evSSB(long double t);

// 地心相对 SSB 的直角坐标位置(AU). 4 位有效数字精度.
Vector3 epSSB(long double t);

// 引力偏转修正
//   z  : 天体的赤道球面坐标 (赤经, 赤纬, 向径)
//   sunEq: 太阳的赤道球面坐标 (赤经, 赤纬, 向径)
//   返回修正后的球面坐标
Vector3 ylpz(const Vector3 &z, const Vector3 &sunEq);

// 严格的恒星视差或光行差改正
//   z: 天体赤道球面坐标 (J2000)
//   v: 地球赤道直角坐标 (J2000)
//   parallax = false → 光行差改正(v 须为 SSB 速度)
//   parallax = true  → 周年视差改正(v 须为 SSB 位置)
Vector3 scGxc(const Vector3 &z, const Vector3 &v, bool parallax);

// 太阳 J2000 球面坐标
//   t: 儒略世纪TD
//   n: 截断项数(精度); 推荐 20
Vector3 sun2000(long double t, int n);

} // namespace star_eph
