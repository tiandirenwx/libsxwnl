// ─────────────────────────────────────────────────────────────────────────
//  ld_compat.h — long double libm 兼容垫片 (通过 CMake -include 全局强制包含)
//
//  背景:
//  arm64 的 `long double` 是 128 位 binary128。实测发现部分平台 (安卓 bionic,
//  鸿蒙同架构) 对 binary128 的 long double 版数学函数存在缺陷 —— 例如
//  `sinl(-94.9)` 会返回 0 (而 `cosl`、double 版 `sin` 均正常)。这会击穿天文
//  级数计算 (llrConv → pty_zty2 → 真太阳时), 使八字排盘在勾选"真太阳时"时
//  日/时柱整体错约 8 小时。
//
//  修复:
//  double (64 位) 精度对本库完全足够 (原版 sxwnl 即 JS double 精度; host/iOS
//  的 long double 本就是 64 位且结果全部正确)。故把 long double 版 libm 统一
//  路由到可靠的 double 版实现, 彻底绕开 binary128 缺陷。
//
//  判据平台无关: 仅当 long double 为 128 位 binary128 (LDBL_MANT_DIG > 64) 时生效,
//  精确命中 arm64 (安卓 / 鸿蒙 / Linux)。x86 的 80 位扩展精度 (libm 正常) 与
//  苹果 / 32 位 ARM 的 64 位 long double 均不受影响。跨平台 src/ 源码零改动。
// ─────────────────────────────────────────────────────────────────────────
#pragma once

// 先纳入标准数学头, 让真实函数声明先于宏定义被处理, 避免宏破坏 <cmath> 声明。
#include <cfloat>
#include <cmath>
#include <math.h>

#if (LDBL_MANT_DIG > 64)

#define sinl(x)      ((long double)::sin((double)(x)))
#define cosl(x)      ((long double)::cos((double)(x)))
#define tanl(x)      ((long double)::tan((double)(x)))
#define asinl(x)     ((long double)::asin((double)(x)))
#define acosl(x)     ((long double)::acos((double)(x)))
#define atanl(x)     ((long double)::atan((double)(x)))
#define atan2l(y, x) ((long double)::atan2((double)(y), (double)(x)))
#define fmodl(a, b)  ((long double)::fmod((double)(a), (double)(b)))
#define sqrtl(x)     ((long double)::sqrt((double)(x)))
#define fabsl(x)     ((long double)::fabs((double)(x)))
#define powl(a, b)   ((long double)::pow((double)(a), (double)(b)))
#define expl(x)      ((long double)::exp((double)(x)))
#define logl(x)      ((long double)::log((double)(x)))
#define log10l(x)    ((long double)::log10((double)(x)))
#define floorl(x)    ((long double)::floor((double)(x)))
#define ceill(x)     ((long double)::ceil((double)(x)))
#define truncl(x)    ((long double)::trunc((double)(x)))
#define roundl(x)    ((long double)::round((double)(x)))
#define hypotl(a, b) ((long double)::hypot((double)(a), (double)(b)))

#endif // LDBL_MANT_DIG > 64
