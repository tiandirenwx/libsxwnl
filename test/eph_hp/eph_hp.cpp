// eph_hp.cpp
// ==========
// 参考实现: 基于 VSOP87B 全项计算地球日心黄经等分量.
//
// 与 src/eph.cpp 中 XL0_calc() 保持算法一致, 仅数据源从截断版 XL0_0
// 替换为 vsop87b_earth_hp.h 中的 XL0_0_HP (2564 项完整版).
// 章动、光行差直接复用 src/eph.h 中的公有函数 nutationLon2 / gxc_sunLon,
// 保证除 VSOP 那一层之外, 其余物理量与生产实现完全对齐.

#include "eph_hp.h"
#include "vsop87d_earth_hp.h"
#include "iau2000a_nut_lon_hp.h"

#include "../../src/eph.h"  // gxc_sunLon 声明 (光行差, 与章动无关, 复用)

#include <cmath>
#include <cstdint>

namespace sxwnl_hp
{

// 与 src/eph.cpp 内部常量对齐:
//   PI     : 圆周率
//   RAD_AS : 每弧度所含角秒 (180 * 3600 / PI)
// 之所以不 #include "src/eph.h" 里的 rad, 是因为它是 static 常量, 未导出.
static constexpr long double PI_HP = 3.14159265358979323846264338327950288L;
static constexpr long double RAD_AS = 180.0L * 3600.0L / PI_HP;

long double XL0_calc_HP(int zn, long double t, int n)
{
    // t: 儒略世纪数 (自 J2000). 内部转成千年数.
    t /= 10.0L;
    long double v = 0, tn = 1, c = 0;
    const long double *F = XL0_0_HP;
    long double n1, n2, N, n0;

    // 位置索引表的偏移 (与 src/eph.cpp XL0_calc 完全一致):
    //   pn = zn * 6 + 1, zn: 0=L, 1=B, 2=R
    int pn = zn * 6 + 1;
    long double N0 = F[pn + 1] - F[pn];  // 该变量 T^0 段的项数三倍数
    for (int i = 0; i < 6; i++, tn *= t)
    {
        n1 = F[pn + i];
        n2 = F[pn + 1 + i];
        n0 = n2 - n1;
        if (!n0)
        {
            continue;
        }
        if (n < 0)
        {
            N = n2;  // 使用该阶所有项
        }
        else
        {
            // 自适应截断: 按项数 n 按比例分配到各时间幂
            // (与 XL0_calc 相同的启发式)
            N = (long long)(3 * n * n0 / N0 + 0.5) + n1;
            if (i)
            {
                N += 3;
            }
            if (N > n2)
            {
                N = n2;
            }
        }
        c = 0;
        for (long long j = (long long)n1; j < (long long)N; j += 3)
        {
            c += F[j] * std::cos((double)(F[j + 1] + t * F[j + 2]));
        }
        v += c * tn;
    }
    v /= F[0];  // 除以振幅倍率 (10^10)

    // 地球特有的长期项修正 (来自 src/eph.cpp XL0_calc xt==0 分支)
    // 这些修正项和 VSOP87 数据本身一起构成 IMCCE 完整方案.
    long double t2 = t * t;
    long double t3 = t2 * t;
    if (zn == 0)
    {
        v += (-0.0728L - 2.7702L * t - 1.1019L * t2 - 0.0996L * t3) / RAD_AS;
    }
    else if (zn == 1)
    {
        v += (+0.0000L + 0.0004L * t + 0.0004L * t2 - 0.0026L * t3) / RAD_AS;
    }
    else if (zn == 2)
    {
        v += (-0.0020L + 0.0044L * t + 0.0213L * t2 - 0.0250L * t3) / 1000000.0L;
    }
    return v;
}

long double E_Lon_HP(long double t, int n)
{
    return XL0_calc_HP(0, t, n);
}

// ============================================================================
// IAU 2000A 全项章动黄经 (基于 IERS Conv 2010, §5.5 + §5.7.2)
// 参考 SOFA 库 nut00a.c 的实现方式 (公有开源, 数值验证过多个平台)
// ============================================================================

// 常量: 弧秒 → 弧度
static constexpr long double DAS2R_HP = 4.848136811095359935899141e-6L;  // π/(180*3600)
// 常量: 微弧秒 → 弧度 (μas)
static constexpr long double DMAS2R_HP = DAS2R_HP / 1.0e6L;
// 常量: 2π
static constexpr long double D2PI_HP = 2.0L * PI_HP;

// 把弧度归化到 [0, 2π)
static inline long double fmod2pi(long double x)
{
    x = std::fmod((double)x, (double)D2PI_HP);
    if (x < 0) x += D2PI_HP;
    return x;
}

// 五个 Delaunay 基本参数 (l, l', F, D, Ω)
// 传入 t: 儒略世纪数 (自 J2000 TT)
// 返回值: 弧度
//
// 数值来源: SOFA v18 nut00a.c (IERS Conv 2010 §5.7.2)
// 关键技巧: 每个基本参数每世纪都完成大量整数圈, 直接展开会丢精度.
//   如 Moon 平近点角每世纪走 1325 圈整数 + "余项" 715923.2178 arcsec.
//   分开处理: 整数圈用 fmod(N·t, 1) · 2π, 余项用多项式.
static inline long double delaunay_l(long double t)
{
    long double as = ((485868.249036L
                      + (715923.2178L
                      + (31.8792L
                      + (0.051635L
                      + (-0.00024470L)
                      * t) * t) * t) * t));
    return std::fmod(as * DAS2R_HP + std::fmod(1325.0L * t, 1.0L) * D2PI_HP,
                     D2PI_HP);
}

static inline long double delaunay_lp(long double t)
{
    long double as = ((1287104.79305L
                      + (1292581.0481L
                      + (-0.5532L
                      + (0.000136L
                      + (-0.00001149L)
                      * t) * t) * t) * t));
    return std::fmod(as * DAS2R_HP + std::fmod(99.0L * t, 1.0L) * D2PI_HP,
                     D2PI_HP);
}

static inline long double delaunay_F(long double t)
{
    long double as = ((335779.526232L
                      + (295262.8478L
                      + (-12.7512L
                      + (-0.001037L
                      + (0.00000417L)
                      * t) * t) * t) * t));
    return std::fmod(as * DAS2R_HP + std::fmod(1342.0L * t, 1.0L) * D2PI_HP,
                     D2PI_HP);
}

static inline long double delaunay_D(long double t)
{
    long double as = ((1072260.70369L
                      + (1105601.2090L
                      + (-6.3706L
                      + (0.006593L
                      + (-0.00003169L)
                      * t) * t) * t) * t));
    return std::fmod(as * DAS2R_HP + std::fmod(1236.0L * t, 1.0L) * D2PI_HP,
                     D2PI_HP);
}

static inline long double delaunay_Om(long double t)
{
    long double as = ((450160.398036L
                      + (-482890.5431L
                      + (7.4722L
                      + (0.007702L
                      + (-0.00005939L)
                      * t) * t) * t) * t));
    return std::fmod(as * DAS2R_HP + std::fmod(-5.0L * t, 1.0L) * D2PI_HP,
                     D2PI_HP);
}

// 8 行星平黄经 (IERS Conv 2010, 弧度直接给出线性系数)
static inline long double planet_lon_Me(long double t) { return fmod2pi(4.402608842L + 2608.7903141574L * t); }
static inline long double planet_lon_Ve(long double t) { return fmod2pi(3.176146697L + 1021.3285546211L * t); }
static inline long double planet_lon_E (long double t) { return fmod2pi(1.753470314L +  628.3075849991L * t); }
static inline long double planet_lon_Ma(long double t) { return fmod2pi(6.203480913L +  334.0612426700L * t); }
static inline long double planet_lon_J (long double t) { return fmod2pi(0.599546497L +   52.9690962641L * t); }
static inline long double planet_lon_Sa(long double t) { return fmod2pi(0.874016757L +   21.3299104960L * t); }
static inline long double planet_lon_U (long double t) { return fmod2pi(5.481293872L +    7.4781598567L * t); }
static inline long double planet_lon_Ne(long double t) { return fmod2pi(5.311886287L +    3.8133035638L * t); }

// 岁差累积 p_A
static inline long double precession_pA(long double t)
{
    return (0.024381750L + 0.00000538691L * t) * t;
}

long double nutationLon_HP(long double t)
{
    // 14 个基本参数 (radians)
    long double F[14];
    F[0]  = delaunay_l(t);
    F[1]  = delaunay_lp(t);
    F[2]  = delaunay_F(t);
    F[3]  = delaunay_D(t);
    F[4]  = delaunay_Om(t);
    F[5]  = planet_lon_Me(t);
    F[6]  = planet_lon_Ve(t);
    F[7]  = planet_lon_E(t);
    F[8]  = planet_lon_Ma(t);
    F[9]  = planet_lon_J(t);
    F[10] = planet_lon_Sa(t);
    F[11] = planet_lon_U(t);
    F[12] = planet_lon_Ne(t);
    F[13] = precession_pA(t);

    // 累加 j=0 段 (常数振幅)
    // 每项 16 个 long double: [A_μas, App_μas, N_0..N_13]
    long double dpsi_j0 = 0;
    const int COLS = NUT_LON_HP_COLS;
    for (int i = 0; i < NUT_LON_HP_J0_COUNT; ++i)
    {
        const long double *row = &NUT_LON_HP_J0[i * COLS];
        long double A   = row[0];
        long double App = row[1];
        long double arg = 0;
        for (int k = 0; k < 14; ++k)
        {
            long double n = row[2 + k];
            if (n != 0.0L)
            {
                arg += n * F[k];
            }
        }
        dpsi_j0 += A * std::sin((double)arg) + App * std::cos((double)arg);
    }

    // 累加 j=1 段 (乘以 t)
    long double dpsi_j1 = 0;
    for (int i = 0; i < NUT_LON_HP_J1_COUNT; ++i)
    {
        const long double *row = &NUT_LON_HP_J1[i * COLS];
        long double A   = row[0];
        long double App = row[1];
        long double arg = 0;
        for (int k = 0; k < 14; ++k)
        {
            long double n = row[2 + k];
            if (n != 0.0L)
            {
                arg += n * F[k];
            }
        }
        dpsi_j1 += A * std::sin((double)arg) + App * std::cos((double)arg);
    }

    // 结果: 单位 μas -> 弧度
    return (dpsi_j0 + t * dpsi_j1) * DMAS2R_HP;
}

long double S_aLon_HP(long double t, int n)
{
    // S_aLon = E_Lon (地球日心黄经) + 章动 + 光行差 + π
    // 相较 XL::S_aLon: E_Lon 用 VSOP87D 全项, 章动用 IAU 2000A 全项, 光行差沿用中精度
    return E_Lon_HP(t, n) + nutationLon_HP(t) + gxc_sunLon(t) + PI_HP;
}

long double S_aLon_t_HP(long double W)
{
    // 与 src/eph.cpp XL::S_aLon_t 相同的两步 Newton 迭代:
    //   1. 用地球平均角速度 v 做初估
    //   2. 中精度 + 全精度 两轮修正
    long double v = 628.3319653318L;
    long double t = (W - 1.75347L - PI_HP) / v;
    t += (W - S_aLon_HP(t, 10)) / v;
    v = 628.332L + 21.0L * std::sin(1.527L + 628.3076L * t)
        + 0.44L * std::sin(1.48L + 1256.6152L * t);
    t += (W - S_aLon_HP(t, -1)) / v;
    return t;
}

}  // namespace sxwnl_hp
