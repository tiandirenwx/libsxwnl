#include "map_projection.h"
#include "const.h"

#include <algorithm>
#include <cmath>

namespace map_projection
{

// 把弧度规约到 [-pi, pi]
static long double rad2rrad(long double v)
{
    const long double TWO_PI = static_cast<long double>(pi2);
    const long double P      = static_cast<long double>(PI);
    v = std::fmod(v, TWO_PI);
    if (v >  P) v -= TWO_PI;
    if (v < -P) v += TWO_PI;
    return v;
}

Projector::Projector()
    : tylx_(0), J0_(0), W0_(0), cosE0_(0), sinE0_(1),
      win_{0.0L, 0.0L, 1.0L, 1.0L}
{}

void Projector::setlx(int tylx, long double J0, long double W0, const Window &win)
{
    tylx_  = tylx;
    J0_    = J0;
    W0_    = W0;
    // 与 JS 中保持一致: cosE0 = sin(W0), sinE0 = cos(W0)
    cosE0_ = std::sin(W0);
    sinE0_ = std::cos(W0);
    win_   = win;
}

long double Projector::mollCZ(long double W) const
{
    constexpr int N = 100;
    const long double P_2 = static_cast<long double>(pi_2);
    long double f = 1.0L;
    if (W < 0) { W = -W; f = -1.0L; }
    if (W > P_2 - 1e-10L) return f;

    if (mollY_.empty())
    {
        mollY_.resize(N + 1);
        for (int i = 0; i < N; ++i)
        {
            long double y0 = 0.0L;
            long double y  = i / 100.0L;
            long double c  = P_2 * std::sin(y * P_2);
            while (std::fabs(y - y0) > 1e-12L)
            {
                y0 = y;
                y += (c - std::asin(y)) / std::sqrt(1.0L - y * y);
                y /= 2.0L;
            }
            mollY_[i] = y;
        }
        mollY_[N] = 1.0L;
    }

    long double n = W / P_2 * N;
    int k = (int)std::floor(n + 0.5L);
    if (k == 0) k = 1;
    if (k >= N) k = N - 1;
    n -= k;
    long double y2 = mollY_[k];
    long double a  = y2 - mollY_[k - 1];
    long double b  = mollY_[k + 1] - y2;
    return f * (y2 + n * (a + b + n * (b - a)) / 2.0L);
}

// 0: 平面正投(斜轴)
Point Projector::toxy0(long double J, long double W) const
{
    J -= J0_ + static_cast<long double>(pi_2);
    long double x = std::cos(W) * std::cos(J);
    long double y = std::cos(W) * std::sin(J);
    long double z = std::sin(W);
    Point p;
    p.x =  x;
    p.y =  cosE0_ * y + sinE0_ * z;
    p.z = -sinE0_ * y + cosE0_ * z;
    return p;
}

// 斜轴方位投影通用前置
static void azimuthal_pre(long double J0, long double cosE0, long double sinE0,
                          long double J, long double W,
                          long double &L, long double &B)
{
    J -= J0;
    long double cosJ = std::cos(J), sinJ = std::sin(J);
    long double cosW = std::cos(W), sinW = std::sin(W);
    L = std::atan2(sinE0 * sinW - cosE0 * cosW * cosJ, sinJ * cosW);
    B = std::acos(cosE0 * sinW + sinE0 * cosW * cosJ);
}

// 1: 斜轴等距方位
Point Projector::toxy1(long double J, long double W) const
{
    const long double P = static_cast<long double>(PI);
    long double L, B;
    azimuthal_pre(J0_, cosE0_, sinE0_, J, W, L, B);
    Point p{0, 0, -1};
    if (B > P - 0.1L) return p;
    B /= P;
    p.x = B * std::cos(L);
    p.y = B * std::sin(L);
    p.z = 1;
    return p;
}

// 2: 斜轴等积方位
Point Projector::toxy2(long double J, long double W) const
{
    const long double P = static_cast<long double>(PI);
    long double L, B;
    azimuthal_pre(J0_, cosE0_, sinE0_, J, W, L, B);
    Point p{0, 0, -1};
    if (B > P - 0.1L) return p;
    B = std::sin(B / 2.0L);
    p.x = B * std::cos(L);
    p.y = B * std::sin(L);
    p.z = 1;
    return p;
}

// 3: 斜轴等角方位
Point Projector::toxy3(long double J, long double W) const
{
    const long double P = static_cast<long double>(PI);
    long double L, B;
    azimuthal_pre(J0_, cosE0_, sinE0_, J, W, L, B);
    Point p{0, 0, -1};
    if (B > P - 1.0L) return p;
    B = std::tan(B / 2.0L);
    p.x = B * std::cos(L);
    p.y = B * std::sin(L);
    p.z = 1;
    return p;
}

// 4: 摩尔威特
Point Projector::toxy4(long double J, long double W) const
{
    const long double P_2 = static_cast<long double>(pi_2);
    W = rad2rrad(W);
    J = rad2rrad(J - J0_);
    Point p;
    p.y = mollCZ(W);
    p.x = std::sqrt(std::max<long double>(0.0L, 1.0L - p.y * p.y)) * J / P_2;
    p.z = 1;
    return p;
}

// 5: 正轴等距圆柱
Point Projector::toxy5(long double J, long double W) const
{
    const long double P = static_cast<long double>(PI);
    Point p;
    p.x = rad2rrad(J - J0_) / P;
    p.y = rad2rrad(W) / P;
    p.z = 1;
    return p;
}

// 6: 正轴等角圆柱
Point Projector::toxy6(long double J, long double W) const
{
    const long double P   = static_cast<long double>(PI);
    const long double P_4 = P / 4.0L;
    Point p;
    if (std::fabs(W) > 1.5L) { p.x = p.y = 0; p.z = -1; return p; }
    long double f = (W < 0) ? -1.0L : 1.0L;
    p.x = rad2rrad(J - J0_) / P;
    p.y = std::log(std::tan(P_4 + std::fabs(W) / 2.0L)) / P * f;
    p.z = 1;
    return p;
}

// 7: 多圆锥
Point Projector::toxy7(long double J, long double W) const
{
    const long double P = static_cast<long double>(PI);
    J = rad2rrad(J - J0_);
    long double C = std::cos(W);
    long double v = std::tan(W);
    Point p;
    if (std::fabs(W) > 0.002L)
    {
        long double s = C * J * v;
        long double l = 1.0L / v;
        p.x = l * std::sin(s);
        p.y = W + l - l * std::cos(s);
    }
    else
    {
        p.x = J;
        p.y = W + J * J * v / 2.0L;
    }
    p.x /= P;
    p.y /= P;
    p.z = 1;
    return p;
}

// 8: 中国灯笼投影
Point Projector::toxy8(long double J, long double W) const
{
    const long double P = static_cast<long double>(PI);
    J = rad2rrad(J - J0_);
    long double C = std::cos(W * 2.0L / 3.0L);
    long double v = std::sin(W * 2.0L / 3.0L) / 5.3L;
    J *= (1.0L - std::fabs(J) / 11.0L / P) * 1.1L;
    W = (W + 0.014830286L * W * W * W) * 1.2L;
    Point p;
    if (std::fabs(W) > 0.002L)
    {
        long double s = C * J * v;
        long double l = 1.0L / v;
        p.x = l * std::sin(s);
        p.y = W + l - l * std::cos(s);
    }
    else
    {
        p.x = J;
        p.y = W + J * J * v / 2.0L;
    }
    p.x /= P;
    p.y /= P;
    p.z = 1;
    return p;
}

Point Projector::toxy(long double J, long double W) const
{
    switch (tylx_)
    {
    case 0: return toxy0(J, W);
    case 1: return toxy1(J, W);
    case 2: return toxy2(J, W);
    case 3: return toxy3(J, W);
    case 4: return toxy4(J, W);
    case 5: return toxy5(J, W);
    case 6: return toxy6(J, W);
    case 7: return toxy7(J, W);
    case 8: return toxy8(J, W);
    default: return Point{J, W, 1};
    }
}

std::vector<long double> Projector::lineArr(const std::vector<long double> &d) const
{
    std::vector<long double> out;
    if (d.empty()) return out;
    out.reserve(d.size());

    long double x1 = win_.cx - win_.rx;
    long double x2 = win_.cx + win_.rx;
    long double y1 = win_.cy - win_.ry;
    long double y2 = win_.cy + win_.ry;

    bool        h  = (d[0] != kMoveTo); // 是否需要插入起点 moveto
    int         dd = 0;
    long double lx = 0, ly = 0;

    for (std::size_t i = 0; i < d.size(); ++i)
    {
        if (d[i] == kMoveTo)
        {
            h  = true;
            dd = 0;
            continue;
        }
        if (i + 1 >= d.size()) break;
        long double J = d[i];
        long double W = d[i + 1];
        ++i; // 跳过 W
        Point G = toxy(J, W);

        if (G.x < x1 || G.x > x2 || G.y < y1 || G.y > y2) G.z = -1;
        if (G.z < 0)
        {
            if (dd % 2 == 1) dd = 2;
            continue;
        }
        if (dd % 2 == 0) { ++dd; h = true; }

        if (tylx_ == 4)
        {
            if (std::fabs(lx - G.x) > 1.0L / 3.0L) h = true;
        }
        if (tylx_ == 5 || tylx_ == 6)
        {
            if (std::fabs(lx - G.x) > 1.0L / 3.0L ||
                std::fabs(ly - G.y) > 1.0L / 3.0L)
                h = true;
        }
        if (tylx_ == 7 || tylx_ == 8)
        {
            long double cx = std::fabs(lx - G.x);
            if (cx > 1.0L / 3.0L || (std::fabs(G.y) > 0.5L && lx * G.x < 0))
                h = true;
        }

        if (h)
        {
            out.push_back(kMoveTo);
            h = false;
        }
        out.push_back(G.x);
        out.push_back(G.y);
        lx = G.x;
        ly = G.y;
    }
    return out;
}

std::vector<long double> Projector::lineSingle(const long double *jw, std::size_t pairCount) const
{
    std::vector<long double> d;
    d.reserve(pairCount * 2);
    for (std::size_t i = 0; i < pairCount; ++i)
    {
        d.push_back(jw[i * 2]);
        d.push_back(jw[i * 2 + 1]);
    }
    return lineArr(d);
}

} // namespace map_projection
