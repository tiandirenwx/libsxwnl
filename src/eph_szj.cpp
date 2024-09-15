#include "eph_szj.h"
#include "eph.h"
#include "JD.h"

SunMoonRiseSet::SunMoonRiseSet()
{
    L = 0;
    fa = 0;
    dt = 0;
    E = 0.409092614;
}
SunMoonRiseSet::~SunMoonRiseSet()
{
}

long double SunMoonRiseSet::getH(long double h, long double w)
{ // h地平纬度,w赤纬,返回时角
    long double c = (sinl(h) - sinl(fa) * sinl(w)) / cosl(fa) / cosl(w);
    if (fabsl(c) > 1)
    {
        return M_PI;
    }
    return acosl(c);
};

void SunMoonRiseSet::Mcoord(long double jd, long double H0, SJ &r)
{                                                        // 章动同时影响恒星时和天体坐标,所以不计算章动。返回时角及赤经纬
    Vector3 z = m_coord((jd + dt) / 36525.0, 40, 30, 8); // 低精度月亮赤经纬
    z = llrConv(z, E);                                   // 转为赤道坐标
    r.H = rad2rrad(pGST(jd, dt) + L - z[0]);             // 得到此刻天体时角
    if (H0)
    {
        r.H0 = getH(0.7275 * cs_rEar / z[2] - 34 * 60 / rad, z[1]); // 升起对应的时角
    }
}

SJ SunMoonRiseSet::Mt(long double jd)
{ // 月亮到中升降时刻计算,传入jd含义与St()函数相同
    dt = dt_T(jd);
    E = hcjj(jd / 36525.0);
    jd -= mod2(0.1726222 + 0.966136808032357 * jd - 0.0366 * dt + L / pi2, 1); // 查找最靠近当日中午的月上中天,mod2的第1参数为本地时角近似值

    SJ r = {};
    long double sv = pi2 * 0.966;
    r.z = r.x = r.s = r.j = r.c = r.h = jd;
    Mcoord(jd, 1, r); // 月亮坐标
    r.s += (-r.H0 - r.H) / sv;
    r.j += (r.H0 - r.H) / sv;
    r.z += (0 - r.H) / sv;
    r.x += (M_PI - r.H) / sv;
    Mcoord(r.s, 1, r);
    r.s += rad2rrad(-r.H0 - r.H) / sv;
    Mcoord(r.j, 1, r);
    r.j += rad2rrad(+r.H0 - r.H) / sv;
    Mcoord(r.z, 0, r);
    r.z += rad2rrad(0 - r.H) / sv;
    Mcoord(r.x, 0, r);
    r.x += rad2rrad(M_PI - r.H) / sv;
    return r;
}

void SunMoonRiseSet::Scoord(long double jd, int xm, SJ &r)
{                                                                              // 章动同时影响恒星时和天体坐标,所以不计算章动。返回时角及赤经纬
    Vector3 z = {XL::E_Lon((jd + dt) / 36525.0, 5) + M_PI - 20.5 / rad, 0, 1}; // 太阳坐标(修正了光行差)
    z = llrConv(z, E);                                                         // 转为赤道坐标
    r.H = rad2rrad(pGST(jd, dt) + L - z[0]);                                   // 得到此刻天体时角

    if (xm == 10 || xm == 1)
    {
        r.H1 = getH(-50 * 60 / rad, z[1]); // 地平以下50分
    }
    if (xm == 10 || xm == 2)
    {
        r.H2 = getH(-6 * 3600 / rad, z[1]); // 地平以下6度
    }
    if (xm == 10 || xm == 3)
    {
        r.H3 = getH(-12 * 3600 / rad, z[1]); // 地平以下12度
    }
    if (xm == 10 || xm == 4)
    {
        r.H4 = getH(-18 * 3600 / rad, z[1]); // 地平以下18度
    }
}

SJ SunMoonRiseSet::St(long double jd)
{ // 太阳到中升降时刻计算,传入jd是当地中午12点时间对应的2000年首起算的格林尼治时间UT
    dt = dt_T(jd);
    E = hcjj(jd / 36525.0);
    jd -= mod2(jd + L / pi2, 1); // 查找最靠近当日中午的日上中天,mod2的第1参数为本地时角近似值

    SJ r = {};
    long double sv = pi2;
    r.z = r.x = r.s = r.j = r.c = r.h = r.c2 = r.h2 = r.c3 = r.h3 = jd;
    r.sm = "";
    Scoord(jd, 10, r);         // 太阳坐标
    r.s += (-r.H1 - r.H) / sv; // 升起
    r.j += (r.H1 - r.H) / sv;  // 降落

    r.c += (-r.H2 - r.H) / sv;  // 民用晨
    r.h += (r.H2 - r.H) / sv;   // 民用昏
    r.c2 += (-r.H3 - r.H) / sv; // 航海晨
    r.h2 += (r.H3 - r.H) / sv;  // 航海昏
    r.c3 += (-r.H4 - r.H) / sv; // 天文晨
    r.h3 += (r.H4 - r.H) / sv;  // 天文昏

    r.z += (0 - r.H) / sv;    // 中天
    r.x += (M_PI - r.H) / sv; // 下中天
    Scoord(r.s, 1, r);
    r.s += rad2rrad(-r.H1 - r.H) / sv;
    if (areAlmostEqual(r.H1, M_PI))
    {
        r.sm += "无升起.";
    }
    Scoord(r.j, 1, r);
    r.j += rad2rrad(+r.H1 - r.H) / sv;
    if (areAlmostEqual(r.H1, M_PI))
    {
        r.sm += "无降落.";
    }

    Scoord(r.c, 2, r);
    r.c += rad2rrad(-r.H2 - r.H) / sv;
    if (areAlmostEqual(r.H2, M_PI))
    {
        r.sm += "无民用晨.";
    }
    Scoord(r.h, 2, r);
    r.h += rad2rrad(+r.H2 - r.H) / sv;
    if (areAlmostEqual(r.H2, M_PI))
    {
        r.sm += "无民用昏.";
    }
    Scoord(r.c2, 3, r);
    r.c2 += rad2rrad(-r.H3 - r.H) / sv;
    if (areAlmostEqual(r.H3, M_PI))
    {
        r.sm += "无航海晨.";
    }
    Scoord(r.h2, 3, r);
    r.h2 += rad2rrad(+r.H3 - r.H) / sv;
    if (areAlmostEqual(r.H3, M_PI))
    {
        r.sm += "无航海昏.";
    }
    Scoord(r.c3, 4, r);
    r.c3 += rad2rrad(-r.H4 - r.H) / sv;
    if (areAlmostEqual(r.H4, M_PI))
    {
        r.sm += "无天文晨.";
    }
    Scoord(r.h3, 4, r);
    r.h3 += rad2rrad(+r.H4 - r.H) / sv;
    if (areAlmostEqual(r.H4, M_PI))
    {
        r.sm += "无天文昏.";
    }

    Scoord(r.z, 0, r);
    r.z += (0 - r.H) / sv;
    Scoord(r.x, 0, r);
    r.x += rad2rrad(M_PI - r.H) / sv;
    return r;
}

void SunMoonRiseSet::calcRTS(long double jd, int n, long double Jdl, long double Wdl, long double sq)
{ // 多天升中降计算,jd是当地起始略日(中午时刻),sq是时区
    int i;
    int64_t c;
    SJ_S rr;
    SJ r;
    if (!rts.size())
    {
        for (i = 0; i < 31; i++)
        {
            rts.push_back({});
        }
    }
    L = Jdl, fa = Wdl, sq /= 24.0; // 设置站点参数
    for (i = 0; i < n; i++)
    {
        rr = rts[i];
        rr.Ms = rr.Mz = rr.Mj = "--:--:--";
    }
    for (i = -1; i <= n; i++)
    {
        if (i >= 0 && i < n)
        { // 太阳
            r = St(jd + i + sq);
            rts[i].s = JD::timeStr(r.s - sq);         // 升
            rts[i].z = JD::timeStr(r.z - sq);         // 中
            rts[i].j = JD::timeStr(r.j - sq);         // 降
            rts[i].c = JD::timeStr(r.c - sq);         // 晨
            rts[i].h = JD::timeStr(r.h - sq);         // 昏
            rts[i].ch = JD::timeStr(r.h - r.c - 0.5); // 光照时间,timeStr()内部+0.5,所以这里补上-0.5
            rts[i].sj = JD::timeStr(r.j - r.s - 0.5); // 昼长
        }
        r = Mt(jd + i + sq); // 月亮
        c = int64(r.s - sq + 0.5) - jd;
        if (c >= 0 && c < n)
        {
            rts[c].Ms = JD::timeStr(r.s - sq);
        }
        c = int64(r.z - sq + 0.5) - jd;
        if (c >= 0 && c < n)
        {
            rts[c].Mz = JD::timeStr(r.z - sq);
        }
        c = int64(r.j - sq + 0.5) - jd;
        if (c >= 0 && c < n)
        {
            rts[c].Mj = JD::timeStr(r.j - sq);
        }
    }
}

void SunMoonRiseSet::setE(long double x)
{
    E = x;
}

long double SunMoonRiseSet::getE() const
{
    return E;
}

void SunMoonRiseSet::setDt(long double x)
{
    dt = x;
}

long double SunMoonRiseSet::getDt() const
{
    return dt;
}

void SunMoonRiseSet::setL(long double x)
{
    L = x;
}

long double SunMoonRiseSet::getL() const
{
    return L;
}

void SunMoonRiseSet::setFa(long double x)
{
    fa = x;
}

long double SunMoonRiseSet::getFa() const
{
    return fa;
}
