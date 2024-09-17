#include "eph_rspl.h"
#include "eph_rsgs.h"
#include <iostream>

EphRspl::EphRspl()
{
    pstSeMax = new SE_MAX_DATA();
    pstNbj = new SE_NBJ_DATA();
    pstSeMax->nasa_r = false; // 为1表示采用NASA的视径比

    pstSeMax->sT = {}; // 地方日食时间表
    pstSeMax->LX = "";
    pstSeMax->sf = 0;
    pstSeMax->sf2 = 0;
    pstSeMax->sf3 = 0;
    pstSeMax->sflx = "";
    pstSeMax->b1 = 0;
    pstSeMax->dur = 0;
    pstSeMax->sun_s = 0;
    pstSeMax->sun_j = 0;
    pstSeMax->P1 = 0;
    pstSeMax->V1 = 0;
    pstSeMax->P2 = 0;
    pstSeMax->V2 = 0;

    // 以下计算南北界用到
    pstNbj->A = {};
    pstNbj->B = {}; // 本半影锥顶点坐标
    pstNbj->P = {}; // t1时刻的日月坐标,g为恒星时
    pstNbj->Q = {}; // t2时刻的日月坐标
    pstNbj->V = {}; // 食界表
    pstNbj->Vc = "";
    pstNbj->Vb = ""; // 食中心类型,本影南北距离*/
}

EphRspl::~EphRspl()
{
    delete pstSeMax;
    delete pstNbj;
}

void EphRspl::secMax(long double jd, long double L, long double fa, long double high)
{
    // 日食的食甚计算(jd为近朔的力学时,误差几天不要紧)
    int i;
    for (i = 0; i < 5; i++)
    {
        pstSeMax->sT[i] = 0; // 分别是:食甚,初亏,复圆,食既,生光
    }
    pstSeMax->LX = "";                     // 类型
    pstSeMax->sf = 0;                      // 食分
    pstSeMax->b1 = 1.0;                    // 月日半径比(食甚时刻)
    pstSeMax->dur = 0;                     // 持续时间
    pstSeMax->P1 = pstSeMax->V1 = 0;       // 初亏方位,P北点起算,V顶点起算
    pstSeMax->P2 = pstSeMax->V2 = 0;       // 复圆方位,P北点起算,V顶点起算
    pstSeMax->sun_s = pstSeMax->sun_j = 0; // 日出日没
    pstSeMax->sf2 = 0;                     // 食分(日出食分)
    pstSeMax->sf3 = 0;                     // 食分(日没食分)
    pstSeMax->sflx = " ";                  // 食分类型
    EphRsgs &rsgsObj = EphRsgs::getInstance();
    rsgsObj.init(jd, 7);
    auto sp = rsgsObj.getRsgsParaData();
    jd = sp.getZjd(); // 食甚初始估值为插值表中心时刻(粗朔)

    _SECXY G = {}, g = {};
    secXY(jd, L, fa, high, G);
    jd -= G.x / 0.2128; // 与食甚的误差在20分钟以内

    long double u, v, dt = 60 / 86400.0, dt2;
    for (i = 0; i < 2; i++)
    {
        secXY(jd, L, fa, high, G);
        secXY(jd + dt, L, fa, high, g);
        u = (g.y - G.y) / dt;
        v = (g.x - G.x) / dt;
        dt2 = -(G.y * u + G.x * v) / (u * u + v * v);
        jd += dt2; // 极值时间
    }

    // 求直线到太阳中心的最小值
    // double x=G.x+dt2*v, y=G.y+dt2*u, rmin=sqrt(x*x+y*y);
    // js v5.10更新计算公式
    long double maxsf = 0, maxjd = jd, rmin, ls, tt;
    for (i = -30; i < 30; i += 6)
    {
        tt = jd + i / 86400.0;
        secXY(tt, L, fa, high, g);
        ls = (g.mr + g.sr - sqrtl(g.x * g.x + g.y * g.y)) / g.sr / 2.0;
        if (ls > maxsf)
        {
            maxsf = ls, maxjd = tt;
        }
    }
    jd = maxjd;
    for (i = -5; i < 5; i += 1)
    {
        tt = jd + i / 86400.0;
        secXY(tt, L, fa, high, g);
        ls = (g.mr + g.sr - sqrtl(g.x * g.x + g.y * g.y)) / g.sr / 2;
        if (ls > maxsf)
        {
            maxsf = ls, maxjd = tt;
        }
    }
    jd = maxjd;
    secXY(jd, L, fa, high, G);
    rmin = sqrtl(G.x * G.x + G.y * G.y);

    pstSeMax->sun_s = sunShengJ(jd - dt_T(jd) + L / pi2, L, fa, -1) + dt_T(jd); // 日出 统一用力学时
    pstSeMax->sun_j = sunShengJ(jd - dt_T(jd) + L / pi2, L, fa, 1) + dt_T(jd);  // 日没 统一用力学时

    if (rmin <= G.mr + G.sr)
    {                         // 偏食计算
        pstSeMax->sT[1] = jd; // 食甚
        pstSeMax->LX = "偏";
        pstSeMax->sf = (G.mr + G.sr - rmin) / G.sr / 2.0; // 食分
        pstSeMax->b1 = G.mr / G.sr;

        secXY(pstSeMax->sun_s, L, fa, high, g);                                   // 日出食分
        pstSeMax->sf2 = (g.mr + g.sr - sqrtl(g.x * g.x + g.y * g.y)) / g.sr / 2.0; // 日出食分
        if (pstSeMax->sf2 < 0)
        {
            pstSeMax->sf2 = 0;
        }

        secXY(pstSeMax->sun_j, L, fa, high, g);                                 // 日没食分
        pstSeMax->sf3 = (g.mr + g.sr - sqrtl(g.x * g.x + g.y * g.y)) / g.sr / 2; // 日没食分
        if (pstSeMax->sf3 < 0)
        {
            pstSeMax->sf3 = 0;
        }

        pstSeMax->sT[0] = lineT(G, v, u, G.mr + G.sr, false); // 初亏
        for (i = 0; i < 3; i++)
        { // 初亏再算3次
            secXY(pstSeMax->sT[0], L, fa, high, g);
            pstSeMax->sT[0] = lineT(g, v, u, g.mr + g.sr, false);
        }

        pstSeMax->P1 = rad2mrad(atan2l(g.x, g.y));                                                     // 初亏位置角
        pstSeMax->V1 = rad2mrad(pstSeMax->P1 - shiChaJ(pGST2(pstSeMax->sT[0]), L, fa, g.sCJ, g.sCW)); // 这里g.sCJ与g.sCW对应的时间与sT[0]还差了一点，所以有一小点误差，不采用真恒星时也误差一点

        pstSeMax->sT[2] = lineT(G, v, u, G.mr + G.sr, true); // 复圆
        for (i = 0; i < 3; i++)
        { // 复圆再算3次
            secXY(pstSeMax->sT[2], L, fa, high, g);
            pstSeMax->sT[2] = lineT(g, v, u, g.mr + g.sr, true);
        }
        pstSeMax->P2 = rad2mrad(atan2l(g.x, g.y));
        pstSeMax->V2 = rad2mrad(pstSeMax->P2 - shiChaJ(pGST2(pstSeMax->sT[2]), L, fa, g.sCJ, g.sCW)); // 这里g.sCJ与g.sCW对应的时间与sT[2]还差了一点，所以有一小点误差，不采用真恒星时也误差一点
    }
    if (rmin <= G.mr - G.sr)
    { // 全食计算
        pstSeMax->LX = "全";
        pstSeMax->sT[3] = lineT(G, v, u, G.mr - G.sr, false); // 食既
        secXY(pstSeMax->sT[3], L, fa, high, g);
        pstSeMax->sT[3] = lineT(g, v, u, g.mr - g.sr, false); // 食既再算1次

        pstSeMax->sT[4] = lineT(G, v, u, G.mr - G.sr, true); // 生光
        secXY(pstSeMax->sT[4], L, fa, high, g);
        pstSeMax->sT[4] = lineT(g, v, u, g.mr - g.sr, true); // 生光再算1次
        pstSeMax->dur = pstSeMax->sT[4] - pstSeMax->sT[3];
    }
    if (rmin <= G.sr - G.mr)
    { // 环食计算
        pstSeMax->LX = "环";
        pstSeMax->sT[3] = lineT(G, v, u, G.sr - G.mr, false); // 食既
        secXY(pstSeMax->sT[3], L, fa, high, g);
        pstSeMax->sT[3] = lineT(g, v, u, g.sr - g.mr, false); // 食既再算1次

        pstSeMax->sT[4] = lineT(G, v, u, G.sr - G.mr, true); // 生光
        secXY(pstSeMax->sT[4], L, fa, high, g);
        pstSeMax->sT[4] = lineT(g, v, u, g.sr - g.mr, true); // 生光再算1次
        pstSeMax->dur = pstSeMax->sT[4] - pstSeMax->sT[3];
    }

    if (pstSeMax->sT[1] < pstSeMax->sun_s && pstSeMax->sf2 > 0)
    {
        pstSeMax->sf = pstSeMax->sf2, pstSeMax->sflx = "#"; // 食甚在日出前，取日出食分
    }
    if (pstSeMax->sT[1] > pstSeMax->sun_j && pstSeMax->sf3 > 0)
    {
        pstSeMax->sf = pstSeMax->sf3, pstSeMax->sflx = "*"; // 食甚在日没后，取日没食分
    }

    for (i = 0; i < 5; i++)
    {
        if (pstSeMax->sT[i] < pstSeMax->sun_s || pstSeMax->sT[i] > pstSeMax->sun_j)
        {
            pstSeMax->sT[i] = 0; // 升降时间之外的日食算值无效，因为地球不是透明的
            //std::cout<< "st value: "<<i<<std::endl;
        }
    }

    pstSeMax->sun_s -= dt_T(jd);
    pstSeMax->sun_j -= dt_T(jd);
}

void EphRspl::nbj(long double jd)
{
    // 南北界计算
    EphRsgs &rsgsObj = EphRsgs::getInstance();
    rsgsObj.init(jd, 7);

    _GJW G = {};
    std::array<long double, 10> &V = pstNbj->V;
    int i;
    for (i = 0; i < 10; i++)
    {
        V[i] = 100;
    }
    pstNbj->Vc = "", pstNbj->Vb = ""; // 返回初始化,纬度值为100表示无解,经度100也是无解,但在以下程序中经度会被转为-PI到+PI

    zb0(jd);
    pp0(G);
    V[0] = G.J, V[1] = G.W, pstNbj->Vc = G.c; // 食中心

    G.J = G.W = 0;
    for (i = 0; i < 2; i++)
    {
        p2p(G.J, G.W, G, true, 1);
    }
    V[2] = G.J, V[3] = G.W; // 本影北界,环食为南界(本影区之内,变差u,v基本不变,所以计算两次足够)
    G.J = G.W = 0;
    for (i = 0; i < 2; i++)
    {
        p2p(G.J, G.W, G, true, -1);
    }
    V[4] = G.J, V[5] = G.W; // 本影南界,环食为北界
    G.J = G.W = 0;
    for (i = 0; i < 3; i++)
    {
        p2p(G.J, G.W, G, false, -1);
    }
    V[6] = G.J, V[7] = G.W; // 半影北界
    G.J = G.W = 0;
    for (i = 0; i < 3; i++)
    {
        p2p(G.J, G.W, G, false, 1);
    }
    V[8] = G.J, V[9] = G.W; // 半影南界

    if (V[3] != 100 && V[5] != 100)
    { // 粗算本影南北距离
        long double x = (V[2] - V[4]) * cosl((V[3] + V[5]) / 2.), y = V[3] - V[5];
        pstNbj->Vb = toFixed(cs_rEarA * sqrtl(x * x + y * y), 0) + "千米";
    }
}

long double EphRspl::lineT(_SECXY G, long double v, long double u, long double r, bool n)
{
    // 已知t1时刻星体位置、速度，求x*x+y*y=r*r时,t的值
    long double b = G.y * v - G.x * u;
    long double A = u * u + v * v;
    long double B = u * b;
    long double C = b * b - r * r * v * v;
    long double D = B * B - A * C;
    if (D < 0)
    {
        return 0;
    }
    D = sqrtl(D);
    if (!n)
    {
        D = -D;
    }
    return G.t + ((-B + D) / A - G.x) / v;
}

void EphRspl::zbXY(_ZB &p, long double L, long double fa)
{

    Vector3 s = {p.S[0], p.S[1], p.S[2]};
    Vector3 m = {p.M[0], p.M[1], p.M[2]};
    s = parallax(s, p.g + L - p.S[0], fa, 0); // 修正了视差的赤道坐标
    m = parallax(m, p.g + L - p.M[0], fa, 0); // 修正了视差的赤道坐标
    //=======视半径========
    p.mr = cs_sMoon / m[2] / rad;
    p.sr = 959.63 / s[2] / rad * cs_AU;
    //=======日月赤经纬差转为日面中心直角坐标(用于日食)==============
    p.x = rad2rrad(m[0] - s[0]) * cosl((m[1] + s[1]) / 2.0);
    p.y = m[1] - s[1];
}

void EphRspl::zb0(long double jd)
{
    // 基本参数计算
    auto deltat = dt_T(jd); // TD-UT
    auto E = hcjj(jd / 36525.0);
    auto zd = nutation2(jd / 36525.0);

    EphRsgs &rsgsObj = EphRsgs::getInstance();

    pstNbj->P.g = pGST(jd - deltat, deltat) + zd[0] * cosl(E + zd[1]); // 真恒星时(不考虑非多项式部分)
    pstNbj->P.S = rsgsObj.sun(jd);
    pstNbj->P.M = rsgsObj.moon(jd);

    long double t2 = jd + 60 / 86400.0;
    pstNbj->Q.g = pGST(t2 - deltat, deltat) + zd[0] * cosl(E + zd[1]);
    pstNbj->Q.S = rsgsObj.sun(t2);
    pstNbj->Q.M = rsgsObj.moon(t2);

    // 转为直角坐标
    Vector3 z1 = {}, z2 = {};
    z1 = llr2xyz(pstNbj->P.S);
    z2 = llr2xyz(pstNbj->P.M);

    long double k = 959.63 / cs_sMoon * cs_AU; // k为日月半径比
    // 本影锥顶点坐标计算
    Vector3 F = {
        (z1[0] - z2[0]) / (1 - k) + z2[0],
        (z1[1] - z2[1]) / (1 - k) + z2[1],
        (z1[2] - z2[2]) / (1 - k) + z2[2]};
    pstNbj->A = xyz2llr(F);
    // 半影锥顶点坐标计算
    Vector3 FF = {
        (z1[0] - z2[0]) / (1 + k) + z2[0],
        (z1[1] - z2[1]) / (1 + k) + z2[1],
        (z1[2] - z2[2]) / (1 + k) + z2[2]};
    pstNbj->B = xyz2llr(FF);
}

void EphRspl::p2p(long double L, long double fa, _GJW &re, bool fAB, int f)
{
    // f取+-1
    _ZB &p = pstNbj->P, &q = pstNbj->Q;
    zbXY(pstNbj->P, L, fa);
    zbXY(pstNbj->Q, L, fa);

    long double u = q.y - p.y, v = q.x - p.x, a = sqrtl(u * u + v * v), r = 959.63 / p.S[2] / rad * cs_AU;
    long double W = p.S[1] + f * r * v / a, J = p.S[0] - f * r * u / a / cosl((W + p.S[1]) / 2.0), R = p.S[2];
    Vector3 AA = fAB ? pstNbj->A : pstNbj->B;

    COORDP pp = lineEar({J, W, R}, AA, p.g);
    re.J = pp.J;
    re.W = pp.W;
}

void EphRspl::pp0(_GJW &re)
{
    // 食中心点计算
    _ZB p = pstNbj->P;
    COORDP pp = lineEar(p.M, p.S, p.g);
    re.J = pp.J;
    re.W = pp.W; // 无解返回值是100

    if (re.W == 100)
    {
        re.c = "";
        return;
    }
    re.c = "全";
    zbXY(p, re.J, re.W);
    if (p.sr > p.mr)
    {
        re.c = "环";
    }
}

void EphRspl::secXY(long double jd, long double L, long double fa, long double high, _SECXY &re)
{
    // 日月xy坐标计算。参数：jd是力学时,站点经纬L,fa,海拔high(千米)
    // 基本参数计算
    auto deltat = dt_T(jd); // TD-UT
    Vector2 zd = nutation2(jd / 36525.0);
    auto gst = pGST(jd - deltat, deltat) + zd[0] * cosl(hcjj(jd / 36525.0) + zd[1]); // 真恒星时(不考虑非多项式部分)

    //=======月亮========
    EphRsgs &rsgsObj = EphRsgs::getInstance();
    auto z = rsgsObj.moon(jd);
    re.mCJ = z[0];
    re.mCW = z[1];
    re.mR = z[2];                          // 月亮视赤经,月球赤纬
    auto mShiJ = rad2rrad(gst + L - z[0]); // 得到此刻月亮时角
    z = parallax(z, mShiJ, fa, high);
    re.mCJ2 = z[0], re.mCW2 = z[1], re.mR2 = z[2]; // 修正了视差的赤道坐标

    //=======太阳========
    z = rsgsObj.sun(jd);
    re.sCJ = z[0];
    re.sCW = z[1];
    re.sR = z[2];                          // 太阳视赤经,太阳赤纬
    auto sShiJ = rad2rrad(gst + L - z[0]); // 得到此刻太阳时角
    z = parallax(z, sShiJ, fa, high);
    re.sCJ2 = z[0], re.sCW2 = z[1], re.sR2 = z[2]; // 修正了视差的赤道坐标

    //=======视半径========
    re.mr = cs_sMoon / re.mR2 / rad;
    re.sr = 959.63 / re.sR2 / rad * cs_AU;
    if (pstSeMax->nasa_r)
    {
        re.mr *= cs_sMoon2 / cs_sMoon; // 0.99925;
    }
    //=======日月赤经纬差转为日面中心直角坐标(用于日食)==============
    re.x = rad2rrad(re.mCJ2 - re.sCJ2) * cosl((re.mCW2 + re.sCW2) / 2.0);
    re.y = re.mCW2 - re.sCW2;
    re.t = jd;
}

SE_MAX_DATA &EphRspl::getMaxData()
{
    return *(this->pstSeMax);
}

SE_NBJ_DATA &EphRspl::getNbjData()
{
    return *(this->pstNbj);
}
