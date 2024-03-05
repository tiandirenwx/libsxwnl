#include "eph_yspl.h"
#include "eph.h"

EphYspl::EphYspl()
{
}
EphYspl::~EphYspl()
{
}
long double EphYspl::lineT(RE0 G, long double v, long double u, long double r, bool n)
{ // 已知t1时刻星体位置、速度，求x*x+y*y=r*r时,t的值
    long double b = G.y * v - G.x * u, A = u * u + v * v, B = u * b, C = b * b - r * r * v * v, D = B * B - A * C;
    if (D < 0)
    {
        return 0;
    }
    D = sqrt(D);
    if (!n)
    {
        D = -D;
    }
    return G.t + ((-B + D) / A - G.x) / v;
}

void EphYspl::lecXY(long double jd, RE0 &re)
{
    // 日月黄经纬差转为日面中心直角坐标(用于月食)
    long double T = jd / 36525.0;
    Vector3 zm = {}, zs = {};

    //=======太阳月亮黄道坐标========
    zs = e_coord(T, -1, -1, -1); // 地球坐标
    zs[0] = rad2mrad(zs[0] + PI + gxc_sunLon(T));
    zs[1] = -zs[1] + gxc_sunLat(T); // 补上太阳光行差
    zm = m_coord(T, -1, -1, -1);    // 月球坐标
    zm[0] = rad2mrad(zm[0] + gxc_moonLon(T));
    zm[1] += gxc_moonLat(T); // 补上月球光行差就可以了

    //=======视半径=======
    re.e_mRad = cs_sMoon / zm[2];                                                // 月亮地心视半径(角秒)
    re.eShadow = (cs_rEarA / zm[2] * rad - (959.63 - 8.794) / zs[2]) * 51 / 50.0;  // 地本影在月球向径处的半径(角秒),式中51/50是大气厚度补偿
    re.eShadow2 = (cs_rEarA / zm[2] * rad + (959.63 + 8.794) / zs[2]) * 51 / 50.0; // 地半影在月球向径处的半径(角秒),式中51/50是大气厚度补偿

    re.x = rad2rrad(zm[0] + PI - zs[0]) * cos((zm[1] - zs[1]) / 2.0);
    re.y = zm[1] + zs[1];

    re.mr = re.e_mRad / rad,
    re.er = re.eShadow / rad,
    re.Er = re.eShadow2 / rad;

    re.t = jd;
}

void EphYspl::lecMax(long double jd)
{ // 月食的食甚计算(jd为近朔的力学时,误差几天不要紧)
    // EphYspl::lT={};
    for (int i = 0; i < 7; i++)
    {
        lT[i] = 0; // 分别是:食甚,初亏,复圆,半影食始,半影食终,食既,生光
    }
    sf = 0;
    LX = "";

    jd = XL::MS_aLon_t2(floor((jd - 4) / 29.5306) * PI * 2 + PI) * 36525; // 低精度的朔(误差10分钟),与食甚相差10分钟左右

    RE0 g = {}, G = {};
    long double u, v;

    // 求极值(平均误差数秒)
    u = -18461 * sin(0.057109 + 0.23089571958 * jd) * 0.23090 / rad; // 月日黄纬速度差
    v = (XL::M_v(jd / 36525.0) - XL::E_v(jd / 36525.0)) / 36525.0;   // 月日黄经速度差
    lecXY(jd, G);
    jd -= (G.y * u + G.x * v) / (u * u + v * v); // 极值时间

    // 精密求极值
    long double dt = 60 / 86400.0;
    lecXY(jd, G);
    lecXY(jd + dt, g); // 精密差分得速度,再求食甚
    u = (g.y - G.y) / dt;
    v = (g.x - G.x) / dt;
    dt = -(G.y * u + G.x * v) / (u * u + v * v);
    jd += dt; // 极值时间

    // 求直线到影子中心的最小值
    long double x = G.x + dt * v, y = G.y + dt * u, rmin = sqrt(x * x + y * y);
    // 注意,以上计算得到了极值及最小距rmin,但没有再次计算极值时刻的半径,对以下的判断造成一定的风险,必要的话可以再算一次。不过必要性不很大，因为第一次极值计算已经很准确了,误差只有几秒
    // 求月球与影子的位置关系
    if (rmin <= G.mr + G.er)
    {               // 食计算
        lT[1] = jd; // 食甚
        LX = "偏";
        sf = (G.mr + G.er - rmin) / G.mr / 2; // 食分

        lT[0] = lineT(G, v, u, G.mr + G.er, 0); // 初亏
        lecXY(lT[0], g);
        lT[0] = lineT(g, v, u, g.mr + g.er, 0); // 初亏再算一次
        //	std::cout<<EphYspl::lT[0]<<std::endl;
        lT[2] = lineT(G, v, u, G.mr + G.er, 1); // 复圆
        lecXY(lT[2], g);
        lT[2] = lineT(g, v, u, g.mr + g.er, 1); // 复圆再算一次
    }
    if (rmin <= G.mr + G.Er)
    {                                           // 半影食计算
        lT[3] = lineT(G, v, u, G.mr + G.Er, 0); // 半影食始
        lecXY(lT[3], g);
        lT[3] = lineT(g, v, u, g.mr + g.Er, 0); // 半影食始再算一次

        lT[4] = lineT(G, v, u, G.mr + G.Er, 1); // 半影食终
        lecXY(lT[4], g);
        lT[4] = lineT(g, v, u, g.mr + g.Er, 1); // 半影食终再算一次
    }
    if (rmin <= G.er - G.mr)
    { // 全食计算
        LX = "全";
        lT[5] = lineT(G, v, u, G.er - G.mr, 0); // 食既
        lecXY(lT[5], g);
        lT[5] = lineT(g, v, u, g.er - g.mr, 0); // 食既再算一次

        lT[6] = lineT(G, v, u, G.er - G.mr, 1); // 生光
        lecXY(lT[6], g);
        lT[6] = lineT(g, v, u, g.er - g.mr, 1); // 生光再算一次
    }
}

std::string EphYspl::getLX() const
{
    return LX;
}

long double EphYspl::getSF() const
{
    return sf;
}

std::array<long double, 7> EphYspl::getLT() const
{
    return lT;
}
