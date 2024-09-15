#include "eph_rsgs.h"
#include "JD.h"

EphRsgs::EphRsgs()
{
    pstPara = new RSGS_PARA();
    //pstPara->Zjd = 0;
    pstPara->setZjd(0);
    //pstPara->Zdt = 0.04;     // 插值点之间的时间间距
    pstPara->setZdt(0.04);
    //pstPara->dT = 0;         // deltatT
    pstPara->setdT(0);
    //pstPara->tanf1 = 0.0046; // 半影锥角
    pstPara->setTanf1(0.0046);
    //pstPara->tanf2 = 0.0045; // 本影锥角
    pstPara->setTanf2(0.0045);
    //pstPara->srad = 0.0046;  // 太阳视半径
    pstPara->setSrad(0.0046);
    //pstPara->bba = 1;        // 贝圆极赤比
    pstPara->setBba(1);
    //pstPara->bhc = 0;        // 黄交线与赤交线的夹角简易作图用
    pstPara->setBhc(0);
    //pstPara->dyj = 23500;    // 地月距
    pstPara->setDyj(23500);
}

EphRsgs::~EphRsgs()
{
    delete pstPara;
}

void EphRsgs::init(long double jd, int n)
{
    // 创建插值表(根数表)
    if (suoN(jd) == suoN(pstPara->getZjd()) && Zs.size() == n * 9)
    {
        return;
    }
    Zs.clear();
    jd = XL::MS_aLon_t2(suoN(jd) * PI * 2) * 36525; // 低精度的朔(误差10分钟)
    pstPara->setZjd(jd);
    pstPara->setdT(dt_T(jd));     // deltat T

    Vector2 zd = nutation2(jd / 36525.0);       // 章动
    long double E = hcjj(jd / 36525.0) + zd[1]; // 黄赤交角
    long double T;
    Vector3 S, M, B;
    std::vector<long double> a(n * 9);
    int i, k;

    for (i = 0; i < n; i++)
    { // 插值点范围不要超过360度(约1个月)
        T = (pstPara->getZjd() + (i - n / 2.0 + 0.5) * pstPara->getZdt()) / 36525.0;

        if (n == 7)
        {
            S = e_coord(T, -1, -1, -1), M = m_coord(T, -1, -1, -1); // 地球坐标及月球坐标,全精度
        }
        if (n == 3)
        {
            S = e_coord(T, 65, 65, 65), M = m_coord(T, -1, 150, 150); // 中精度
        }
        if (n == 2)
        {
            S = e_coord(T, 20, 20, 20), M = m_coord(T, 30, 30, 30); // 低精度
        }

        S[0] = S[0] + zd[0] + gxc_sunLon(T) + PI;
        S[1] = -S[1] + gxc_sunLat(T); // 补上太阳光行差及章动
        M[0] = M[0] + zd[0] + gxc_moonLon(T);
        M[1] = M[1] + gxc_moonLat(T); // 补上月球光行差及章动
        S = llrConv(S, E);
        M = llrConv(M, E);
        S[2] *= cs_AU; // 转为赤道坐标
        if (i && S[0] < a[0])
        {
            S[0] += pi2; // 确保插值数据连续
        }
        if (i && M[0] < a[3])
        {
            M[0] += pi2; // 确保插值数据连续
        }

        k = i * 9;
        a[k + 0] = S[0], a[k + 1] = S[1], a[k + 2] = S[2]; // 存入插值表
        a[k + 3] = M[0], a[k + 4] = M[1], a[k + 5] = M[2];

        // 贝塞尔坐标的z轴坐标计算,得到a[k+6,7,8]交点赤经,贝赤交角,真恒星时
        S = llr2xyz(S), M = llr2xyz(M);
        B = xyz2llr({S[0] - M[0], S[1] - M[1], S[2] - M[2]});
        B[0] = PI / 2 + B[0];
        B[1] = PI / 2 - B[1];
        if (i && B[0] < a[6])
        {
            B[0] += pi2; // 确保插值数据连续
        }

        a[k + 6] = B[0], a[k + 7] = B[1], a[k + 8] = pGST(T * 36525 - pstPara->getdT(), pstPara->getdT()) + zd[0] * cos(E); // 真恒星时
    }
    Zs = a;
    // 一些辅助参数的计算
    auto p = a.size() - 9;
    pstPara->setDyj((a[2] + a[p + 2] - a[5] - a[p + 5]) / 2.0 / cs_rEar); // 地月平均距离
    pstPara->setTanf1((cs_k0 + cs_k) / pstPara->getDyj());                     // tanf1半影锥角
    pstPara->setTanf2((cs_k0 - cs_k2) / pstPara->getDyj());                    // tanf2本影锥角
    pstPara->setSrad(cs_k0 / ((a[2] + a[p + 2]) / 2.0 / cs_rEar));
    pstPara->setBba(sinl((a[1] + a[p + 1]) / 2.0));
    pstPara->setBba(cs_ba * (1 + (1 - cs_ba2) * pstPara->getBba() * pstPara->getBba() / 2.0));
    pstPara->setBhc(-atan(tanl(E) * sinl((a[6] + a[p + 6]) / 2.0))); // 黄交线与赤交线的夹角
}

Vector3 EphRsgs::chazhi(long double jd, int xt)
{
    // 日月坐标快速计算(贝赛尔插值法),计算第p个根数开始的m个根数
    int p = xt * 3, m = 3; // 计算第p个根数开始的m个根数
    int i, N = Zs.size() / 9;
    std::vector<long double> &B = Zs;
    Vector3 z;
    int w = B.size() / N;                                               // 每节点个数
    long double t = (jd - pstPara->getZjd()) / pstPara->getZdt() + N / 2.0 - 0.5; // 相对于第一点的时间距离

    if (N == 2)
    {
        for (i = 0; i < m; i++, p++)
        {
            z[i] = B[p] + (B[p + w] - B[p]) * t;
        }
        return z;
    }
    auto c = int64(t + 0.5);
    if (c <= 0)
    {
        c = 1;
    }
    if (c > N - 2)
    {
        c = N - 2; // 确定c,并对超出范围的处理
    }
    t -= c, p += c * w; // c插值中心,t为插值因子,t再转为插值中心在数据中的位置
    for (i = 0; i < m; i++, p++)
    {
        z[i] = B[p] + (B[p + w] - B[p - w] + (B[p + w] + B[p - w] - B[p] * 2) * t) * t / 2;
    }
    return z;
}

Vector3 EphRsgs::cd2bse(Vector3 z, Vector3 I)
{ // 赤道转贝塞尔坐标
    Vector3 r{z[0] - I[0], z[1], z[2]};
    r = llrConv(r, -I[1]);
    return llr2xyz(r);
}
Vector3 EphRsgs::bse2cd(Vector3 z, Vector3 I)
{ // 贝塞尔转赤道坐标
    Vector3 r = xyz2llr(z);
    r = llrConv(r, I[1]);
    r[0] = rad2mrad(r[0] + I[0]);
    return r;
}
Vector3 EphRsgs::bse2db(Vector3 z, Vector3 I, bool f)
{ // 贝赛尔转地标(p点到原点连线与地球的交点,z为p点直角坐标),f=1时把地球看成椭球
    Vector3 r = xyz2llr(z);
    r = llrConv(r, I[1]);
    r[0] = rad2rrad(r[0] + I[0] - I[2]);
    if (f)
    {
        r[1] = atanl(tanl(r[1]) / cs_ba2);
    }
    return r;
}
Vector3 EphRsgs::bseXY2db(long double x, long double y, Vector3 I, bool f)
{
    // 贝赛尔转地标(过p点垂直于基面的线与地球的交点,p坐标为(x,y,任意z)),f=1时把地球看成椭球
    long double b = f ? cs_ba : 1;
    COORDP F = lineEar2(x, y, 2, x, y, 0, b, 1, I); // 求中心对应的地标
    return {F.J, F.W, 0};
}

Vector3 EphRsgs::bseM(long double jd)
{ // 月亮的贝塞尔坐标
    Vector3 a = cd2bse(chazhi(jd, 1), chazhi(jd, 2));
    a[0] /= cs_rEar, a[1] /= cs_rEar, a[2] /= cs_rEar;
    return a;
}

// 以下计算日食总体情况

_VXY EphRsgs::Vxy(long double x, long double y, long double s, long double vx, long double vy)
{ // 地球上一点的速度，用贝塞尔坐标表达，s为贝赤交角
    _VXY r = {};
    long double h = 1 - x * x - y * y;
    if (h < 0)
    {
        h = 0; // 越界置0，使速度场连续，置零有助于迭代时单向收敛
    }
    else
    {
        h = sqrtl(h);
    }
    r.vx = pi2 * (sinl(s) * h - cosl(s) * y);
    r.vy = pi2 * x * cosl(s);
    r.Vx = vx - r.vx;
    r.Vy = vy - r.vy;
    r.V = sqrtl(r.Vx * r.Vx + r.Vy * r.Vy);
    return r;
}
_RSM EphRsgs::rSM(long double mR)
{ // rm,rs单位千米
    _RSM re = {};
    re.r1 = cs_k + pstPara->getTanf1() * mR;  // 半影半径
    re.r2 = cs_k2 - pstPara->getTanf2() * mR; // 本影半径
    re.ar2 = std::abs(re.r2);
    re.sf = cs_k2 / mR / cs_k0 * (pstPara->getDyj() + mR); // 食分
    return re;
}
Vector3 EphRsgs::qrd(long double jd, long double dx, long double dy, bool fs)
{ // 求切入点
    long double ba2 = pstPara->getBba() * pstPara->getBba();
    Vector3 M = bseM(jd);
    long double x = M[0], y = M[1];
    _RSM B = rSM(M[2]);
    long double r = 0;
    if (fs)
    {
        r = B.r1;
    }
    auto d = 1 - (1.0 / ba2 - 1) * y * y / (x * x + y * y) / 2.0 + r;
    auto t = (d * d - x * x - y * y) / (dx * x + dy * y) / 2;
    x += t * dx, y += t * dy, jd += t;

    auto c = (1 - ba2) * r * x * y / d / d / d;
    x += c * y;
    y -= c * x;
    Vector3 re = bse2db({x / d, y / d, 0}, bse(jd), 1);
    // re[0] +=0.275/radd; //转为deltatT为66秒的历书经度
    re[2] = jd;
    return re;
}

_FEATURE EphRsgs::feature(long double jd)
{
    // 日食的基本特征
    jd = pstPara->getZjd(); // 低精度的朔(误差10分钟)

    long double tg = 0.04;
    long double jd1 = jd - tg / 2.0;
    _FEATURE re = {};

    long double ls;
    Vector3 a = bseM(jd - tg);
    Vector3 b = bseM(jd);
    Vector3 c = bseM(jd + tg);
    auto vx = (c[0] - a[0]) / tg / 2.0;
    auto vy = (c[1] - a[1]) / tg / 2.0;
    auto vz = (c[2] - a[2]) / tg / 2.0;
    auto ax = (c[0] + a[0] - 2 * b[0]) / tg / tg;
    auto ay = (c[1] + a[1] - 2 * b[1]) / tg / tg;
    auto v = sqrtl(vx * vx + vy * vy), v2 = v * v;

    // 影轴在贝塞尔面扫线的特征参数
    re.jdSuo = jd;        // 朔
    re.dT = pstPara->getdT();  // deltat T
    re.ds = pstPara->getBhc(); // 黄交线与赤交线的夹角
    re.vx = vx;           // 影速x
    re.vy = vy;           // 影速y
    re.ax = ax;
    re.ay = ay;
    re.v = v;
    re.k = vy / vx; // 斜率

    auto t0 = -(b[0] * vx + b[1] * vy) / v2;
    re.jd = jd + t0;                         // 中点时间
    re.xc = b[0] + vx * t0;                  // 中点坐标x
    re.yc = b[1] + vy * t0;                  // 中点坐标y
    re.zc = b[2] + vz * t0 - 1.37 * t0 * t0; // 中点坐标z
    re.D = (vx * b[1] - vy * b[0]) / v;
    re.d = fabsl(re.D); // 直线到圆心的距离
    re.I = bse(re.jd); // 中心点的贝塞尔z轴的赤道坐标及恒星时，(J,W,g)

    // 影轴交点判断
    COORDP F = lineEar2(re.xc, re.yc, 2, re.xc, re.yc, 0, cs_ba, 1, re.I); // 求中心对应的地标
    // 四个关键点的影子半径计算
    _RSM Bc, Bp, B2, B3;
    long double dt, t2, t3, t4, t5, t6;
    Bc = Bp = B2 = B3 = rSM(re.zc); // 中点处的影子半径
    if (F.W != 100)
    {
        Bp = rSM(re.zc - F.R2);
    }
    if (re.d < 1)
    {
        dt = sqrtl(1 - re.d * re.d) / v;
        t2 = t0 - dt, t3 = t0 + dt;                // 中心始终参数
        B2 = rSM(t2 * vz + b[2] - 1.37 * t2 * t2); // 中心线始影半径
        B3 = rSM(t3 * vz + b[2] - 1.37 * t3 * t3); // 中心线终影半径
    }
    ls = 1;
    dt = 0;
    if (re.d < ls)
    {
        dt = sqrtl(ls * ls - re.d * re.d) / v;
    }
    t2 = t0 - dt, t3 = t0 + dt; // 偏食始终参数,t2,t3
    ls = 1 + Bc.r1;
    dt = 0;
    if (re.d < ls)
    {
        dt = sqrtl(ls * ls - re.d * re.d) / v;
    }
    t4 = t0 - dt, t5 = t0 + dt; // 偏食始终参数,t4,t5
    t6 = -b[0] / vx;            // 视午参数l6

    if (re.d < 1)
    {
        re.gk1 = qrd(t2 + jd, vx, vy, 0); // 中心始
        re.gk2 = qrd(t3 + jd, vx, vy, 0); // 中心终
    }
    else
    {
        re.gk1 = {0, 0, 0};
        re.gk2 = {0, 0, 0};
    }

    re.gk3 = qrd(t4 + jd, vx, vy, 1); // 偏食始
    re.gk4 = qrd(t5 + jd, vx, vy, 1); // 偏食终
    re.gk5 = bseXY2db(t6 * vx + b[0], t6 * vy + b[1], bse(t6 + jd), 1);
    re.gk5[2] = t6 + jd; // 地方视午日食

    // 日食类型、最大食地标、食分、太阳地平坐标
    Vector3 lls;
    if (F.W == 100)
    { // 无中心线
        // 最大食地标及时分
        lls = bse2db({re.xc, re.yc, 0}, re.I, 0);
        re.zxJ = lls[0], re.zxW = lls[1];                    // 最大食地标
        re.sf = (Bc.r1 - (re.d - 0.9972)) / (Bc.r1 - Bc.r2); // 0.9969是南北极区的平半径
        // 类型判断
        if (re.d > 0.9972 + Bc.r1)
        {
            re.lx = "N";
        } // 无食,半影没有进入
        else if (re.d > 0.9972 + Bc.ar2)
        {
            re.lx = "P";
        } // 偏食,本影没有进入
        else
        {
            if (Bc.sf < 1)
                re.lx = "A0";
            else
                re.lx = "T0";
        } // 中心线未进入,本影部分进入(无中心，所以只是部分地入)
    }
    else
    { // 有中心线
        // 最大食地标及时分
        re.zxJ = F.J, re.zxW = F.W; // 最大食地标
        re.sf = Bp.sf;              // 
        // 类型判断
        if (re.d > 0.9966 - Bp.ar2)
        {
            if (Bp.sf < 1)
            {
                re.lx = "A1";
            }
            else
            {
                re.lx = "T1";
            }
        } // 中心进入,但本影没有完全进入
        else
        { // 本影全进入有中心日食
            if (Bp.sf >= 1)
            {
                re.lx = "H";
                if (B2.sf > 1)
                {
                    re.lx = "H2"; // 全环食,全始
                }
                if (B3.sf > 1)
                {
                    re.lx = "H3"; // 全环食,全终
                }
                if (B2.sf > 1 && B3.sf > 1)
                {
                    re.lx = "T"; // 全食
                }
            }
            else
            {
                re.lx = "A"; // 环食
            }
        }
    }
    re.Sdp = CD2DP(sun(re.jd), re.zxJ, re.zxW, re.I[2]); // 太阳在中心点的地平坐标

    // 食带宽度和时延
    if (F.W != 100)
    {
        re.dw = std::abs(2 * Bp.r2 * cs_rEar) / sinl(re.Sdp[1]);   // 食带宽度
        _VXY llls = Vxy(re.xc, re.yc, re.I[1], re.vx, re.vy); // 求地表影速
        re.tt = 2 * std::abs(Bp.r2) / llls.V;                     // 时延
    }
    else
    {
        re.dw = re.tt = 0;
    }
    return re;
}

// 界线图
void EphRsgs::push(Vector3 z, std::vector<long double> &p)
{
    p.push_back(z[0]); // 保存第一食甚线A或B根
    p.push_back(z[1]);
}

Vector4 EphRsgs::nanbei(Vector3 M, long double vx0, long double vy0, long double h, long double r, Vector3 I)
{
    // vx0,vy0为影足速度(也是整个影子速度),h=1计算北界,h=-1计算南界
    long double x = M[0] - vy0 / vx0 * r * h, y = M[1] + h * r, z;
    long double vx, vy, v, sinA, cosA;

    for (int i = 0, js = 0; i < 3; i++)
    {
        z = 1 - x * x - y * y;
        if (z < 0)
        {
            if (js)
            {
                break;
            }
            z = 0;
            js++;
        } // z小于0则置0，如果两次小于0，可能不收敛造成的，故不再迭代了
        z = sqrt(z);
        x -= (x - M[0]) * z / M[2];
        y -= (y - M[1]) * z / M[2];
        vx = vx0 - pi2 * (sinl(I[1]) * z - cosl(I[1]) * y);
        vy = vy0 - pi2 * cosl(I[1]) * x;
        v = sqrtl(vx * vx + vy * vy);
        sinA = h * vy / v, cosA = h * vx / v;
        x = M[0] - r * sinA, y = M[1] + r * cosA;
    }
    long double X = M[0] - cs_k * sinA, Y = M[1] + cs_k * cosA;
    COORDP p = lineEar2(X, Y, M[2], x, y, 0, cs_ba, 1, I);
    return {p.J, p.W, x, y};
}

bool EphRsgs::mDian(Vector3 M, long double vx0, long double vy0, bool AB, long double r, Vector3 I, std::vector<long double> &A)
{ // 日出日没食甚
    long double R;
    NODE p;
    Vector3 a = M;
    _VXY c = {};
    for (int i = 0; i < 2; i++)
    { // 迭代求交点
        c = Vxy(a[0], a[1], I[1], vx0, vy0);
        p = lineOvl(M[0], M[1], c.Vy, -c.Vx, 1, pstPara->getBba());
        if (!p.n)
        {
            break;
        }
        if (AB)
        {
            a = {p.A[0], p.A[1], 0}, R = p.R1;
        }
        else
        {
            a = {p.B[0], p.B[1], 0}, R = p.R2;
        }
    }
    if (p.n && R <= r)
    {                                      // 有交点
        a = bse2db({a[0], a[1], 0}, I, 1); // 转为地标
        A.push_back(a[0]);                 // 保存第一食甚线A或B根
        A.push_back(a[1]);
        return true;
    }
    return false;
}

std::string EphRsgs::jieX3(long double jd)
{
    // 界线表
    long double k, ls;
    Vector4 p;
    int i;
    _FEATURE re = feature(jd); // 求特征参数

    long double t = int64(re.jd * 1440) / 1440.0 - 3 / 24.0;
    long double N = 360, dt = 1 / 1440.0;
    std::string s = "", s2;

    for (i = 0; i < N; i++, t += dt)
    {
        long double vx = re.vx + re.ax * (t - re.jdSuo);
        long double vy = re.vy + re.ay * (t - re.jdSuo);
        Vector3 M = bseM(t);  // 此刻月亮贝塞尔坐标(其x和y正是影足)
        _RSM B = rSM(M[2]);   // 本半影等
        long double r = B.r1; // 半影半径
        Vector3 I = bse(t);   // 贝塞尔坐标参数
        s2 = JD::JD2str(t + J2000) + " ", k = 0;
        // 南北界
        p = nanbei(M, vx, vy, +1, r, I);
        //if (p[1] != 100)
        if(notEqual(p[1],100))
        {
            s2 += rad2str2(p[0]) + "  " + rad2str2(p[1]) + " |", k++;
        }
        else
        {
            s2 += "-------------------|"; // 半影北界
        }
        p = nanbei(M, vx, vy, +1, B.r2, I);
        //if (p[1] != 100)
        if(notEqual(p[1],100))
        {
            s2 += rad2str2(p[0]) + "  " + rad2str2(p[1]) + " |", k++;
        }
        else
        {
            s2 += "-------------------|"; // 本影北界
        }
        Vector3 pp = bseXY2db(M[0], M[1], I, 1);
        p = {pp[0], pp[1], pp[2], 0};
        //if (p[1] != 100)
        if(notEqual(p[1],100))
        {
            s2 += rad2str2(p[0]) + "  " + rad2str2(p[1]) + " |", k++;
        }
        else
        {
            s2 += "-------------------|"; // 中心线
        }
        p = nanbei(M, vx, vy, -1, B.r2, I);
        //if (p[1] != 100)
        if(notEqual(p[1],100))
        {
            s2 += rad2str2(p[0]) + "  " + rad2str2(p[1]) + " |", k++;
        }
        else
        {
            s2 += "-------------------|"; // 本影南界
        }
        p = nanbei(M, vx, vy, -1, r, I);
        //if (p[1] != 100)
        if(notEqual(p[1],100))
        {
            s2 += rad2str2(p[0]) + " " + rad2str2(p[1]) + " ", k++;
        }
        else
        {
            s2 += "------------------- "; // 半影南界
        }
        if (k)
        {
            s += s2 + "\n";
        }
    }
    return "\033[41;37;1m 时间(力学时) 半影北界限 本影北界线 中心线 本影南界线 半影南界线\033[0m(伪本影南北界应互换)\n\n\n\n" + s;
}

RSGS_PARA &EphRsgs::getRsgsParaData()
{
    return *(this->pstPara);
}


