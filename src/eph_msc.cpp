#include <cstring>
#include <cmath>
#include "eph_msc.h"
#include "JD.h"

EphMSC::EphMSC()
{
    pstMscData = new S_EPH_MSC_DATA();
    std::memset(pstMscData, 0, sizeof(S_EPH_MSC_DATA));
}

EphMSC::~EphMSC()
{
    delete pstMscData;
}

void EphMSC::calc(long double T, long double L,long double fa,long double high)
{
    //sun_moon类的成员函数。参数：T是力学时,站点经纬L,fa,海拔high(千米)
    //基本参数计算
    pstMscData->T = T, pstMscData->L = L, pstMscData->fa = fa;
    pstMscData->dt = dt_T(T);			//TD-UT
    pstMscData->jd = T - pstMscData->dt;	//UT
    T /= 36525.0;
    Vector2 zd = nutation2(T);
    pstMscData->dL = zd[0];			//黄经章
    pstMscData->dE = zd[1];			//交角章动
    pstMscData->E = hcjj(T) + pstMscData->dE;	//真黄赤交角
    pstMscData->gst = pGST(pstMscData->jd, pstMscData->dt) + pstMscData->dL * cos(pstMscData->E);	//真恒星时(不考虑非多项式部分)
    Vector3 z;

    //=======月亮========
    //月亮黄道坐标
    z = m_coord(T, -1, -1, -1);	//月球坐标
    z[0] = rad2mrad(z[0] + gxc_moonLon(T) + pstMscData->dL);
    z[1] += gxc_moonLat(T);		//补上月球光行差及章动
    pstMscData->mHJ = z[0];
    pstMscData->mHW = z[1];
    pstMscData->mR = z[2];			//月球视黄经,视黄纬,地月质心距

    //月球赤道坐标
    z = llrConv(z, pstMscData->E);	//转为赤道坐标
    pstMscData->mCJ = z[0];
    pstMscData->mCW = z[1];			//月球视赤经,月球赤纬

    //月亮时角计算
    pstMscData->mShiJ = rad2mrad(pstMscData->gst + L - z[0]);	//得到此刻天体时角
    if (pstMscData->mShiJ > PI)
    {
        pstMscData->mShiJ -= pi2;
    }

    //修正了视差的赤道坐标
    z = parallax(z, pstMscData->mShiJ, fa, high);	//视差修正
    pstMscData->mCJ2 = z[0], pstMscData->mCW2 = z[1], pstMscData->mR2 = z[2];

    //月亮时角坐标
    z[0] += PI / 2 - pstMscData->gst - L;	//转到相对于地平赤道分点的赤道坐标(时角坐标)

    //月亮地平坐标
    z = llrConv(z, PI / 2 - fa);	//转到地平坐标(只改经纬度)
    z[0] = rad2mrad(PI / 2 - z[0]);
    pstMscData->mDJ = z[0];
    pstMscData->mDW = z[1];			//方位角,高度角
    if (z[1] > 0)
    {
        z[1] += MQC(z[1]);		//大气折射修正
    }
    pstMscData->mPJ = z[0];
    pstMscData->mPW = z[1];			//方位角,高度角

    //=======太阳========
    //太阳黄道坐标
    z = e_coord(T, -1, -1, -1);	//地球坐标
    z[0] = rad2mrad(z[0] + PI + gxc_sunLon(T) + pstMscData->dL);	//补上太阳光行差及章动
    z[1] = -z[1] + gxc_sunLat(T);	//z数组为太阳地心黄道视坐标
    pstMscData->sHJ = z[0];
    pstMscData->sHW = z[1];
    pstMscData->sR = z[2];			//太阳视黄经,视黄纬,日地质心距

    //太阳赤道坐标
    z = llrConv(z, pstMscData->E);	//转为赤道坐标
    pstMscData->sCJ = z[0];
    pstMscData->sCW = z[1];			//太阳视赤经,视赤纬

    //太阳时角计算
    pstMscData->sShiJ = rad2mrad(pstMscData->gst + L - z[0]);	//得到此刻天体时角
    if (pstMscData->sShiJ > PI)
    {
        pstMscData->sShiJ -= pi2;
    }

    //修正了视差的赤道坐标
    z = parallax(z, pstMscData->sShiJ, fa, high);	//视差修正
    pstMscData->sCJ2 = z[0], pstMscData->sCW2 = z[1], pstMscData->sR2 = z[2];

    //太阳时角坐标
    z[0] += PI / 2 - pstMscData->gst - L;	//转到相对于地平赤道分点的赤道坐标

    //太阳地平坐标
    z = llrConv(z, PI / 2 - fa);
    z[0] = rad2mrad(PI / 2 - z[0]);
    //z[1] -= 8.794/rad/z[2]*cos(z[1]); //直接在地平坐标中视差修正(这里把地球看为球形,精度比parallax()稍差一些)
    pstMscData->sDJ = z[0];
    pstMscData->sDW = z[1];			//方位角,高度角

    if (z[1] > 0)
    {
        z[1] += MQC(z[1]);		//大气折射修正
    }
    pstMscData->sPJ = z[0];
    pstMscData->sPW = z[1];			//方位角,高度角

    //=======其它========
    //时差计算
    long double t = T / 10, t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t;
    long double Lon = (1753470142 + 6283319653318 * t + 529674 * t2 + 432 * t3 - 1124 * t4 - 9 * t5) / 1000000000 + PI - 20.5 / rad;	//修正了光行差的太阳平黄经
    Lon = rad2mrad(Lon - (pstMscData->sCJ - pstMscData->dL * cos(pstMscData->E)));	//(修正了光行差的平黄经)-(不含dL*cos(E)的视赤经)
    if (Lon > PI)
    {
        Lon -= pi2;				//得到时差,单位是弧度
    }
    pstMscData->sc = Lon / pi2;		//时差(单位:日)

    //真太阳与平太阳
    pstMscData->pty = pstMscData->jd + L / pi2;	//平太阳时
    pstMscData->zty = pstMscData->jd + L / pi2 + pstMscData->sc;	//真太阳时

    //视半径
    // pstMscData->mRad = moonRad(pstMscData->mR,pstMscData->mDW);  //月亮视半径(角秒)
    pstMscData->mRad = cs_sMoon / pstMscData->mR2;	//月亮视半径(角秒)
    pstMscData->sRad = 959.63 / pstMscData->sR2;	//太阳视半径(角秒)
    pstMscData->e_mRad = cs_sMoon / pstMscData->mR;	//月亮地心视半径(角秒)
    pstMscData->eShadow = (cs_rEarA / pstMscData->mR * rad - (959.63 - 8.794) / pstMscData->sR) * 51 / 50;	//地本影在月球向径处的半径(角秒),式中51/50是大气厚度补偿
    pstMscData->eShadow2 = (cs_rEarA / pstMscData->mR * rad + (959.63 + 8.794) / pstMscData->sR) * 51 / 50;	//地半影在月球向径处的半径(角秒),式中51/50是大气厚度补偿
    pstMscData->mIll = moonIll(T);	//月亮被照面比例

    //中心食计算
    if (std::abs(rad2rrad(pstMscData->mCJ - pstMscData->sCJ)) < 50.0 / 180.0 * PI)
    {
        COORDP pp = lineEar({ pstMscData->mCJ, pstMscData->mCW, pstMscData->mR }
                , { pstMscData->sCJ, pstMscData->sCW, pstMscData->sR * cs_AU }
                , pstMscData->gst);
        pstMscData->zx_J = pp.J;
        pstMscData->zx_W = pp.W;		//无解返回值是100
    }
    else

    {
        pstMscData->zx_J = pstMscData->zx_W = 100;
    }
}

std::string EphMSC::toString(bool flag)
{
    std::string s;
    s = "-------------------------------------------\n";
    s = s + "平太阳 " + JD::timeStr(pstMscData->pty) + " 真太阳 " + JD::timeStr(pstMscData->zty) + "\n";
    s = s + "时差 " + m2fm(pstMscData->sc * 86400, 2, 1) + " 月亮被照亮 " + to_str(pstMscData->mIll * 100, 2) + "% ";
    s = s + "\n";

    s = s + "-------------------------------------------\n表一       月亮            太阳\n";
    s = s + "视黄经 " + rad2str(pstMscData->mHJ, 0) + "  " + rad2str(pstMscData->sHJ, 0) + "\n";
    s = s + "视黄纬 " + rad2str(pstMscData->mHW, 0) + "  " + rad2str(pstMscData->sHW, 0) + "\n";
    s = s + "视赤经 " + rad2str(pstMscData->mCJ, 1) + "  " + rad2str(pstMscData->sCJ, 1) + "\n";
    s = s + "视赤纬 " + rad2str(pstMscData->mCW, 0) + "  " + rad2str(pstMscData->sCW, 0) + "\n";
    s = s + "距离    " + to_str(pstMscData->mR, 2) + "千米    " + to_str(pstMscData->sR, 8) + "AU" + "\n";

    s = s + "-------------------------------------------\n表二       月亮            太阳\n";
    s = s + "方位角 " + rad2str(pstMscData->mPJ, 0) + "  " + rad2str(pstMscData->sPJ, 0) + "\n";
    s = s + "高度角 " + rad2str(pstMscData->mPW, 0) + "  " + rad2str(pstMscData->sPW, 0) + "\n";
    s = s + "时角   " + rad2str(pstMscData->mShiJ, 0) + "  " + rad2str(pstMscData->sShiJ, 0) + "\n";
    s = s + "视半径   " + m2fm(pstMscData->mRad, 2, 0) + "       " + m2fm(pstMscData->sRad, 2, 0) + " (观测点)\n";

    if (flag)
    {
        s = s + "-------------------------------------------\n";
        s = s + "力学时" + JD::JD2str(pstMscData->T + J2000);
        s = s + " ΔT=" + to_str(pstMscData->dt * 86400, 1) + "秒\n";
        s = s + "黄经章 " + to_str(pstMscData->dL / pi2 * 360 * 3600, 2) + "\" ";
        s = s + "交角章 " + to_str(pstMscData->dE / pi2 * 360 * 3600, 2) + "\" ";
        s = s + "ε=" + rad2str(pstMscData->E, 0) + "\n";
        s = s + "-------------------------------------------\n";

    }
    return s;
}

S_EPH_MSC_DATA& EphMSC::getData()
{
    return *(this->pstMscData);
}
