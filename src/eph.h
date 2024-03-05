#pragma once
#include <cstdint>
#include <cstring>
#include <string>
#include <array>
#include "const.h"
#include "eph_data.h"

typedef std::array<long double,2> Vector2;
typedef std::array<long double,3> Vector3;
typedef std::array<long double,4> Vector4;

struct COORDP
{
	long double x;
	long double y;
	long double z;
	long double R1;
	long double R2;
	long double D;
	long double X;
	long double J;
	long double W;
};

struct NODE
{
	Vector3 A;
	Vector3 B;
	long double R1;
	long double R2;
	int n;
};

//====================================
/*****
ecFast()函数返回参数说明
r.jdSuo 朔时刻
r.lx    日食类型
*****/
typedef struct
{
	long double jd;
	long double jdSuo;
	long double ac;
	std::string lx;
} _ECFAST;

//=================================数学工具=========================================

// 对超过0-2PI的角度转为0-2PI;
long double rad2mrad(long double v);
// 对超过-PI到PI的角度转为-PI到PI;
long double rad2rrad(long double v);

// 将弧度转为字串,精确到分
std::string rad2str2(long double d);
// 将弧度转为字串,保留2位
std::string rad2str(long double d, bool tim);
//===============角度格式化==================
std::string rad2strE(long double d, bool flag, int ext);

std::string toFixed(long double z, int n);

std::string fill_str(std::string s, int n, std::string c);
std::string dtoa_milo2(long double value, int precision, bool rightpad);

std::string to_str(long double value, int precision = 4, int length = 0, bool right_align = false);


// 秒转为分秒,fx为小数点位数,fs为1转为"分秒"格式否则转为"角度分秒"格式
std::string m2fm(long double v, int fx, int fs);

// 临界余数(a与最近的整倍数b相差的距离);
long double mod2(long double a, long double b);
// 球面转直角坐标;
Vector3 llr2xyz(Vector3 JW);
// 直角坐标转球;
Vector3 xyz2llr(Vector3 xyz);
// 球面坐标旋转;
Vector3 llrConv(Vector3 JW, long double E);
// 赤道坐标转为地平坐标;
Vector3 CD2DP(Vector3 z, long double L, long double fa, long double gst);
// 求角度差;
long double j1_j2(long double J1, long double W1, long double J2, long double W2);
// 日心球面转地心球面,Z星体球面坐标,A地球球面坐标;
// 本含数是通用的球面坐标中心平移函数,行星计算中将反复使用;
Vector3 h2g(Vector3 z, Vector3 a);
// 视差角(不是视差);
long double shiChaJ(long double gst, long double L, long double fa, long double J, long double W);

long double pGST2(long double jd);
long double pGST(long double T, long double dt);
long double gxc_moonLon(long double t);
long double gxc_moonLat(long double t);
int64_t suoN(long double jd);

Vector2 nutation2(long double t);
long double hcjj(long double t);
Vector3 e_coord(long double t, int n1, int n2, int n3);
Vector3 m_coord(long double t, int n1, int n2, int n3);
Vector3 p_coord(int xt, long double t, int n1, int n2, int n3);
long double gxc_sunLon(long double t);
long double gxc_sunLat(long double t);

long double sunShengJ(long double jd, long double L, long double fa, long double sj);

Vector3 parallax(Vector3 z, long double H, long double fa, long double high);

long double MQC(long double h);
long double MQC2(long double ho);

long double pty_zty(long double t);


// 星历计算
long double XL0_calc(int xt, int zn, long double t, int n);
long double XL1_calc(int zn, long double t, int n);


// 物件XL : 日月黄道平分点坐标、视坐标、速度、已知经度反求时间等方面的计算
namespace XL
{
	//=====================
	// 星历函数(日月球面坐标计算)

	long double E_Lon(long double t, int n); // 地球经度计算,返回Date分点黄经,传入世纪数、取项数
	long double M_Lon(long double t, int n); // 月球经度计算,返回Date分点黄经,传入世纪数,n是项数比例
	// 地球速度,t是世纪数,误差小于万分3													//=========================
	long double E_v(long double t);

	// 月球速度计算,传入世经数
	long double M_v(long double t);

	//=========================
	// 月日视黄经的差值
	long double MS_aLon(long double t, long double Mn, long double Sn);

	// 太阳视黄经
	long double S_aLon(long double t, long double n);

	//=========================
	// 已知地球真黄经求时间
	long double E_Lon_t(long double W);

	// 已知真月球黄经求时间
	long double M_Lon_t(long double W);

	// 已知月日视黄经差求时间
	long double MS_aLon_t(long double W);

	// 已知太阳视黄经反求时间
	long double S_aLon_t(long double W);

	// 已知月日视黄经差求时间,高速低精度,误差不超过600秒(只验算了几千年)
	long double MS_aLon_t2(long double W);
	// 已知太阳视黄经反求时间,高速低精度,最大误差不超过600秒
	long double S_aLon_t2(long double W);

	long double moonIll(long double t);

	// 转入地平纬度及地月质心距离,返回站心视半径(角秒)
	long double moonRad(long double r, long double h);

	// 求月亮近点时间和距离,t为儒略世纪数力学时
	Vector2 moonMinR(long double t, bool min);

	Vector2 moonNode(long double t, bool asc);

	// 地球近远点
	Vector2 earthMinR(long double t, bool min);

	Vector2 daJu(int xt, long double t, bool dx);

	Vector2 xingHR(int xt, long double t, bool f);

	Vector4 xingSP(int xt, long double t, int n, long double w0, long double ts, long double tp);

	Vector2 xingHY(int xt, long double t);

	Vector4 xingMP(int xt, long double t, int n, long double E, Vector4 g);

	//========行星天象及星历=============
	long double xingJJ(int xt, long double t, int jing);

	std::string xingX(int xt, long double jd, double L, double fa);

	Vector3 xingLiu0(int xt, long double t, int n, long double gxs);

	long double xingLiu(int xt, long double t, bool sn);

};

//========日月食计算使用的一些函数=============
COORDP lineEar(Vector3 P, Vector3 Q, long double gst);
COORDP lineEll(long double x1, long double y1, long double z1, long double x2, long double y2, long double z2, long double e, long double r);
COORDP lineEar2(long double x1, long double y1, long double z1, long double x2, long double y2, long double z2, long double e, long double r, Vector3 I);
NODE lineOvl(long double x1, long double y1, long double dx, long double dy, long double r, long double ba);
NODE cirOvl(long double R, long double ba, long double R2, long double x0, long double y0);
_ECFAST ecFast(long double jd);

long double moonIll(long double t);
long double moonRad(long double r, long double h);

//=================================deltat T计算=====================================
long double dt_T(long double t);

long double pty_zty2(long double t);

// 考虑真太阳时,与pty_zty2实现一模一样，这里先这样，后面去掉
inline long double mst_ast(long double t)
{
	long double L = (1753470142 + 628331965331.8 * t + 5296.74 * t * t) / (1000000000 - 0.0) + PI;
	Vector3 z{};

	long double E = (84381.4088 - 46.836051 * t) / rad;
	z[0] = XL0_calc(0, 0, t, 5) + PI;
	z[1] = 0;
	z[2] = 0;

	z = llrConv(z, E);
	L = rad2rrad(L - z[0]);
	return L / (PI * 2);
}

// 精气
inline long double qi_accurate(long double W, bool astFlag = false, long double longitude = 120.0)
{
	long double d = XL::S_aLon_t(W);
	long double t = d * 36525;
	if (astFlag)
	{
		return t - dt_T(t) + mst_ast(d) + longitude / 360.0;
	}
	else
	{
		return t - dt_T(t) + 8.0 / 24.0;
	}
}

inline long double so_accurate(long double W)
{
	long double t = XL::MS_aLon_t(W) * 36525;
	return t - dt_T(t) + 8.0 / 24.0;
} // 精朔

inline long double qi_accurate2(long double jd, bool astFlag = false, long double longitude = 120.0)
{
	long double d = PI / 12.0;
	long double w = floor((jd + 293) / 365.2422 * 24) * d;

	long double a = qi_accurate(w, astFlag, longitude);
	if ((a - jd) > 5)
	{
		return qi_accurate(w - d, astFlag, longitude);
	}
	else if ((a - jd) < -5)
	{
		return qi_accurate(w + d, astFlag, longitude);
	}
	else
	{
		return a;
	}
}

// 精朔
inline long double so_accurate2(long double jd)
{
	return so_accurate(floor((jd + 8) / 29.5306) * PI * 2);
}
