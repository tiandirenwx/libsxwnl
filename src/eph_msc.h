#pragma once
#include <string>
#include "eph.h"

typedef struct EphMoonSunCalcData
{
	long double T;	 // TD力学时
	long double L;	 // 地理纬度
	long double fa;	 // 地理经度
	long double dt;	 // 力学-世界时时差
	long double jd;	 // UT世界时
	long double dL;	 // 黄经章动
	long double dE;	 // 黄纬章动
	long double E;	 // 交角章动
	long double gst; // 真恒星时

	long double mHJ;   // 月球视黄经
	long double mHW;   // 月球视黄纬
	long double mR;	   // 地月质心距
	long double mCJ;   // 月球视赤经
	long double mCW;   // 月球视赤纬
	long double mShiJ; // 月球时角

	long double mCJ2; // 时差修正后的赤道坐标
	long double mCW2;
	long double mR2;
	long double mDJ; // 高度角
	long double mDW; // 方位角
	long double mPJ; // 大气折射修正后的高度角
	long double mPW; // 大气折射修正后的方位角

	long double sHJ;  // 太阳视黄经
	long double sHW;  // 太阳视黄纬
	long double sCJ;  // 太阳视赤经
	long double sCW;  // 太阳视赤纬
	long double sCJ2; // 时差修正后的赤道坐标
	long double sCW2;
	long double sR2;
	long double sShiJ; // 太阳时角

	long double sDJ; // 高度角
	long double sDW; // 方位角
	long double sR;
	long double sPJ; // 方位角
	long double sPW; // 高度角
	long double sc;	 // 时差

	long double pty;	  // 平恒星时
	long double zty;	  // 真恒星时
	long double mRad;	  // 月亮视半径
	long double sRad;	  // 太阳视半径
	long double e_mRad;	  // 月亮地心视半径
	long double eShadow;  // 地本影在月球向径处的半径(角秒)
	long double eShadow2; // 地半影在月球向径处的半径(角秒)
	long double mIll;	  // 月面被照亮比例
	long double zx_J;	  // 中心食坐标
	long double zx_W;
} S_EPH_MSC_DATA;

class EphMSC
{
public:
    const long double Eps = 1e-12;
	EphMSC();
	~EphMSC();
	// 参数：T是力学时,站点经纬L,fa,海拔high(千米)
	void calc(long double T, long double L,long double fa,long double high);
	std::string toString(bool flag);
	S_EPH_MSC_DATA &getData();

private:
	S_EPH_MSC_DATA *pstMscData;
};
