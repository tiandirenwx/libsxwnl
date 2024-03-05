#pragma once
#include <string>
#include "JD.h"
#include "geo.h"
/*
日食时间周期定义：
  食年 346.620 eclipse_year
  沙罗 6585.32 saros
  朔望 29.5306 syzygy
*/
#define ECLPSE_YEAR 346.620
#define SAROS       6585.32
#define SYZYGY      29.5306


std::string rysCalc(Time &d, bool is_utc, bool nasa_r,const JINGWEI &jw = {116.3833333, 39.900, "默认", "北京"});
std::string rs_search(int Y,int M,int n,bool fs);

// 不同周期下的日食概略推算
// fs 功能选择（function select）1
// jd0 初始儒略日
// step 步长
std::string rs2_calc(uint8_t fs,long double jd0, long double step = SYZYGY);

// 打印界线图
std::string rs2_jxb();

// 升降相关计算
std::string shengjiang(int y, int m, int d,const JINGWEI &jw = {116.3833333, 39.900, "默认", "北京"});
std::string shengjiang2(int y,const JINGWEI &jw = {116.3833333, 39.900, "默认", "北京"});
std::string shengjiang3(int y);