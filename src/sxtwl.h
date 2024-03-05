#pragma once
#include <cstdint>
#include <functional>
#include "day.h"

struct SLunarDay
{
	SLunarDay() {};
	int year;
	int monthIdx;
	int dayIdx;
	int lyGanIdx;
	int lyZhiIdx;
	bool isLeap;
	bool isNext;
};

struct SSolarDay
{
	SSolarDay() {};
	SSolarDay(int year, int month, int day): year(year),month(month),day(day)
	{

	}
	int year;
	int month;
	int day;
};

namespace sxtwl
{
	Day *fromSolar(int year, uint8_t month, int day);
	Day *fromLunar(int year, uint8_t month, int day, bool isRun = false,bool isSpec = false);
    Day *fromMoslem(int year, uint8_t month, int day);
    SLunarDay solar2Lunar(int _year, uint8_t _month, int _day, int _hour,bool isLunarNewYear = false );
	SLunarDay solar2Lunar2(const Time &t);
	//通过四柱获取年月日, 返回的是儒略日列表
	std::vector<double> siZhu2Year(GZ year, GZ  yue, GZ  ri, GZ  shi, int fromYear, int  toYear);
	//获取时辰上的那个天干
    GZ  getShiGz(uint8_t dayTg,  uint8_t hour, bool isZaoWanZiShi = true);
	//获取一年中的润月(不存，则返回0)
	uint8_t getRunMonth(int By);
	//获取一月中的阴日数量 
	uint8_t getLunarMonthNum(int By, uint8_t month, bool isRun = false,bool isSpec = false);
	//儒略日数转公历
	Time JD2DD(double jd);
	//公历转儒略日
	double toJD(Time& time);
};
