#pragma once

#include <cstdint>
#include <algorithm>
#include <vector>
#include "JD.h"
#include "const.h"
#include "SSQ.h"

//static SSQ SSQPtr;

struct GZ
{
    GZ(){};
    GZ(uint8_t tg, uint8_t dz) : tg(tg), dz(dz)
    {
    }
    uint8_t tg;
    uint8_t dz;
};

class Day
{
private:
    int d0; //儒略日整数，对应12点
    double jdF;  //儒略日小数部分
    int bd0;   //公历月首,中午
    int bdn;   //本月天数(公历)

    int solar_year_;       //公历年
    uint8_t solar_month_;  //公历月
    int solar_day_;        //公历日

    int moslem_year_;      //回历年
    uint8_t moslem_month_; //回历月
    int moslem_day_;       //回历日

    int lunar_month_;           //阴历月的月
    uint8_t lunar_day_;         //阴历月的日
    int lunar_total_days_;      //该阴历月的总天数
    bool is_lunar_leap_month_;  //是不是阴历的润月
    bool is_lunar_spec_month_;   //是不是农历的下一个重月
    int cur_dz_;  
    int cur_xz_;
    int cur_lq_;
    int cur_mz_;
    int cur_xs_;

    int lunar_lichun_year_;  //以立春为界， 定农历纪年(10进制,1984年起算)
    int lunar_jun_year_;     //以正月初一为界，农历纪年(10进制,1984年起算)
    int huangdi_year_;       //黄帝纪年

    uint8_t week_;      //当前日的星期,星期几
    uint8_t week0_;     //本月第一天的星期
    uint8_t weeki_;     //本日所在的周序号
    uint8_t weekn_;     //本月的总周数

    int8_t yx_idx_;    // 月相索引
    double yxjd_;   	// 月相时刻(儒略日)
   
    int8_t jieling_;    //节令值
    uint8_t xingzuo_;   //星座
    long double jieqi_jd_;   //节气最体的时间

    uint8_t  lunar_day_12jian_;

    GZ *gz_lichun_year_;  //干支纪年(立春)
    GZ *gz_jan_year_;  //干支纪年(正月 春节)
    
    GZ *gz_month_; //月天干地支
    GZ *gz_day_;   // 日天干地支

private:
    explicit Day(int d,double jf = 0);
    void checkSSQ();
    void checkLunarData();
    void checkSolarData();
    void checkExtSolarData();
    void checkMoslemData();
    void checkJQData();
    void checkYXData();
    static Day *getLunarDateSpecial(int year, uint8_t month, int day, bool isRun, bool isSpec);
    static Day *getLunarDate(int year, uint8_t month, int day, bool isRun, bool isSpec);

public:
    virtual ~Day();
    Day *after(int day);
    Day *before(int day);
    
    double getJulianDate();
    // 获取阴历日期
    int getLunarDay();
    //获取阴历月
    uint8_t getLunarMonth();
    int getHuangdiYear();   //获取黄帝纪年
    //获取阴历年 chineseNewYearBoundary是否以春节为界
    int getLunarYear(bool chineseNewYearBoundary = true);
    //获取阴历年干支 chineseNewYearBoundary是否以春节为界
    GZ getYearGZ(bool chineseNewYearBoundary = false);
    GZ getMonthGZ();
    GZ getDayGZ();
    // 注意非早晚子时的时候，day要算第二天
	GZ getHourGZ(uint8_t hour, bool isZaoWanZiShi = true);
    bool isLunarLeap();
    bool isSpecNextMonth();
    uint8_t getLunarDay12Jian();

    //获取公历时间
    int getSolarYear();
    uint8_t  getSolarMonth();
    int getSolarDay();

    //获取回历
    int getMoslemYear();
    uint8_t  getMoslemMonth();
    int getMoslemDay();

    uint8_t getWeek();
    // 处于该月的第几周
    uint8_t getWeekIndex();
    uint8_t getFirstWeekDayOfMonth();
    uint8_t getTotalWeekNumsOfMonth();

    //查找月相
    bool hasYueXiang();
    uint8_t getYueXiang();
    double getYueXiangJD();

    //是否有节气
    bool hasJieQi();
    // 获取节气
    uint8_t getJieQi() ;
    double getJieQiJD();
    // 获取星座
    uint8_t  getConstellation();


public:
    static Day *fromDate(const Time &t);

    static Day *fromJulianDay(double jd);
    
    static Day *fromSolar(int _year, uint8_t _month, int _day);
    
    static Day *fromLunar(int year, uint8_t month, int day, bool isRun = false, bool isSpec = false);
    
    static Day *fromMoslem(int _year, uint8_t _month, int _day);

    //static Day *solarToLunar(int _year, uint8_t _month, int _day,int hour = 12);
    
    //static Day *lunarToSolar(int _year, uint8_t _month, int _day, bool _isRun = false, bool _isSpec = false);
};
