#pragma once

#include <cstdint>
#include <algorithm>
#include <vector>
#include <string>
#include "JD.h"
#include "const.h"
#include "SSQ.h"
#include "festival_zh.h"

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
    long double jdF;  //儒略日小数部分
    int bd0;   //公历月首,中午
    int bdn;   //本月天数(公历)

    int solar_year_;       //公历年
    uint8_t solar_month_;  //公历月
    int solar_day_;        //公历日

    int moslem_year_;      //回历年
    uint8_t moslem_month_; //回历月
    int moslem_day_;       //回历日

    int lunar_month_;           //阴历月的月名索引(SSQ.ym 中的值, 0=正 1=二 ... 11=腊)
    int lunar_month_idx_;       //阴历月在 SSQ 表内的序号(mk), -1 = 未计算
    uint8_t lunar_day_;         //阴历月的日
    int lunar_total_days_;      //该阴历月的总天数
    bool is_lunar_leap_month_;  //是不是阴历的润月
    bool is_lunar_spec_month_;   //是不是农历的下一个重月
    int lunar_month_style_;      //阴历月显示风格(LunarMonthNameStyle)
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
    long double yxjd_;   	// 月相时刻(儒略日)
   
    int8_t jieling_;    //节令值(天文口径, qi_accurate 精确交气时刻所在日)
    int8_t lipu_jq_idx_; //历谱口径节气索引(整日表+QB定气修正, 对应 sxwnl ob.Ljq); -2未算 -1无 0-23
    uint8_t xingzuo_;   //星座
    long double jieqi_jd_;   //节气最体的时间

    uint8_t  lunar_day_12jian_;

    GZ *gz_lichun_year_;  //干支纪年(立春)
    GZ *gz_jan_year_;  //干支纪年(正月 春节)
    
    GZ *gz_month_; //月天干地支
    GZ *gz_day_;   // 日天干地支

    // 节日聚合缓存
    bool                festivals_cached_ = false;
    festival::DayInfo   festivals_cached_value_;

private:
    explicit Day(int d,long double jf = 0);
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
    
    long double getJulianDate();
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
    //阴历月显示风格(LunarMonthNameStyle), 供显示层统一翻译月名
    int getLunarMonthStyle();
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
    long double getYueXiangJD();

    //是否有节气(天文口径: qi_accurate 精确交气时刻所在日)
    bool hasJieQi();
    // 获取节气
    uint8_t getJieQi() ;
    long double getJieQiJD();

    // 历谱口径节气(整日表 + QB 定气修正, 对应权威 sxwnl 网页版日历标注/ob.Ljq)
    // 与天文口径 getJieQi*() 的差异: 古代(1645年以前)因 QB 修正, 两者所在日可能相差 1 天;
    // 1645年及以后二者一致。月历"节气日"标注应使用本口径以对齐权威 sxwnl。
    bool hasLiPuJieQi();
    uint8_t getLiPuJieQi();
    // 获取星座
    uint8_t  getConstellation();

    // ── 数九/三伏/入梅/出梅 用到的距气天数 ──
    int getCurDz();  // 距冬至天数(负数表示尚未到)
    int getCurXz();  // 距夏至
    int getCurLq();  // 距立秋
    int getCurMz();  // 距芒种
    int getCurXs();  // 距小暑

    // 当前阴历月的总天数(29 或 30)
    int getLunarMonthDays();

    // 下一阴历月是否为"正月"(除夕/小年判定)
    bool isNextLunarMonthZheng();

    // ── 节日 ──
    // 计算并返回该日完整节日聚合 (首次调用后缓存)
    festival::DayInfo getFestivalInfo();

    // ── 派生显示字段(单一事实源, 各上层模块/UI 直接复用) ──
    std::string getLunarMonthName();   // "正月"/"闰二月"
    std::string getLunarDayName();     // "初一".."三十"
    std::string getJieQiName();        // "" 或 "冬至" (天文口径)
    std::string getJieQiTimeStr();     // "" 或 "HH:MM:SS" (天文口径)
    std::string getLiPuJieQiName();    // "" 或 "冬至" (历谱口径, 整日表+QB, 对齐权威 sxwnl 网页版)
    std::string getYueXiangName();     // "" 或 "朔/上弦/望/下弦"
    std::string getYueXiangTimeStr();  // "" 或 "HH:MM:SS"
    std::string getConstellationName();// "白羊座"
    std::string getShengXiao();        // "马"


public:
    static Day *fromDate(const Time &t);

    static Day *fromJulianDay(long double jd);
    
    static Day *fromSolar(int _year, uint8_t _month, int _day);
    
    static Day *fromLunar(int year, uint8_t month, int day, bool isRun = false, bool isSpec = false);
    
    static Day *fromMoslem(int _year, uint8_t _month, int _day);

    //static Day *solarToLunar(int _year, uint8_t _month, int _day,int hour = 12);
    
    //static Day *lunarToSolar(int _year, uint8_t _month, int _day, bool _isRun = false, bool _isSpec = false);
};
