#pragma once

#include <memory>
#include <string>
#include <vector>
#include "day.h"

// ═══════════════════════════════════════════════════════════════
//  MonthCalendar —— 公历/农历/回历 综合月历类
//
//  对应上游 sxwnl/lunar.js 中的 Lunar.yueLiCalc(By, Bm)。
//
//  设计原则: 不重复 Day 类已有的数据与逻辑。每个 CalendarDay 仅持有一
//  个 Day 智能指针, 所有字段(公历/农历/回历/干支/节气/月相/星座/距气/
//  节日) 通过 day->getXxx() 直接获取, 这样 sxwnl 内核中"日历单日"逻辑
//  只有一份且经过验证。
//
//  用法:
//    sxwnl::MonthCalendar cal(2026, 6);
//    for (const auto& cd : cal.days()) {
//        cd.day->getSolarDay();
//        cd.day->getLunarMonthName();
//        cd.day->getJieQiName();
//        cd.day->getFestivalInfo();   // 已自动缓存
//    }
// ═══════════════════════════════════════════════════════════════

namespace sxwnl
{

struct CalendarDay
{
    std::shared_ptr<Day> day;       // 所有字段经此查询; 经过 Day 已验证的逻辑
    int                  indexInMonth;  // 0-based, 公历月内第几天
};

class MonthCalendar
{
public:
    // 一次性计算整月. solarMonth ∈ [1, 12]
    // throws std::invalid_argument on bad input
    MonthCalendar(int solarYear, int solarMonth);

    int year()           const { return year_; }
    int month()          const { return month_; }
    uint8_t firstWeek()  const { return w0_; }          // 月首星期 0..6
    uint8_t totalWeeks() const { return totalWeeks_; }
    int dayCount()       const { return (int)days_.size(); }
    int firstJulianDay() const { return d0_; }          // J2000 起算

    // 该月所属"农历年"信息(春节为界), 由首日 Day 派生
    GZ                  yearGZ()    const { return yearGZ_; }
    const std::string&  shengXiao() const { return shengXiao_; }
    const std::string&  nianHao()   const { return nianHao_; }

    const std::vector<CalendarDay>& days() const { return days_; }
    const CalendarDay& at(int idx) const { return days_.at(idx); }

private:
    int       year_, month_;
    uint8_t   w0_, totalWeeks_;
    int       d0_;
    GZ        yearGZ_;
    std::string shengXiao_;
    std::string nianHao_;
    std::vector<CalendarDay> days_;
};

} // namespace sxwnl
