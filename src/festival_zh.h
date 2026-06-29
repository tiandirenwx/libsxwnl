#pragma once
#include <cstdint>
#include <string>

class Day;  // forward decl: 上层 Day 入口 (见 day.h)

// ═══════════════════════════════════════════════════════════════
//  festival_zh —— 中文公历/农历节日 + 周节日 + 杂节(数九/三伏/入梅/出梅)
//
//  数据与上游 sxwnl/lunar.js 中 oba.sFtv/wFtv 与 obb.getDayName 一致。
//  采用纯函数 + 静态数据表的设计, 不依赖天文计算。
//
//  两种调用方式 (语义完全等价):
//    1. festival::computeAll(day)            —— 直接传 Day, 内部抽取上下文
//    2. festival::computeAll(ctx)            —— 直接传 POD 上下文(便于单元测试)
// ═══════════════════════════════════════════════════════════════

namespace festival
{

// 节日分级: 与 lunar.js 一致, 用单字符标识
// '#' 放假节日(A 类)、'I' 重要节日(B 类)、'.' 一般节日(C 类)
constexpr char TIER_HOLIDAY = '#';
constexpr char TIER_MAJOR   = 'I';
constexpr char TIER_MINOR   = '.';

// ── 输入: 已经计算好的某一日上下文 ───────────────────────────
struct DayContext
{
    // 公历
    int year;          // 公历年(天文纪年, 可为负或 0)
    int month;         // 公历月 1-12
    int day;           // 公历日 1-31
    int weekday;       // 0-6, 0=星期日
    int weekiInMonth;  // 本日所在的周序号(0 起算)
    int firstWeekday;  // 本月 1 日的星期 (0-6)
    int weekTotal;     // 本月的总周数

    // 农历
    int  lunarMonth;        // 1-12, 1=正月, 12=腊月
    int  lunarDay;          // 1-30
    int  lunarMonthDays;    // 当前农历月天数(29 或 30)
    bool isLunarLeap;       // 当前为闰月
    bool nextMonthIsZheng;  // 下一农历月是否为"正月"(用于除夕/小年判定)

    // 节气
    int jieQiIdx;      // -1 表示当日无节气; 0-23 与 Jqmc 对应

    // 距关键节气的天数(整数)
    int curDz;         // 距冬至
    int curXz;         // 距夏至
    int curLq;         // 距立秋
    int curMz;         // 距芒种
    int curXs;         // 距小暑

    // 日干支索引
    int dayGanIdx;     // 0-9 (甲乙丙丁戊己庚辛壬癸)
    int dayZhiIdx;     // 0-11
};

// ── 输出: 该日聚合后的节日信息 ────────────────────────────────
struct DayInfo
{
    std::string holiday;   // A 类节日(置顶/置红)
    std::string major;     // B 类节日
    std::string minor;     // C 类节日
    std::string misc;      // 杂节: 数九/三伏/入梅/出梅
    bool        isOffDay;  // 法定放假(国家节假日 + 周末)

    DayInfo() : isOffDay(false) {}
};

// ── 主入口 ────────────────────────────────────────────────────
// 输入已计算好的 DayContext, 返回完整节日聚合
DayInfo computeAll(const DayContext &ctx);

// 直接对 Day 求值. 等价于 computeAll(fromDay(d)).
DayInfo computeAll(Day &d);

// 由 Day 对象抽取节日计算所需的上下文
DayContext fromDay(Day &d);

// ── 分项接口 (可单独调用, 便于扩展/测试) ──────────────────────
// 1. 公历定日节日
void appendSolarFixed(const DayContext &ctx, DayInfo &out);
// 2. 公历周节日(本月第 N 个星期 X)
void appendSolarWeekly(const DayContext &ctx, DayInfo &out);
// 3. 农历节日 + 节气作为节日
void appendLunarFixed(const DayContext &ctx, DayInfo &out);
// 4. 杂节: 数九/三伏/入梅/出梅
void appendMisc(const DayContext &ctx, DayInfo &out);

} // namespace festival
