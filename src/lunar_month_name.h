#pragma once
//
// 阴历月名统一生成器 —— 八字 / 年历 / 月历 / 日详情各处的"单一事实源"。
//
// 设计目的:
//   历史上帝王多次更改月建别名, 使得同一年可能出现两个"数字月名"相同、却需不同
//   写法的农历月。过去各显示层(capi/bazi/lunar_ob/day)各自用 (isLeap, isSpec, year)
//   重复推断名称, 逻辑分散且易错(如 690 年"建寅"被误显示成"拾壹月")。
//
//   现统一为: SSQ.calcY 计算出每月的
//       · 数字月名 nameNum(1..12: 正=1 … 十=10, 冬=11, 腊=12)
//       · 显示风格 style  (LunarMonthNameStyle, 见 SSQ.h)
//       · 是否闰月 isLeap
//   显示层只需调用 lunarMonthName() 即可得到最终字符串, 不再各自推断。
//
//   与上游寿星天文历 lunar.js 的对应关系:
//       · 普通月           -> Ymc[](正/二/…/十/冬/腊)          (JS ymc: 十一/十二/…)
//       · 连续同名月(CE)    -> SYmc[](拾贰月等)                  (JS mc='拾贰')
//       · 689-701 建寅      -> "一月"                            (JS mc='一')
//       · 秦汉年终置闰(闰)  -> "后九月"  (y∈[-220,-104])          (JS ns[i+3]='后九')
//       · 春秋战国年终置闰  -> "十三月"  (y∈[-721,-220))          (JS ns[i+3]='十三')
//
// 说明: is_spec(重月回环标记)只服务于农历↔公历互转, 与显示彻底解耦, 此处不使用。
//
#include <string>
#include "SSQ.h"        // LunarMonthNameStyle
#include "sx_lang_zh.h" // Ymc / SYmc / BDLeapYueName

namespace sxwnl
{

// year   : 公历年(仅用于古历"后九月/十三月"区间判定)
// nameNum : 数字月名 1..12
// style  : LunarMonthNameStyle
// isLeap : 是否闰月
// 古历年终置闰月名: 春秋/战国"十三月"(y∈[-721,-220)), 秦汉"后九月"(y∈[-220,-104])。
// 空串表示 year 不在古历置闰区间。
inline std::string ancientLeapMonthName(int year)
{
    if (year >= -721 && year < -220) return std::string(BDLeapYueName[0]) + "月"; // 十三月
    if (year >= -220 && year <= -104) return std::string(BDLeapYueName[1]) + "月"; // 后九月
    return {};
}

inline std::string lunarMonthName(int year, int nameNum, int style, bool isLeap)
{
    // 古历年终置闰月(一年至多一个): 优先按显示风格判定, 与 leap 下标解耦。
    // 之所以不只依赖 isLeap: leap_month_==0 兼作"本年无闰"哨兵, 当置闰月落在年历
    // 首月(下标 0, 如公元前 255 年"十三月")时 isLeap 会被漏判, 故以风格为准。
    if (style == LUNAR_MONTH_ANCIENT_LEAP)
    {
        std::string name = ancientLeapMonthName(year);
        if (!name.empty()) return name;
    }
    if (isLeap)
    {
        std::string name = ancientLeapMonthName(year);
        if (!name.empty()) return name;
    }

    std::string prefix = isLeap ? "闰" : "";
    std::string base;
    switch (style)
    {
    case LUNAR_MONTH_YI:
        base = "一"; // 689-701 建寅, 区别于建子改称的"正"
        break;
    case LUNAR_MONTH_SYMC:
        base = (nameNum >= 1 && nameNum <= 12) ? SYmc[nameNum - 1] : std::to_string(nameNum);
        break;
    case LUNAR_MONTH_NORMAL:
    default:
        base = (nameNum >= 1 && nameNum <= 12) ? Ymc[nameNum - 1] : std::to_string(nameNum);
        break;
    }
    return prefix + base + "月";
}

} // namespace sxwnl
