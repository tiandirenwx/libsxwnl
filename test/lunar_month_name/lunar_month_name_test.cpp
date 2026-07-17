// 阴历月名(尤其古历"十三月/后九月"与普通闰月)验证测试
//
// 背景 (结合 JS lunar.js 与 C++ 实现):
//   * JS lunar.js 的 calcY 中, 古历(YY∈[-721,-104])把年终置闰月的月名直接写成
//     字符串 "十三"(春秋/战国, ns[i+3]='十三') 或 "后九"(秦汉, ns[i+3]='后九'),
//     且该分支 return 前不设置 this.leap —— 即"十三/后九"是"月名", 一年至多一个。
//   * C++ SSQ 排月序把月序存为数字, 无法直接存字符串, 故显示层(capi lunar_month_str /
//     enum_lunar_window / bazi getLunar)按"该月是年终置闰月(isLeap) + 农历年落在
//     对应区间"来重建 "十三月"(y∈[-721,-220)) / "后九月"(y∈[-220,-104])。
//   * 关键: 只有那个 isLeap 的置闰月才是十三/后九, 不能对全年套用; 重月(isSpec,
//     如"后一个十月/冬月")仍走正常/ SYmc 逻辑, 不应被误标成十三/后九。
//
// 本测试通过导出的 C API (sxwnl_get_day_info) 端到端验证真实显示逻辑。
//
// 编译(CMake): make lunar_month_name_test  ;  运行: ./build/bin/lunar_month_name_test

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "sxwnl_capi.h"
#include "SSQ.h"
#include "const.h"

namespace {

int g_pass = 0;
int g_fail = 0;

void check(bool cond, const std::string &msg) {
    if (cond) {
        ++g_pass;
    } else {
        ++g_fail;
        std::printf("  [FAIL] %s\n", msg.c_str());
    }
}

// 取某公历日的农历月名/日名; 无效日期返回空串
std::string monthName(int y, int m, int d) {
    SxwnlDayInfo info;
    if (sxwnl_get_day_info(y, m, d, &info) != 0) return {};
    return std::string(info.lunar_month_name);
}
std::string dayName(int y, int m, int d) {
    SxwnlDayInfo info;
    if (sxwnl_get_day_info(y, m, d, &info) != 0) return {};
    return std::string(info.lunar_day_name);
}

// ─────────────────────────────────────────────
// Part A: 现代已知日期 —— 普通月 / 普通闰月 端到端
// ─────────────────────────────────────────────
void testModern() {
    std::printf("== Part A: 现代普通月/闰月 ==\n");
    struct Case { int y, m, d; const char *mon; const char *day; };
    const std::vector<Case> cases = {
        {2024, 2, 10, "正月",   "初一"}, // 甲辰年春节
        {2024, 1,  1, "冬月",   "二十"}, // 癸卯年冬月(十一月)
        {2023, 2, 20, "二月",   "初一"}, // 闰二月之前的正常二月
        {2023, 3, 22, "闰二月", "初一"}, // 2023 闰二月
        {2020, 4, 23, "四月",   "初一"},
        {2020, 5, 23, "闰四月", "初一"}, // 2020 闰四月
        {2017, 7, 23, "闰六月", "初一"}, // 2017 闰六月
    };
    for (const auto &c : cases) {
        std::string mn = monthName(c.y, c.m, c.d);
        std::string dn = dayName(c.y, c.m, c.d);
        std::printf("  %04d-%02d-%02d -> %s%s (期望 %s%s)\n",
                    c.y, c.m, c.d, mn.c_str(), dn.c_str(), c.mon, c.day);
        check(mn == c.mon, std::string("月名 ") + c.mon);
        check(dn == c.day, std::string("日名 ") + c.day);
    }
}

// 统计某公历年内出现的各农历月名(按天计数)
std::map<std::string, int> collectYear(int year) {
    std::map<std::string, int> cnt;
    for (int m = 1; m <= 12; ++m) {
        for (int d = 1; d <= 31; ++d) {
            std::string mn = monthName(year, m, d);
            if (!mn.empty()) cnt[mn]++;
        }
    }
    return cnt;
}

// ─────────────────────────────────────────────
// Part B: 古历特殊月 —— 十三月 / 后九月
//   * 出现时必须是正确朝代名 (十三: y∈[-721,-220); 后九: y∈[-220,-104])
//   * 单年内最多一个置闰月 (≤31 天), 绝不能全年皆为十三/后九
// ─────────────────────────────────────────────
void testAncientRange(int y0, int y1, const char *expectSpecial, const char *forbidSpecial) {
    std::printf("== Part B: 古历 [%d..%d] 期望特殊名=%s, 禁止=%s ==\n",
                y0, y1, expectSpecial, forbidSpecial);
    int yearsWithSpecial = 0;
    for (int year = y0; year <= y1; ++year) {
        auto cnt = collectYear(year);
        int specialDays = 0, normalDays = 0, forbidDays = 0;
        std::string names;
        for (const auto &kv : cnt) {
            names += kv.first + "(" + std::to_string(kv.second) + ") ";
            if (kv.first == expectSpecial) specialDays += kv.second;
            else if (kv.first == forbidSpecial) forbidDays += kv.second;
            else normalDays += kv.second;
        }
        std::printf("  %d: %s\n", year, names.c_str());
        // 禁止出现别朝代的特殊名
        check(forbidDays == 0,
              std::string("年 ") + std::to_string(year) + " 不应出现 " + forbidSpecial);
        // 特殊名若出现, 必须被限制在单个置闰月内 (≤31 天), 不能覆盖全年
        check(specialDays <= 31,
              std::string("年 ") + std::to_string(year) + " 十三/后九月天数应≤31, 实为 "
              + std::to_string(specialDays));
        // 出现特殊名的年份, 普通月必须占绝大多数 (证明不是全年套用)
        if (specialDays > 0) {
            ++yearsWithSpecial;
            check(normalDays >= specialDays * 5,
                  std::string("年 ") + std::to_string(year)
                  + " 普通月天数应远多于特殊月");
        }
    }
    // 该区间内至少应有若干年出现特殊名, 证明"十三/后九"确实被正确生成
    check(yearsWithSpecial > 0,
          std::string("区间 [") + std::to_string(y0) + ".." + std::to_string(y1)
          + "] 应至少有一年出现 " + expectSpecial);
    std::printf("  -> 出现 %s 的年份数: %d\n", expectSpecial, yearsWithSpecial);
}

// ─────────────────────────────────────────────
// Part C: 纪年字符串 → 天文年 (对齐上游 tools.js year2Ayear)
//   * B/b/* 及 公元前/前  => 公元前记法, 天文 = 1 - n
//   * 前导 '-'            => 直接天文负年 (如 "-220" => 天文 -220)
//   * 公元/纯数字          => 直接天文年
//   历史 bug: '-' 曾被当公元前, "-220" 误解析为 -219 => 生肖差一位(龙误作蛇)
// ─────────────────────────────────────────────
void testYearParse() {
    std::printf("== Part C: 纪年字符串 -> 天文年 ==\n");
    struct YC { const char *s; int astro; };
    const std::vector<YC> cases = {
        {"-220",       -220}, // 关键回归: 直接天文年, 不是公元前220(-219)
        {"-1",           -1},
        {"220",         220},
        {"2024",       2024},
        {"2024年",     2024},
        {"B221",       -220}, // 公元前221 = 天文 -220
        {"b221",       -220},
        {"*221",       -220},
        {"公元前221",  -220},
        {"前221",      -220},
        {"公元2024",   2024},
        {"B1",            0}, // 公元前1年 = 天文 0
    };
    for (const auto &c : cases) {
        int got = sxwnl_year_str_to_astro(c.s);
        std::printf("  \"%s\" -> %d (期望 %d)\n", c.s, got, c.astro);
        check(got == c.astro, std::string("纪年解析 ") + c.s);
    }
}

// ─────────────────────────────────────────────
// Part D: 生肖端到端 —— 天文年 -220 应为龙(辰), -219 为蛇(巳)
//   通过 sxwnl_get_day_info 取立春后某日的年生肖, 印证 Part C 修复的实际效果
// ─────────────────────────────────────────────
void testShengXiao() {
    std::printf("== Part D: 生肖 (立春后取样) ==\n");
    struct SC { int astroYear; const char *sx; };
    const std::vector<SC> cases = {
        {-220, "龙"}, // 鸿蒙正确值
        {-219, "蛇"}, // 旧 Android bug 的错误值, 对应天文 -219
        {2024, "龙"}, // 甲辰
        {2025, "蛇"}, // 乙巳
    };
    for (const auto &c : cases) {
        SxwnlDayInfo info;
        // 取公历 6/1, 稳在立春之后、次年立春之前
        int rc = sxwnl_get_day_info(c.astroYear, 6, 1, &info);
        std::string sx = (rc == 0) ? std::string(info.sheng_xiao) : std::string("<err>");
        std::printf("  天文年 %d 6-1 -> 生肖 %s (期望 %s)\n", c.astroYear, sx.c_str(), c.sx);
        check(sx == c.sx, std::string("生肖 年") + std::to_string(c.astroYear));
    }
}

// ─────────────────────────────────────────────
// Part E: 具体已知日期锚点 (端到端: 年干支 + 月名 + 日名)
//   -220-10-19 (天文年) = 庚辰年 后九月 初一  —— 秦汉年终置闰月
// ─────────────────────────────────────────────
void testKnownDay() {
    std::printf("== Part E: 已知日期锚点 ==\n");
    struct DC { int y, m, d; const char *gz, *mon, *day; };
    const std::vector<DC> cases = {
        {-220, 10, 19, "庚辰", "后九月", "初一"},
        {-220, 10, 18, "庚辰", "九月",   "廿九"}, // 后九月前一日
        {-220, 11, 17, "庚辰", "后九月", "三十"}, // 后九月末日
        // 回归: "十三月"恰好落在年历首月(下标0), leap_month_==0 与"本年无闰"哨兵
        //       冲突, 曾被误显示为"腊月"。现以显示风格承载, 应正确显示"十三月"。
        {-256, 11, 25, "甲辰", "十三月", "初一"},
        {-256, 12, 24, "甲辰", "十三月", "三十"}, // 十三月末日
        {-256, 12, 25, "甲辰", "正月",   "初一"}, // 次月恢复正月
    };
    for (const auto &c : cases) {
        SxwnlDayInfo info;
        int rc = sxwnl_get_day_info(c.y, c.m, c.d, &info);
        std::string gz  = (rc == 0) ? std::string(info.year_gz) : "<err>";
        std::string mon = (rc == 0) ? std::string(info.lunar_month_name) : "<err>";
        std::string day = (rc == 0) ? std::string(info.lunar_day_name) : "<err>";
        std::printf("  %d-%02d-%02d -> %s年 %s%s (期望 %s年 %s%s)\n",
                    c.y, c.m, c.d, gz.c_str(), mon.c_str(), day.c_str(),
                    c.gz, c.mon, c.day);
        check(gz == c.gz && mon == c.mon && day == c.day,
              std::string("锚点 ") + std::to_string(c.y) + "-"
              + std::to_string(c.m) + "-" + std::to_string(c.d));
    }
}

// ─────────────────────────────────────────────
// Part F: 年历 API (sxwnl_get_year_calendar) —— 重月显示规则回归
//   * 回归 "-220 首月误显示拾壹月(应为冬月)" bug。
//   * 古历区间(y∈[-721,-104])重月: 显示走正常农历名(十月/冬月), 但 is_spec 标记
//     必须保留(供农历↔公历"重月"转换回环)。特殊名只保留 后九月/十三月。
//   * CE 纪年因月建别名变更产生的连续同名月(如 240 年): 仍显示 "拾贰月"(SYmc),
//     与上游 lunar.js 一致。
// ─────────────────────────────────────────────
struct YCScan {
    int total = 0;
    std::vector<SxwnlYearCalMonth> m;
};
YCScan yearCal(int year) {
    YCScan s;
    s.m.resize(16);
    s.total = sxwnl_get_year_calendar(year, s.m.data(), 16);
    if (s.total < 0) s.total = 0;
    s.m.resize(s.total);
    return s;
}
int countName(const YCScan &s, const char *name) {
    int c = 0;
    for (const auto &e : s.m)
        if (std::string(e.month_name) == name) c++;
    return c;
}
const SxwnlYearCalMonth *findName(const YCScan &s, const char *name) {
    for (const auto &e : s.m)
        if (std::string(e.month_name) == name) return &e;
    return nullptr;
}

void testYearCalendar() {
    std::printf("== Part F: 年历 API 月名/标记 ==\n");

    // -220: 关键回归 —— 首月应"冬月"(非"拾壹月"), 且 is_spec 标记保留
    {
        YCScan s = yearCal(-220);
        check(s.total > 0, "年历 -220 应有数据");
        if (s.total > 0) {
            std::string first = s.m[0].month_name;
            std::printf("  -220 首月: %s (is_spec=%d) [期望 冬月, is_spec=1]\n",
                        first.c_str(), s.m[0].is_spec);
            check(first == "冬月", "年历 -220 首月应为 冬月(非 拾壹月)");
            check(s.m[0].is_spec == 1, "年历 -220 首月应保留 is_spec 标记(转换回环)");
        }
        check(countName(s, "拾壹月") == 0, "年历 -220 不应出现 拾壹月");
        check(countName(s, "拾月") == 0, "年历 -220 不应出现 拾月");
        const SxwnlYearCalMonth *hj = findName(s, "后九月");
        check(hj != nullptr && hj->is_leap == 1, "年历 -220 应含闰 后九月");
    }

    // -221: 两个"十月"(显示同名), 恰一个 is_spec; 含"冬月"; 不出现 拾/拾壹
    {
        YCScan s = yearCal(-221);
        int shiCnt = countName(s, "十月");
        int specShi = 0;
        for (const auto &e : s.m)
            if (std::string(e.month_name) == "十月" && e.is_spec) specShi++;
        std::printf("  -221 十月出现 %d 次, 其中 is_spec %d 个 [期望 2 / 1]\n",
                    shiCnt, specShi);
        check(shiCnt == 2, "年历 -221 应有两个 十月");
        check(specShi == 1, "年历 -221 两个十月中恰有一个 is_spec(转换回环)");
        check(countName(s, "冬月") >= 1, "年历 -221 应含 冬月");
        check(countName(s, "拾月") == 0 && countName(s, "拾壹月") == 0,
              "年历 -221 不应出现 拾/拾壹月");
    }

    // 240 (CE): 连续同名月 → 保留 拾贰月(SYmc); 相邻 腊月 存在
    {
        YCScan s = yearCal(240);
        const SxwnlYearCalMonth *sp = findName(s, "拾贰月");
        std::printf("  240 含拾贰月: %s\n", sp ? "是" : "否");
        check(sp != nullptr && sp->is_spec == 1, "年历 240 应含 拾贰月(is_spec=1)");
        check(countName(s, "腊月") >= 1, "年历 240 应含 腊月");
    }
}

// ─────────────────────────────────────────────
// Part G: 重月端到端 (逐日月名/日名)
//   * 古历秦汉"后十月/后冬月": 显示正常农历名(十月/冬月), 绝不出现 拾/拾壹。
//   * CE 纪年连续十二月: 后一个显示 拾贰月。
// ─────────────────────────────────────────────
void testRepeatedMonths() {
    std::printf("== Part G: 重月端到端 (古历正常名 / CE 拾贰月) ==\n");
    struct DC { int y, m, d; const char *mon, *day; };
    const std::vector<DC> cases = {
        {-221, 10,  1, "冬月",   "初一"}, // 冬月(前一个)
        {-221, 10, 31, "十月",   "初一"}, // 后一个十月 —— 正常名
        {-221, 11, 29, "冬月",   "初一"}, // 后一个冬月 —— 正常名(回归 拾壹月 bug 的日级视角)
        { 23,  12, 31, "拾贰月", "初一"}, // CE 后一个十二月
        { 24,   1, 12, "拾贰月", "十三"}, // 同一 拾贰月 内
    };
    for (const auto &c : cases) {
        std::string mn = monthName(c.y, c.m, c.d);
        std::string dn = dayName(c.y, c.m, c.d);
        std::printf("  %d-%02d-%02d -> %s%s (期望 %s%s)\n",
                    c.y, c.m, c.d, mn.c_str(), dn.c_str(), c.mon, c.day);
        check(mn == c.mon && dn == c.day,
              std::string("重月锚点 ") + std::to_string(c.y) + "-"
              + std::to_string(c.m) + "-" + std::to_string(c.d));
    }
}

// ─────────────────────────────────────────────
// Part H: 与上游 lunar.js 逐年逐月对齐 (JS 忠实复刻为 oracle)
//
//   背景: 历史上帝王多次更改月建别名(见 SSQ.cpp 注释), 涉及:
//     · 8.01.15~23.12.02  建子为十二(顺推)
//     · 237.04.12~239.12.13 建子为十二(顺推)
//     · 23/239 年底连续两个十二月 → 拾贰月
//     · 689.12.18~700.11.15 建子为正, 建寅为"一月"   ← 690 拾壹月 bug 所在
//     · 761.12.02~762.03.30 建子为正(顺推), 762 出现重复的四月/五月
//   这里把 lunar.js calcY 的"名称转换(月建别名)"循环忠实复刻为 oracle, 用 SSQ 的
//   HS/ZQ 实数据驱动, 逐年比对 sxwnl_get_year_calendar 的月名 (多重集, 归一化
//   冬↔十一 / 腊↔十二)。多重集比对可容忍"连续同名月"的 拾贰 落在前/后月的差异
//   (C++ 落后一月, JS 落前一月, 两者月名集合一致)。
//   注: oracle 仅适用于非古历分支(年份 > -104); 古历十三月/后九月见 Part B。
// ─────────────────────────────────────────────
static const char *JS_ymc[12] = {
    "十一","十二","正","二","三","四","五","六","七","八","九","十"};

// 忠实复刻 lunar.js calcY 的"月名(月建别名)"生成, 用 SSQ 的实数据(HS 朔 / ZQ 气 /
// calc 朔气解算)驱动, 逐月还原上游会显示的月名字符串。返回该年年历窗口内每月月名
// (未加"月"字; 古历"后九/十三"直接给中文名, 便于与 C++ 归一化后比对)。
//
// 覆盖 SSQ.cpp 中出现的全部特殊情形:
//   · [-721,-104] 古历分支: 春秋/战国"十三"、秦汉"后九"年终置闰; 月建别名累积推算
//   · 8.01~23.12 / 237.04~239.12  建子为十二(顺推)
//   · 23/239 年底连续两个十二月 → 后一个"拾贰"
//   · 689.12~700.11  建子为正、建寅为"一"     ← 690"拾壹月"bug 区
//   · 761.12~762.03  建子为正(顺推)
std::vector<std::string> jsYearMonthNames(int year) {
    auto int2d = [](long double v) { return (long)std::floor((double)v); };
    SSQ &ssq = SSQ::getInstance();
    ssq.calcY((int)int2d((year - 2000.0) * 365.2422 + 180));
    std::vector<long double> A = ssq.getZhongQi(); // 25 气
    std::vector<int> B = ssq.getHS();              // 15 朔

    // 与 SSQ.calcY 一致地确定年份 YY (用于判定古历分支)
    int YY = (int)int2d((A[0] + 10 + 180) / 365.2422) + 2000;

    std::vector<std::string> out;

    // ── 古历分支 [-721,-104]: 十三月 / 后九月 ──────────────────────
    if (YY >= -721 && YY <= -104) {
        long ns[12] = {0,0,0,0,0,0,0,0,0,0,0,0};
        for (int i = 0; i < 3; i++) {
            long yy = YY + i - 1;
            if (yy >= -721) { // 春秋历, 年终置闰名"十三", 年首建正(建子=十一顺推, 首月正=建寅)
                ns[i]     = ssq.calc((long double)(1457698 - J2000) + int2d(0.342 + (yy + 721) * 12.368422) * 29.5306L, SType);
                ns[i + 3] = 1;  // 闰月月序(→ yueIdx[1]=腊, 但 isLeap 时显示"十三")
                ns[i + 6] = 2;  // 月建基准
            }
            if (yy >= -479) { // 战国历
                ns[i]     = ssq.calc((long double)(1546083 - J2000) + int2d(0.500 + (yy + 479) * 12.368422) * 29.5306L, SType);
                ns[i + 3] = 1;
                ns[i + 6] = 2;
            }
            if (yy >= -220) { // 秦汉历, 年终置闰名"后九", 年首建亥(十月)
                ns[i]     = ssq.calc((long double)(1640641 - J2000) + int2d(0.866 + (yy + 220) * 12.369000) * 29.5306L, SType);
                ns[i + 3] = 9;
                ns[i + 6] = 11;
            }
        }
        for (int i = 0; i < 14; i++) {
            if (B[i + 1] > A[24]) break;               // 与年历窗口一致
            int nn;
            for (nn = 2; nn >= 0; nn--) if (B[i] >= ns[nn]) break;
            int f1 = (int)int2d((B[i] - ns[nn] + 15) / 29.5306L); // 该月积数
            if (f1 < 12) {
                int v2 = (int)(((f1 + ns[nn + 6]) % 12 + 12) % 12);
                out.push_back(JS_ymc[v2]);
            } else {
                out.push_back(ns[nn + 3] == 9 ? "后九" : "十三"); // 年终置闰月
            }
        }
        return out;
    }

    // ── 现代分支(气朔置闰法) ──────────────────────────────────────
    std::vector<int> ym(14);
    for (int i = 0; i < 14; i++) ym[i] = i;
    if (B[13] <= A[24]) {
        int i;
        for (i = 1; B[i + 1] > A[2 * i] && i < 13; i++);
        for (; i < 14; i++) ym[i]--;
    }
    for (int i = 0; i < 14; i++) {
        if (B[i + 1] > A[24]) break; // 与年历窗口一致
        long double Dm = (long double)B[i] + J2000;
        int v2 = ym[i];
        std::string mc = JS_ymc[((v2 % 12) + 12) % 12];
        if (Dm >= 1724360 && Dm <= 1729794) mc = JS_ymc[(((v2 + 1) % 12) + 12) % 12];
        else if (Dm >= 1807724 && Dm <= 1808699) mc = JS_ymc[(((v2 + 1) % 12) + 12) % 12];
        else if (Dm >= 1999349 && Dm <= 1999467) mc = JS_ymc[(((v2 + 2) % 12) + 12) % 12];
        else if (Dm >= 1973067 && Dm <= 1977052) { if (v2 % 12 == 0) mc = "正"; if (v2 == 2) mc = "一"; }
        if (Dm == 1729794 || Dm == 1808699) mc = "拾贰";
        out.push_back(mc);
    }
    return out;
}

// 归一化: 去"闰"/"月", 冬→十一, 腊→十二
std::string normMonth(std::string s) {
    size_t p;
    while ((p = s.find("闰")) != std::string::npos) s.erase(p, 3);
    while ((p = s.find("月")) != std::string::npos) s.erase(p, 3);
    if (s == "冬") s = "十一";
    if (s == "腊") s = "十二";
    return s;
}

// 逐年(连续, 无遗漏)比对某区间, 汇总统计。special!=nullptr 时额外要求:
// 该区间内至少有一年 JS oracle 产出了 special 月名(证明特殊月确被覆盖到)。
void sweepOracle(const char *label, int y0, int y1, const char *special) {
    int pass = 0, fail = 0, withSpecial = 0;
    for (int y = y0; y <= y1; ++y) {
        if (y == 0) continue; // 无公元 0 年
        std::vector<std::string> js = jsYearMonthNames(y);
        std::vector<SxwnlYearCalMonth> m(20);
        int n = sxwnl_get_year_calendar(y, m.data(), 20);
        std::vector<std::string> a, b;
        for (int i = 0; i < n; ++i) a.push_back(normMonth(m[i].month_name));
        for (auto &s : js) b.push_back(normMonth(s));
        std::sort(a.begin(), a.end());
        std::sort(b.begin(), b.end());
        bool ok = (a == b);
        if (!ok) {
            std::string sa, sb;
            for (int i = 0; i < n; ++i) { sa += m[i].month_name; sa += " "; }
            for (auto &s : js) { sb += s + "月 "; }
            std::printf("  [DIFF] %d C++=[%s] JS=[%s]\n", y, sa.c_str(), sb.c_str());
        }
        if (ok) ++pass; else ++fail;
        check(ok, std::string("年历月名集合对齐 JS  年") + std::to_string(y));
        if (special) {
            for (auto &s : js) if (s == special) { ++withSpecial; break; }
        }
    }
    std::printf("  [%s] 逐年比对 [%d..%d]: 共 %d 年, 通过 %d, 失败 %d",
                label, y0, y1, pass + fail, pass, fail);
    if (special) {
        std::printf(", 其中出现\"%s\"的年份 %d 个", special, withSpecial);
        check(withSpecial > 0,
              std::string("区间 [") + std::to_string(y0) + ".." + std::to_string(y1)
              + "] 应至少一年出现 " + special);
    }
    std::printf("\n");
}

void testJsOracle() {
    std::printf("== Part H: 逐年与 lunar.js oracle 全覆盖对齐 (月名多重集) ==\n");
    // 古历: 每一年都比对 (十三月区间 + 后九月区间, 含改历分界)
    sweepOracle("春秋战国-十三月", -721, -221, "十三");
    sweepOracle("秦汉-后九月",     -220, -104, "后九");
    // CE 月建别名改动区(含 8-23 / 237-240 连续十二月 + 689-701 建寅一月 + 761-762)
    sweepOracle("CE改历-前段",       -103, 260, "拾贰");
    sweepOracle("武周/中唐改历",       680, 770, nullptr);
    // 现代广覆盖(每年)
    sweepOracle("中古-近代",          261, 679, nullptr);
    sweepOracle("近代-现代",          771, 2100, nullptr);
}

// ─────────────────────────────────────────────
// Part I: 690/701/762 修复回归 (端到端逐日) + 240/24 保持
//   690 建寅: 应"一月"(旧 bug: 壹月/拾壹月)
//   701 建丑(700.11.15 后已恢复): 应"腊月"(旧 bug: 拾贰月)
//   762 恢复旧历后的重复月: 正常"四月/五月"(旧 bug: 肆月/伍月)
//   240/24 连续十二月: 保持 拾贰月 (C++ 既定行为, 客户端依赖)
// ─────────────────────────────────────────────
void testReformMonthsFixed() {
    std::printf("== Part I: 690/701/762 月建别名修复回归 ==\n");
    struct DC { int y, m, d; const char *mon; };
    const std::vector<DC> cases = {
        {690, 1,  5, "正月"},   // 建子改称 正
        {690, 1, 20, "腊月"},   // 建丑(不变)
        {690, 2, 20, "一月"},   // 建寅改称 一 (核心回归: 非 壹月/拾壹月)
        {690, 3, 20, "二月"},
        {699, 2, 20, "一月"},   // 期末年份同样为 一月
        {700, 2, 20, "一月"},
        {701, 1, 20, "腊月"},   // 恢复后: 建丑 = 腊月 (非 拾贰月)
        {701, 2, 20, "正月"},   // 恢复后: 建寅 = 正月
        {762, 3,  5, "四月"},   // 恢复旧历: 重复四月(前)
        {762, 4,  5, "五月"},   // 重复五月(前)
        {762, 5,  5, "四月"},   // 重复四月(后) —— 正常名, 非 肆月
        {762, 6,  5, "五月"},   // 重复五月(后) —— 正常名, 非 伍月
        {240, 1,  5, "腊月"},   // CE 连续十二月: 前一个
        {240, 1, 20, "拾贰月"}, // 后一个 → 拾贰月 (保持)
        { 24, 1, 10, "拾贰月"}, // 跨年延续的 拾贰月
    };
    for (const auto &c : cases) {
        std::string mn = monthName(c.y, c.m, c.d);
        std::printf("  %d-%02d-%02d -> %s (期望 %s)\n", c.y, c.m, c.d, mn.c_str(), c.mon);
        check(mn == c.mon, std::string("改历月名 ") + std::to_string(c.y) + "-"
              + std::to_string(c.m) + "-" + std::to_string(c.d) + " 应为 " + c.mon);
    }
}

} // namespace

int main() {
    std::printf("=== 阴历月名验证测试 ===\n");

    testModern();
    // 春秋/战国: 年终置闰月称"十三月"
    testAncientRange(-600, -585, "十三月", "后九月");
    // 秦汉: 年终置闰月称"后九月"
    testAncientRange(-200, -185, "后九月", "十三月");
    testYearParse();
    testShengXiao();
    testKnownDay();
    testYearCalendar();
    testRepeatedMonths();
    testJsOracle();
    testReformMonthsFixed();

    std::printf("\n=== 结果: PASS=%d FAIL=%d ===\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
