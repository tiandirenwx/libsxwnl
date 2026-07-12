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

#include <cstdio>
#include <map>
#include <string>
#include <vector>

#include "sxwnl_capi.h"

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

    std::printf("\n=== 结果: PASS=%d FAIL=%d ===\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
