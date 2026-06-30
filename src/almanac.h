#pragma once

// ═══════════════════════════════════════════════════════════════════
//  老黄历模块 (Almanac)
//
//  架构契约:
//    本模块是"无算法的查表器" -- 不重算干支/农历/节气, 全部由 libsxwnl
//    在 DayContext 中提供. 这样保证全项目"干支/农历"只有一个真源,
//    永不出现日历显示甲子, 老黄历显示癸亥的边界 bug.
//
//  数据来源:
//    查表数据移植自 lunar-javascript (https://github.com/6tail/lunar-javascript,
//    MIT). 原项目用模板字符串组装, 本模块直接转为中文字面量以减少运行时开销.
//
//  扩展方式:
//    新增黄历项 -> 在 DayAlmanac 增字段 -> almanac_data.h 加数据表
//    -> almanac.cpp 加 fillXxx() 函数 -> query() 中调用.
// ═══════════════════════════════════════════════════════════════════

#include <string>
#include <vector>

namespace sxwnl::almanac {

// ── 输入: 完全由 libsxwnl 算好 ─────────────────────────────
struct DayContext {
    int    year_gan   = 0;   // 0..9
    int    year_zhi   = 0;   // 0..11
    int    month_gan  = 0;
    int    month_zhi  = 0;   // 按节气分月 (libsxwnl 已处理)
    int    day_gan    = 0;
    int    day_zhi    = 0;
    int    week_day   = 0;   // 0=Sunday, 1=Monday, ... 6=Saturday
    int    lunar_month = 1;  // 1..12
    int    lunar_day   = 1;  // 1..30
    bool   is_leap_month = false;
    double julian_day  = 0;  // 12:00 对应的儒略日 (含 J2000 偏移)

    // 24 节气索引: 0=春分, ... 与 sxwnl Jqmc 一致.
    //   today_jie_qi:    当日是否为节气, -1 表示否
    //   tomorrow_jie_qi: 明日是否为节气, -1 表示否
    //   用于检测四绝(立春/夏/秋/冬前一日) 与 四离(春分/秋分/夏至/冬至前一日)
    int    today_jie_qi    = -1;
    int    tomorrow_jie_qi = -1;
};

// ── 一条择日典籍语录 (董公择日要诀 / 玉匣记 / 通胜经 ...) ──
//   通用结构: 数据由 almanac_data_*.h 提供, 业务方按 source 区分.
struct AlmanacQuote {
    std::string source;   // 典籍来源 "董公择日要诀"
    std::string title;    // 段标题   "正月·开子日"
    std::string body;     // 原文     "甲子自死之金... 用之大吉"
    std::string luck;     // 整体基调 "吉"/"凶"/"平"/"混"/""  (启发式自动推断)
};

// ── 一个神煞 (天德/月厌/白虎入中宫/三合 ...) ──────────────
struct ShenSha {
    std::string name;     // "天德" "月厌大祸" ...
    bool is_lucky;        // 吉神=true / 凶神=false
    int  weight;          // 权重 (大煞=3, 中煞=2, 一般=1)
};

// ── 一个吉时 (按日支查的 5 吉时 + 按日干查的贵人时) ────────
struct LuckyHour {
    std::string name;     // "福德"/"凤辇"/"宝光"/"太乙"/"少微"/"贵人"
    int  zhi;             // 0..11
};

// ── 用事择吉建议 (单条) ────────────────────────────────
struct EventAdvice {
    std::string event;    // "动土"/"安床"/"嫁娶" ...
    bool        suitable; // 是否适合 (综合命中规则后的结论)
    std::string reason;   // 一句话理由 (命中的关键神煞)
};

// ── 输出: 一日的老黄历完整数据 ─────────────────────────────
struct DayAlmanac {
    // ── 二十八宿 ──
    std::string xiu;            // "角"
    std::string xiu_zheng;      // "木"  (七政: 木金土日月火水)
    std::string xiu_animal;     // "蛟"
    std::string xiu_luck;       // "吉" / "凶"
    std::string xiu_gong;       // "东"  (四宫: 东青龙/南朱雀/西白虎/北玄武)

    // ── 黄道黑道 ──
    std::string twelve_god;     // "青龙" / "白虎" 等
    std::string huang_hei;      // "黄道" / "黑道"
    bool        is_huang_dao = false;

    // ── 冲煞 ──
    std::string chong_sheng_xiao; // "马"
    std::string chong_gan_zhi;    // "戊午"  (本日六十甲子中冲日的具体干支)
    std::string sha;              // "南"

    // ── 五吉神方位 ──
    std::string xi_shen_fang;     // 喜神   "东北"
    std::string yang_gui_fang;    // 阳贵神 "西南"
    std::string yin_gui_fang;     // 阴贵神 "东北"
    std::string fu_shen_fang;     // 福神   "东南"
    std::string cai_shen_fang;    // 财神   "东北"

    // ── 彭祖百忌 ──
    std::string peng_zu_gan;      // "甲不开仓 财物耗散"
    std::string peng_zu_zhi;      // "子不问卜 自惹祸殃"

    // ── 择日典籍语录 (董公择日要诀 / 玉匣记 / 通胜经 ... 可扩展) ──
    //   每条 quote 是一段独立典籍解读. 同一天不同典籍各自一条.
    std::vector<AlmanacQuote> quotes;

    // ── 神煞 (吉/凶, 按权重排序) ─────────────────────────
    std::vector<ShenSha> shen_sha;

    // ── 宜忌 (由 shen_sha 投票推得) ──────────────────────
    std::vector<std::string> yi;       // ["嫁娶","出行","安葬"]
    std::vector<std::string> ji;       // ["动土","上梁"]

    // ── 吉曜时法 (5 个按日支 + 1 个按日干, 共最多 6 条) ──
    std::vector<LuckyHour> lucky_hours;

    // ── 6 类用事择吉 (动土/竖柱/上梁/建屋/安灶/安床) ─────
    std::vector<EventAdvice> event_advices;

    // ── 特殊提示 (节气切换/转煞/三煞方位等) ──────────────
    //   每条都是与当日节气有关的"今日特别要注意的事".
    std::vector<std::string> notes;
};

DayAlmanac query(const DayContext& ctx);

} // namespace sxwnl::almanac
