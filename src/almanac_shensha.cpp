// ═══════════════════════════════════════════════════════════════════
//  神煞计算 + 宜忌推导 + 吉时 + 用事择吉
//
//  设计:
//    所有规则都基于 DayContext 中由 libsxwnl 提供的"日柱/月支/日支/
//    农历月日"等已计算字段, 本文件零额外日历运算.
//
//  分模块:
//    [A] compute_shen_sha    — 35+ 神煞投票
//    [B] derive_yi_ji        — 神煞 → 宜忌字典聚合
//    [C] compute_lucky_hours — 5吉时 + 贵人时
//    [D] compute_events      — 6 类用事择吉
//    [E] compute_notes       — 节气/月份特别提示
// ═══════════════════════════════════════════════════════════════════
#include "almanac.h"
#include "almanac_data.h"
#include "almanac_dong_gong.h"
#include "almanac_shensha.h"

#include <algorithm>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace sxwnl::almanac {

namespace {

using namespace data;

// ─── 辅助 ──────────────────────────────────────────────────

inline int gz60_of_day(const DayContext& c) {
    return gz60_idx(c.day_gan, c.day_zhi);
}

inline void push_shensha(std::vector<ShenSha>& v, std::string n, bool lucky, int w) {
    for (auto& s : v) if (s.name == n) return;   // 去重
    v.push_back({std::move(n), lucky, w});
}

inline bool has_shensha(const std::vector<ShenSha>& v, std::string_view n) {
    for (auto& s : v) if (s.name == n) return true;
    return false;
}

// ─── [A] 35+ 神煞投票 ─────────────────────────────────────

void compute_shen_sha(const DayContext& c, std::vector<ShenSha>& out) {
    const int gz   = gz60_of_day(c);
    const int sea  = season_of(c.month_zhi);
    const int grp  = month_group_of(c.month_zhi);
    const int mz   = c.month_zhi;
    const int dz   = c.day_zhi;
    const int dg   = c.day_gan;

    // —— 固定日柱集合 ——
    if (kJinShen.has(gz))           push_shensha(out, "进神",            true,  2);
    if (kShiLing.has(gz))           push_shensha(out, "十灵日",          true,  2);
    if (kXuanNuSanQi.has(gz))       push_shensha(out, "玄女三奇",        true,  3);
    if (kShiQuanFuGui.has(gz))      push_shensha(out, "十全富贵",        true,  3);
    if (kDaZang.has(gz))            push_shensha(out, "大葬日",          true,  1);
    if (kMingFei.has(gz))           push_shensha(out, "鸣吠日",          true,  1);
    if (kXiaoZang.has(gz))          push_shensha(out, "小葬日",          true,  1);
    if (kMingFeiDui.has(gz))        push_shensha(out, "鸣吠对日",        true,  1);
    if (kWuHe.has(gz))              push_shensha(out, "五合吉日",        true,  2);

    if (kWuLi.has(gz))              push_shensha(out, "五离日",          false, 2);
    if (kBaiHuRuZhongGong.has(gz))  push_shensha(out, "白虎入中宫",      false, 3);
    if (kJiuTuGui.has(gz))          push_shensha(out, "九土鬼",          false, 2);

    // —— 季节 × 日柱 ——
    if (kTianShe.has(sea, gz))         push_shensha(out, "天赦日",      true,  3);
    if (kTianDiZhuanSha.has(sea, gz))  push_shensha(out, "天地转煞",    false, 3);
    if (kZhengSiFei.has(sea, gz))      push_shensha(out, "正四废",      false, 3);

    // —— 月组 × 日柱 ——
    if (kShaGong.has(grp, gz))       push_shensha(out, "煞贡",          true,  2);
    if (kZhiXing.has(grp, gz))       push_shensha(out, "直星",          true,  2);
    if (kRenZhuan.has(grp, gz))      push_shensha(out, "人专",          true,  2);
    if (kHuoXing.has(grp, gz))       push_shensha(out, "火星日",        false, 1);

    // —— 月支 → 日支 ——
    if (kHuangSha[mz]    == dz)      push_shensha(out, "黄沙日",        false, 1);
    if (kXiaoHongSha[mz] == dz)      push_shensha(out, "小红沙",        false, 2);
    if (kTianZei[mz]     == dz)      push_shensha(out, "天贼日",        false, 1);
    if (kYueYan[mz]      == dz)      push_shensha(out, "月厌大祸",      false, 3);
    if (kBingXiao[mz]    == dz)      push_shensha(out, "冰消瓦解",      false, 2);
    if (kShouSi[mz]      == dz)      push_shensha(out, "受死日",        false, 3);
    if (kWangWang[mz]    == dz)      push_shensha(out, "往亡日",        false, 2);
    if (kSanHe[mz]       == dz)      push_shensha(out, "三合吉日",      true,  3);
    if (kLiuHe[mz]       == dz)      push_shensha(out, "六合吉日",      true,  2);
    if (kYueDe[mz]       == dg)      push_shensha(out, "月德",          true,  3);

    // —— 日干 → 日支 ——
    if (kTianYiGuiRen[dg].a == dz || kTianYiGuiRen[dg].b == dz)
        push_shensha(out, "天乙贵人", true, 3);
    if (kWenChang[dg] == dz)         push_shensha(out, "文昌",          true,  2);

    // —— 三奇 (日干属甲戊庚/乙丙丁/壬癸辛) ——
    //   传统三奇要看年/月/日柱齐全, 此处仅单日提示
    if (dg == 0 || dg == 4 || dg == 6) push_shensha(out, "天上三奇(甲戊庚)", true, 1);
    if (dg == 1 || dg == 2 || dg == 3) push_shensha(out, "地下三奇(乙丙丁)", true, 1);
    if (dg == 7 || dg == 8 || dg == 9) push_shensha(out, "人中三奇(壬癸辛)", true, 1);

    // —— 重日 (日支 巳/亥) ——
    if (dz == 5 || dz == 11)         push_shensha(out, "重日",          false, 1);

    // —— 复日 (农历月 → 日干) ——
    if (c.lunar_month >= 1 && c.lunar_month <= 12 &&
        kFuRi[c.lunar_month] == dg)  push_shensha(out, "复日",          false, 1);

    // —— 绝烟火 (按月份-1 mod 4) ——
    if (c.lunar_month >= 1 && c.lunar_month <= 12) {
        const int q = (c.lunar_month - 1) % 4;
        if (kJueYanHuo.has(q, gz))   push_shensha(out, "绝烟火",        false, 2);
    }

    // —— 杨公忌 / 赤口 (农历日) ——
    if (is_yang_gong_ji(c.lunar_month, c.lunar_day))
        push_shensha(out, "杨公忌日", false, 3);
    if (is_chi_kou(c.lunar_month, c.lunar_day))
        push_shensha(out, "赤口日",   false, 1);

    // —— 十恶大败 (年干组 × 农历月 × 日柱) ——
    const int ygg = year_gan_group_of(c.year_gan);
    for (auto& it : kShiE) {
        if (it.year_gan_group == ygg && it.lunar_month == c.lunar_month && it.gz60 == gz) {
            push_shensha(out, "十恶大败", false, 3);
            break;
        }
    }

    // —— 四绝日 / 四离日 (节气前一日) ——
    //   Jqmc 索引: 0冬至 3立春 6春分 9立春 12夏至 15立秋 18秋分 21立冬
    if (c.tomorrow_jie_qi == 3 || c.tomorrow_jie_qi == 9 ||
        c.tomorrow_jie_qi == 15 || c.tomorrow_jie_qi == 21) {
        push_shensha(out, "四绝日", false, 3);
    }
    if (c.tomorrow_jie_qi == 0 || c.tomorrow_jie_qi == 6 ||
        c.tomorrow_jie_qi == 12 || c.tomorrow_jie_qi == 18) {
        push_shensha(out, "四离日", false, 2);
    }

    // —— 金神七煞 (二十八宿 角亢奎娄牛鬼星, 通过当前日宿名判定) ——
    //   依赖 fill_xiu 已写入 out.xiu, 此处用枚举名直接匹配
    //   (注: 此函数运行前 out.xiu 未必有值, 因此放到 query 序列后段补判)

    // —— 按权重排序: 权重大的优先 ——
    std::stable_sort(out.begin(), out.end(),
        [](const ShenSha& a, const ShenSha& b) { return a.weight > b.weight; });
}

// 二十八宿驱动的"金神七煞"补判: 在 fill_xiu 之后运行
void compute_jin_shen_qi_sha(const DayAlmanac& a, std::vector<ShenSha>& out) {
    static const std::set<std::string_view> kJinShenXiu = {
        "角", "亢", "奎", "娄", "牛", "鬼", "星"
    };
    if (kJinShenXiu.count(a.xiu)) {
        push_shensha(out, "金神七煞", false, 3);
    }
}

// ─── [B] 神煞 → 宜忌字典 ──────────────────────────────────

struct YiJiEffect {
    std::vector<std::string_view> yi;
    std::vector<std::string_view> ji;
};

inline const std::unordered_map<std::string, YiJiEffect>& effects() {
    static const std::unordered_map<std::string, YiJiEffect> m = {
        // 吉神
        {"天赦日",    {{"祈福","祭祀","求嗣","解冤","赦罪","动土","修造"}, {}}},
        {"月德",      {{"祭祀","祈福","出行","嫁娶","入宅","求嗣"},        {}}},
        {"天乙贵人",  {{"上官赴任","求贵","见贵","求财","订盟"},            {}}},
        {"文昌",      {{"入学","求嗣","祭祀","求文"},                        {}}},
        {"玄女三奇",  {{"修方","偷修","破土","入宅"},                        {}}},
        {"煞贡",      {{"修造","动土","开张","出行"},                        {}}},
        {"直星",      {{"修造","嫁娶","入宅"},                                {}}},
        {"人专",      {{"嫁娶","订盟","纳采"},                                {}}},
        {"进神",      {{"修造","祭祀","求财"},                                {}}},
        {"十灵日",    {{"祈福","求贵","出行"},                                {}}},
        {"十全富贵",  {{"造葬","嫁娶","上官","开市"},                        {}}},
        {"五合吉日",  {{"嫁娶","订盟","纳采"},                                {}}},
        {"三合吉日",  {{"嫁娶","开市","结盟","出行","入宅"},                {}}},
        {"六合吉日",  {{"嫁娶","订盟","会亲"},                                {}}},
        {"大葬日",    {{"安葬","破土"},                                       {}}},
        {"鸣吠日",    {{"安葬"},                                              {}}},
        {"小葬日",    {{"安葬"},                                              {}}},
        {"鸣吠对日",  {{"安葬"},                                              {}}},
        // 凶神
        {"白虎入中宫",{{}, {"嫁娶","入宅","起造","动土","修造"}}},
        {"金神七煞",  {{}, {"出兵","行船","起造","婚姻"}}},
        {"天地转煞",  {{}, {"动土","破土","修造","葬埋"}}},
        {"正四废",    {{}, {"嫁娶","开市","上官","出行","入宅"}}},
        {"月厌大祸",  {{}, {"嫁娶","出行","订盟","纳采"}}},
        {"受死日",    {{}, {"嫁娶","出行","求医","开市","入宅"}}},
        {"往亡日",    {{}, {"赴任","出行","嫁娶","求谋"}}},
        {"冰消瓦解",  {{}, {"修造","嫁娶","入宅","开市"}}},
        {"小红沙",    {{}, {"嫁娶","出行","入宅"}}},
        {"黄沙日",    {{}, {"出行"}}},
        {"天贼日",    {{}, {"竖柱","迁宅","动土","开仓库"}}},
        {"九土鬼",    {{}, {"上任","出行","起造","动土","交易"}}},
        {"火星日",    {{}, {"修造","上盖","起灶","裁衣"}}},
        {"绝烟火",    {{}, {"造葬","移徙","进宅","嫁娶"}}},
        {"杨公忌日",  {{}, {"嫁娶","开市","入宅","出行","动土"}}},
        {"赤口日",    {{}, {"会亲","订盟","出行"}}},
        {"五离日",    {{}, {"嫁娶","订盟"}}},
        {"十恶大败",  {{}, {"嫁娶","开市","上官","入宅","出行","动土"}}},
        {"复日",      {{}, {"安葬","破土"}}},
        {"重日",      {{}, {"安葬"}}},
        {"四绝日",    {{}, {"开市","上官","出行","入宅","嫁娶","动土"}}},
        {"四离日",    {{}, {"嫁娶","出行","会客"}}},
    };
    return m;
}

void derive_yi_ji(const std::vector<ShenSha>& ss,
                  std::vector<std::string>& yi,
                  std::vector<std::string>& ji) {
    // (event_name → score). 吉神 +weight, 凶神 -weight.
    std::unordered_map<std::string, int> score;
    const auto& dict = effects();

    for (const auto& s : ss) {
        auto it = dict.find(s.name);
        if (it == dict.end()) continue;
        const int delta = s.weight;
        if (s.is_lucky) {
            for (auto e : it->second.yi) score[std::string(e)] += delta;
        } else {
            for (auto e : it->second.ji) score[std::string(e)] -= delta;
        }
    }

    // 拆成 yi / ji 两个列表 (score>0=宜, score<0=忌, 同名优先看冲突最大者)
    std::vector<std::pair<std::string,int>> sorted(score.begin(), score.end());
    std::stable_sort(sorted.begin(), sorted.end(),
        [](auto& a, auto& b) {
            if (std::abs(a.second) != std::abs(b.second))
                return std::abs(a.second) > std::abs(b.second);
            return a.first < b.first;
        });

    constexpr int kMaxEach = 5;
    for (auto& [name, sc] : sorted) {
        if (sc > 0 && (int)yi.size() < kMaxEach) yi.push_back(name);
        else if (sc < 0 && (int)ji.size() < kMaxEach) ji.push_back(name);
    }
}

// ─── [C] 吉曜时法 ─────────────────────────────────────────
//   按日支查 5 个吉时: 福德/凤辇/宝光/太乙/少微
//   贵人时按日干算 (天乙贵人在 时支 = 贵人方位之一)
//
//   规律: 日支每+1, 5吉时全部 +2 (mod 12). 用基准日支 子(0) 推所有.
//
//   基准 (子日): 福德=子, 凤辇=午, 宝光=丑, 太乙=申, 少微=酉
//
//   注: 文档表头列了 6 列(福德/凤辇/宝光/太乙/少微/贵人), 但表里只有 5
//        个数值. "贵人"是按日干独立算的, 不在按日支的表里.
void compute_lucky_hours(const DayContext& c, std::vector<LuckyHour>& out) {
    if (c.day_zhi < 0 || c.day_zhi > 11) return;
    constexpr int base[5] = { 0, 6, 1, 8, 9 };  // 子日 福/凤/宝/太/少
    constexpr const char* names[5] = { "福德", "凤辇", "宝光", "太乙", "少微" };
    const int shift = (c.day_zhi * 2) % 12;
    for (int i = 0; i < 5; ++i)
        out.push_back({ names[i], (base[i] + shift) % 12 });

    // 贵人时: 天乙贵人方 (按日干)
    if (c.day_gan >= 0 && c.day_gan < 10) {
        const auto p = kTianYiGuiRen[c.day_gan];
        out.push_back({ "贵人(阳)", p.a });
        out.push_back({ "贵人(阴)", p.b });
    }
}

// ─── [D] 6 类用事择吉 ─────────────────────────────────────
//
//   原理: 综合命中神煞与典籍专项吉日清单, 给出 suitable=true/false
//         + reason 一句话.
//
//   typedef: void check_event(name, suitable_gz_set, ban_shensha, advice)

namespace event_data {
    // 安床吉日(按农历月) - 直接从董公诀文档表
    //   每月最多 12 个日柱, 用 GanZhiSet60 表示
    inline constexpr std::array<GanZhiSet60, 13> kAnChuang = {{
        {},                                       // [0] unused
        make_gz_set({{1,3},{3,9},{9,9},{1,1},{3,1},{9,1},{5,3},{7,3}}),                  // 正月
        make_gz_set({{0,2},{2,2},{8,2},{1,7},{3,7},{5,7},{1,11},{3,11}}),                // 二月
        make_gz_set({{0,0},{2,0},{6,0},{8,0},{9,3},{1,3},{5,3},{7,3},{1,5},{3,5}}),      // 三月
        make_gz_set({{0,0},{2,0},{6,0},{3,1},{1,3},{7,3},{6,6},{0,4},{2,4},{2,6}}),      // 四月
        make_gz_set({{0,2},{2,2},{6,2},{8,2},{8,4},{0,4},{2,4}}),                        // 五月
        make_gz_set({{0,2},{2,2},{8,2},{1,11},{3,11}}),                                  // 六月
        make_gz_set({{0,0},{2,0},{8,0},{6,0},{2,4},{6,4}}),                              // 七月
        make_gz_set({{1,1},{3,1},{7,1},{9,1},{5,1},{0,4},{2,4},{8,4},{6,4}}),            // 八月
        make_gz_set({{0,6},{6,6},{2,6},{7,3},{7,9},{3,9},{9,9},{1,11},{3,11}}),          // 九月
        make_gz_set({{0,0},{2,0},{6,0},{1,7},{7,7}}),                                    // 十月
        make_gz_set({{0,2},{2,2},{6,2},{1,5},{5,5},{7,7},{1,11},{3,11}}),                // 冬月
        make_gz_set({{0,0},{2,0},{8,0},{0,2},{2,2}}),                                    // 腊月
    }};

    // 动土吉日 (固定日柱集 - 文档"动土吉日"段)
    inline constexpr GanZhiSet60 kDongTu = make_gz_set({
        {5,5},{0,10},{1,11},{4,2},{5,3},{8,6},{0,8},{1,9},{4,0},{6,2},
        {1,7},{5,11},{8,2},{9,3},{2,6},{4,8},{5,9},{8,0},{5,7},{6,8},{7,9}
    });

    // 竖柱吉日
    inline constexpr GanZhiSet60 kShuZhu = make_gz_set({
        {1,5},{7,1},{0,2},{1,11},{1,9},{5,9},{8,0},{5,7},{6,8},{4,0},
        {1,7},{5,11},{5,3},{0,8},{5,1},{6,2},{9,3},{4,8},{8,10}
    });

    // 上梁吉日 (较多, 36 个)
    inline constexpr GanZhiSet60 kShangLiang = make_gz_set({
        {0,0},{1,1},{3,3},{4,4},{5,5},{6,6},{7,7},{8,8},{0,10},{2,0},
        {4,2},{6,4},{8,6},{0,8},{2,10},{4,0},{6,2},{0,6},{2,8},{3,9},
        {4,10},{5,11},{6,0},{7,1},{8,2},{9,3},{1,5},{3,7},{5,9},{7,11},
        {9,1},{1,3},{3,5},{5,7},{7,9},{9,11}
    });

    // 建屋吉日 (较多, 34 个)
    inline constexpr GanZhiSet60 kJianWu = make_gz_set({
        {0,0},{3,3},{4,4},{5,5},{7,7},{8,8},{9,9},{2,0},{3,1},{5,3},
        {6,4},{9,7},{0,8},{5,9},{2,10},{4,0},{6,2},{3,9},{9,5},{1,7},
        {5,11},{7,1},{8,2},{9,3},{0,4},{1,5},{4,8},{6,10},{7,11},{9,1},
        {1,3},{2,4},{6,8},{7,9}
    });
}  // namespace event_data

void compute_events(const DayContext& c,
                    const std::vector<ShenSha>& ss,
                    std::vector<EventAdvice>& out) {
    const int gz = gz60_of_day(c);

    auto reason_from = [&](std::string event, bool yi) -> std::string {
        // 找一条最高权重支持/反对该事件的神煞作 reason
        const auto& dict = effects();
        for (const auto& s : ss) {
            auto it = dict.find(s.name);
            if (it == dict.end()) continue;
            const auto& list = yi ? it->second.yi : it->second.ji;
            for (auto& e : list) {
                if (e == event) return s.name;
            }
        }
        return yi ? "本日干支契合" : "本日干支不利";
    };

    auto add = [&](const std::string& name, bool suitable, std::string reason) {
        out.push_back({name, suitable, std::move(reason)});
    };

    // 动土
    {
        bool ji = has_shensha(ss, "动土") /* placeholder */
            || has_shensha(ss, "天地转煞") || has_shensha(ss, "天贼日")
            || has_shensha(ss, "杨公忌日") || has_shensha(ss, "白虎入中宫");
        bool yi = event_data::kDongTu.has(gz) || has_shensha(ss, "天赦日") || has_shensha(ss, "月德");
        if (yi && !ji)       add("动土", true,  reason_from("动土", true));
        else if (ji)         add("动土", false, reason_from("动土", false));
    }
    // 竖柱
    {
        bool yi = event_data::kShuZhu.has(gz) || has_shensha(ss, "月德");
        bool ji = has_shensha(ss, "天地转煞") || has_shensha(ss, "杨公忌日");
        if (yi && !ji)       add("竖柱", true,  reason_from("竖柱", true));
        else if (ji)         add("竖柱", false, reason_from("竖柱", false));
    }
    // 上梁
    {
        bool yi = event_data::kShangLiang.has(gz);
        bool ji = has_shensha(ss, "杨公忌日") || has_shensha(ss, "白虎入中宫");
        if (yi && !ji)       add("上梁", true,  reason_from("上梁", true));
        else if (ji)         add("上梁", false, reason_from("上梁", false));
    }
    // 建屋
    {
        bool yi = event_data::kJianWu.has(gz) || has_shensha(ss, "天赦日");
        bool ji = has_shensha(ss, "白虎入中宫") || has_shensha(ss, "天地转煞");
        if (yi && !ji)       add("建屋", true,  reason_from("建屋", true));
        else if (ji)         add("建屋", false, reason_from("建屋", false));
    }
    // 安灶 (建除规则: 开/定/平/收 大吉; 建/破 妨家长)
    //   建除星索引 (day_zhi - month_zhi + 12) % 12
    //   建=0 除=1 满=2 平=3 定=4 执=5 破=6 危=7 成=8 收=9 开=10 闭=11
    {
        const int jian = ((c.day_zhi - c.month_zhi) + 12) % 12;
        const bool yi = (jian == 10 || jian == 4 || jian == 3 || jian == 9);
        const bool ji = (jian == 0 || jian == 6 || c.day_gan == 2 /* 丙 */);
        if (yi && !ji) add("安灶", true,  "建除值开/定/平/收");
        else if (ji)   add("安灶", false, "建破丙日不宜安灶");
    }
    // 安床
    if (c.lunar_month >= 1 && c.lunar_month <= 12) {
        const bool yi = event_data::kAnChuang[c.lunar_month].has(gz);
        const bool ji = has_shensha(ss, "杨公忌日") || has_shensha(ss, "受死日")
                     || has_shensha(ss, "正四废");
        if (yi && !ji) add("安床", true,  std::string(kDongGongMonthName[c.month_zhi]) + "安床吉日");
        else if (ji)   add("安床", false, reason_from("出行", false));  // 借用"出行"凶神
    }
    // 嫁娶 / 出行 / 安葬 已经通过 yi/ji 字段表达, 此处不重复
}

// ─── [E] 节气/月份特别提示 ───────────────────────────────

void compute_notes(const DayContext& c, std::vector<std::string>& notes) {
    // 三煞方位提示 (每月不同)
    //   寅卯辰(春)三煞在北; 巳午未(夏)在西; 申酉戌(秋)在南; 亥子丑(冬)在东
    //   (与官方董公诀对应的"三煞方"略)
    const int sea = season_of(c.month_zhi);
    constexpr const char* san_sha_dir[4] = { "北方亥子丑", "西方申酉戌", "南方巳午未", "东方寅卯辰" };
    notes.push_back(std::string("三煞在") + san_sha_dir[sea] + "方, 忌修造动土");
}

}  // anonymous namespace

// ─── 对外入口 (在 almanac.cpp 的 query() 中调用) ─────────

void fill_shen_sha(const DayContext& c, DayAlmanac& out) {
    compute_shen_sha(c, out.shen_sha);
    compute_jin_shen_qi_sha(out, out.shen_sha);
    derive_yi_ji(out.shen_sha, out.yi, out.ji);
    compute_lucky_hours(c, out.lucky_hours);
    compute_events(c, out.shen_sha, out.event_advices);
    compute_notes(c, out.notes);
}

}  // namespace sxwnl::almanac
