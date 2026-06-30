// ═══════════════════════════════════════════════════════════════════
//  老黄历查表实现
// ═══════════════════════════════════════════════════════════════════
#include "almanac.h"
#include "almanac_data.h"
#include "almanac_dong_gong.h"

#include <cstdint>
#include <string>
#include <string_view>

namespace sxwnl::almanac {

namespace {

using namespace data;

// ── 容错索引访问: 越界则返回空串/默认值 ───────────────────
template <typename T, std::size_t N>
inline std::string_view safe_at(const std::array<T, N>& arr, int idx) {
    return (idx >= 0 && idx < static_cast<int>(N)) ? std::string_view(arr[idx])
                                                   : std::string_view{};
}

inline std::string sv2str(std::string_view sv) { return std::string(sv); }

// ── 1) 二十八宿 ────────────────────────────────────────────
void fill_xiu(const DayContext& c, DayAlmanac& out) {
    if (c.day_zhi < 0 || c.day_zhi > 11) return;
    if (c.week_day < 0 || c.week_day > 6) return;

    const int idx = kXiuByDayZhiWeekday[c.day_zhi % 4][c.week_day];
    const XiuAttr& a = kXiu28[idx];
    out.xiu        = sv2str(a.name);
    out.xiu_zheng  = sv2str(a.zheng);
    out.xiu_animal = sv2str(a.animal);
    out.xiu_luck   = sv2str(a.luck);
    out.xiu_gong   = sv2str(a.gong);
}

// ── 2) 黄道黑道 12 神 ──────────────────────────────────────
//   起例:  寅月  子日为青龙;
//           卯月  寅日为青龙;
//           ...
//           即 青龙日支 = ((month_zhi - 2) * 2 + 12) % 12 .
//   再以 day_zhi 与青龙日支的距离, 在 12 神序列中取出当日值神.
void fill_twelve_god(const DayContext& c, DayAlmanac& out) {
    if (c.month_zhi < 0 || c.month_zhi > 11) return;
    if (c.day_zhi   < 0 || c.day_zhi   > 11) return;

    const int qing_long_zhi = ((c.month_zhi - 2) * 2 + 24) % 12;  // +24 防负数
    const int god_idx       = ((c.day_zhi - qing_long_zhi) + 12) % 12;
    const TwelveGod& g = kTwelveGod[god_idx];
    out.twelve_god    = sv2str(g.name);
    out.is_huang_dao  = g.is_huang_dao;
    out.huang_hei     = g.is_huang_dao ? "黄道" : "黑道";
}

// ── 3) 冲煞 ────────────────────────────────────────────────
void fill_chong_sha(const DayContext& c, DayAlmanac& out) {
    if (c.day_gan < 0 || c.day_gan > 9)  return;
    if (c.day_zhi < 0 || c.day_zhi > 11) return;

    const ChongShaInfo& cs = kChongSha[c.day_zhi];
    out.chong_sheng_xiao = sv2str(cs.sheng_xiao);
    out.sha              = sv2str(cs.sha);

    // 冲日干支: 冲日干 = (day_gan + 4) % 10; 冲日支 = (day_zhi + 6) % 12
    const int chong_gan = (c.day_gan + 4) % 10;
    const int chong_zhi = (c.day_zhi + 6) % 12;
    out.chong_gan_zhi = sv2str(kGanName[chong_gan]) + sv2str(kZhiName[chong_zhi]);
}

// ── 4) 五吉神方位 ──────────────────────────────────────────
void fill_positions(const DayContext& c, DayAlmanac& out) {
    if (c.day_gan < 0 || c.day_gan > 9) return;
    out.xi_shen_fang   = sv2str(kPosXiShen     [c.day_gan]);
    out.yang_gui_fang  = sv2str(kPosYangGuiShen[c.day_gan]);
    out.yin_gui_fang   = sv2str(kPosYinGuiShen [c.day_gan]);
    out.fu_shen_fang   = sv2str(kPosFuShen     [c.day_gan]);
    out.cai_shen_fang  = sv2str(kPosCaiShen    [c.day_gan]);
}

// ── 5) 彭祖百忌 ────────────────────────────────────────────
void fill_peng_zu(const DayContext& c, DayAlmanac& out) {
    out.peng_zu_gan = sv2str(safe_at(kPengZuGan, c.day_gan));
    out.peng_zu_zhi = sv2str(safe_at(kPengZuZhi, c.day_zhi));
}

// ── 6) 择日典籍语录 (董公择日要诀, 通用接口可扩展) ─────────
//
// luck 启发式: 同一段内可能既述"甲子大吉"又述"壬子大凶", 故采用
//   "吉/凶/混/平/空" 五档. 规则:
//     含"大吉"且不含"大凶/凶" -> "吉"
//     含"大凶"或多次"凶"且不含"大吉" -> "凶"
//     两者皆有                  -> "混"
//     仅含"次吉/小吉"           -> "平"
//     都没匹配                  -> ""
inline bool contains(std::string_view s, std::string_view kw) {
    return s.find(kw) != std::string_view::npos;
}

inline std::string infer_luck(std::string_view body) {
    const bool good = contains(body, "大吉") || contains(body, "上吉") ||
                      contains(body, "全吉") || contains(body, "顺利") ||
                      contains(body, "顺遂") || contains(body, "诸事称心");
    const bool bad  = contains(body, "大凶") || contains(body, "更凶") ||
                      contains(body, "百事不宜") || contains(body, "诸事不宜") ||
                      contains(body, "百事皆忌") || contains(body, "百事忌") ||
                      contains(body, "百事不利") || contains(body, "主凶") ||
                      contains(body, "凶不可言");
    if (good && bad) return "混";
    if (good)        return "吉";
    if (bad)         return "凶";
    if (contains(body, "次吉") || contains(body, "小吉") || contains(body, "可用"))
        return "平";
    return {};
}

void fill_quotes(const DayContext& c, DayAlmanac& out) {
    if (c.month_zhi < 0 || c.month_zhi > 11) return;
    if (c.day_zhi   < 0 || c.day_zhi   > 11) return;

    // 董公择日要诀 (按 [月支][日支] 查表; 建除星由 (day_zhi - month_zhi) 推得)
    std::string_view body = kDongGongQuotes[c.month_zhi][c.day_zhi];
    if (!body.empty()) {
        const int  jian_idx = ((c.day_zhi - c.month_zhi) + 12) % 12;
        AlmanacQuote q;
        q.source = "董公择日要诀";
        q.title  = sv2str(safe_at(kDongGongMonthName, c.month_zhi))
                 + "·" + sv2str(safe_at(kJian12NameDongGong, jian_idx))
                 + sv2str(safe_at(kZhiName, c.day_zhi)) + "日";
        q.body   = sv2str(body);
        q.luck   = infer_luck(body);
        out.quotes.push_back(std::move(q));
    }
    // 未来扩展点: 玉匣记 / 通胜经 / 协纪辨方 ...
}

} // anonymous namespace

// 声明: 实现在 almanac_shensha.cpp
void fill_shen_sha(const DayContext& ctx, DayAlmanac& out);

DayAlmanac query(const DayContext& ctx) {
    DayAlmanac r;
    fill_xiu        (ctx, r);
    fill_twelve_god (ctx, r);
    fill_chong_sha  (ctx, r);
    fill_positions  (ctx, r);
    fill_peng_zu    (ctx, r);
    fill_quotes     (ctx, r);
    fill_shen_sha   (ctx, r);  // 必须在 fill_xiu 之后(金神七煞依赖 out.xiu)
    return r;
}

} // namespace sxwnl::almanac
