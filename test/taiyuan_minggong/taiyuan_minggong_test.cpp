// ─────────────────────────────────────────────────────────────────────────────
// 胎元 / 命宫 —— 规格符合性测试
//
// 判据全部来自命理典籍与公开算例, 不从实现反抄:
//
//  胎元 (《三命通会》《渊海子平》):
//    "人生月后紧接着这个月的天干与生月后第三个月的地支相配, 就为胎元"
//    即 月柱天干进一位、地支进三位。
//
//  命宫 (《三命通会》·数字速查法):
//    寅一 卯二 辰三 巳四 午五 未六 申七 酉八 戌九 亥十 子十一 丑十二;
//    月数 + 时数 不满十四者以十四减之, 满十四者以廿六减之, 余数即命宫地支数;
//    命宫天干再由年干起五虎遁配出。
//    论月以"中气"为准 —— 出生已过本月中气者, 作次月计。
//
// ── A. 典籍算例逐条对照 ──────────────────────────────────────────────
//    A1 命宫: 巳月酉时→卯, 申月丑时→申, 辰月未时→午 (观星阁/算命网/苏民峰同源例)
//    A2 命宫干支完整例: 戊寅年 酉月 辰时, 白露后五日未过秋分 → 丙辰
//    A3 胎元: 月柱庚午→辛酉 等
//
// ── B. 纯函数全域性质 ────────────────────────────────────────────────
//    B1 胎元遍历 60 甲子: 干 +1、支 +3, 且胎元本身必为合法干支组合(干支同阴阳)
//    B2 命宫遍历 (年干10 × 月支12 × 时支12 × 过中气2) = 2880 组:
//       - 地支恒落在 [0,11], 天干恒落在 [0,9]
//       - 天干严格满足五虎遁: gan == ((年干%5)*2+2 + (命宫支-寅)) % 10
//       - 命宫干支必同阴阳(可构成六十甲子中的合法组合)
//    B3 命宫对"过中气"的响应: 过中气等价于月支进一位后未过中气
//    B4 非法入参一律返回空串, 不崩溃
//
// ── C. 与真实排盘联动 ────────────────────────────────────────────────
//    C1 中气必须严格落在 (月首节, 月末节) 开区间内 —— 三种历法 × 多年采样
//    C2 胎元与该盘月柱的关系恒成立 (端到端穿过 BaziBase)
//    C3 命宫与该盘 年干/月支/时支/过中气标志 的关系恒成立
//    C4 过中气标志的独立复核: 在节后逐时扫描, isAfterZhongQi() 必须单调翻转一次
//    C5 平气法下中气恰在两节正中 (平气按定义 24 等分, 节与中等距)
//
// 编译: make taiyuan_minggong_test
// 运行: ./build/bin/taiyuan_minggong_test [-v] [-range 1900 2100]
// ─────────────────────────────────────────────────────────────────────────────

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include "bazi.h"
#include "bazi_analysis.h"
#include "JD.h"
#include "sx_lang_zh.h"

namespace {

// ═══════════════════════════════════════════════════════════════════════
//  最小测试框架
// ═══════════════════════════════════════════════════════════════════════
bool gVerbose = false;

struct Group {
    const char *name;
    int pass = 0, fail = 0;
    std::string ctx;
    std::vector<std::string> msgs;

    void check(bool ok, const std::string &what) {
        if (ok) {
            pass++;
            if (gVerbose) std::printf("    [ok] %s %s\n", ctx.c_str(), what.c_str());
        } else {
            fail++;
            if (msgs.size() < 20) msgs.push_back(ctx + " | " + what);
        }
    }
};

std::vector<Group *> gGroups;

Group &group(const char *name) {
    static std::vector<Group *> owned;
    auto *g = new Group{name};
    owned.push_back(g);
    gGroups.push_back(g);
    return *g;
}

// 地支索引 → 命宫本数(寅一…丑十二)
int gongNum(int zhi) { return (zhi - 2 + 12) % 12 + 1; }

// 名称 → 索引
int ganIdx(const std::string &s) {
    for (int i = 0; i < 10; i++) if (s == Gan[i]) return i;
    return -1;
}
int zhiIdx(const std::string &s) {
    for (int i = 0; i < 12; i++) if (s == Zhi[i]) return i;
    return -1;
}

const char *lifaName(LiFaType lf) {
    switch (lf) {
    case XianDaiNongLifa_DingQiFa:    return "定气";
    case YuWuWeiZiPingLifa_DingDongZhi: return "平气定冬至";
    case YuWuWeiZiPingLifa_DingXiaZhi:  return "平气定夏至";
    }
    return "?";
}

// ═══════════════════════════════════════════════════════════════════════
//  A. 典籍算例
// ═══════════════════════════════════════════════════════════════════════
void testA() {
    Group &g = group("A 典籍算例");

    // A1 —— 只验地支, 年干任取(地支与年干无关)
    struct MgCase { const char *month, *hour, *want; const char *src; };
    const MgCase mg[] = {
        {"巳", "酉", "卯", "巳月酉时: 4+8=12 <14, 14-12=2 → 卯"},
        {"申", "丑", "申", "申月丑时: 7+12=19 >14, 26-19=7 → 申"},
        {"辰", "未", "午", "辰月未时: 3+6=9 <14, 14-9=5 → 午"},
        {"酉", "辰", "辰", "酉月辰时: 8+3=11 <14, 14-11=3 → 辰"},
    };
    for (const auto &c : mg) {
        int mz = zhiIdx(c.month), hz = zhiIdx(c.hour), want = zhiIdx(c.want);
        for (int yg = 0; yg < 10; yg++) {
            g.ctx = c.src;
            auto r = bazi::mingGong(yg, mz, hz, /*afterZhongQi=*/false);
            g.check(r[1] == want,
                    std::string("年干") + Gan[yg] + " 命宫支 得" +
                    (r[1] >= 0 ? Zhi[r[1]] : "?") + " 期望" + c.want);
        }
    }

    // A2 —— 完整干支例 (苏民峰/壹读 同源):
    //   戊寅年 酉月 辰时, 生于白露后五日, 未过秋分中气 → 命宫 丙辰
    g.ctx = "戊寅年 酉月 辰时 未过秋分";
    g.check(bazi::mingGongStr(ganIdx("戊"), zhiIdx("酉"), zhiIdx("辰"), false) == "丙辰",
            "得 " + bazi::mingGongStr(ganIdx("戊"), zhiIdx("酉"), zhiIdx("辰"), false) +
            " 期望 丙辰");

    //   丙年 命宫在辰 → 壬辰 ("丙辛之岁从庚起": 寅庚 卯辛 辰壬)
    //   构造: 月支酉(8) + 时支辰(3) = 11 → 14-11 = 3 → 辰
    g.ctx = "丙年 命宫辰";
    g.check(bazi::mingGongStr(ganIdx("丙"), zhiIdx("酉"), zhiIdx("辰"), false) == "壬辰",
            "得 " + bazi::mingGongStr(ganIdx("丙"), zhiIdx("酉"), zhiIdx("辰"), false) +
            " 期望 壬辰");

    //   壬年 命宫在午 → 丙午 ("丁壬壬寅顺行流": 寅壬 卯癸 辰甲 巳乙 午丙)
    //   构造: 月支辰(3) + 时支未(6) = 9 → 14-9 = 5 → 午
    g.ctx = "壬年 命宫午";
    g.check(bazi::mingGongStr(ganIdx("壬"), zhiIdx("辰"), zhiIdx("未"), false) == "丙午",
            "得 " + bazi::mingGongStr(ganIdx("壬"), zhiIdx("辰"), zhiIdx("未"), false) +
            " 期望 丙午");

    // A3 —— 胎元
    struct TyCase { const char *mg, *mz, *want; };
    const TyCase ty[] = {
        {"庚", "午", "辛酉"},  // 百度百科/苏民峰例
        {"丙", "寅", "丁巳"},
        {"癸", "亥", "甲寅"},  // 干支双回绕
        {"壬", "戌", "癸丑"},
    };
    for (const auto &c : ty) {
        g.ctx = std::string("胎元 月柱") + c.mg + c.mz;
        auto got = bazi::taiYuanStr(ganIdx(c.mg), zhiIdx(c.mz));
        g.check(got == c.want, "得 " + got + " 期望 " + c.want);
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  B. 纯函数全域性质
// ═══════════════════════════════════════════════════════════════════════
void testB() {
    Group &g = group("B 纯函数全域性质");

    // B1 胎元: 遍历 60 甲子
    for (int i = 0; i < 60; i++) {
        int mgan = i % 10, mzhi = i % 12;
        g.ctx = std::string("月柱") + Gan[mgan] + Zhi[mzhi];
        auto ty = bazi::taiYuan(mgan, mzhi);
        g.check(ty[0] == (mgan + 1) % 10, "干进一位");
        g.check(ty[1] == (mzhi + 3) % 12, "支进三位");
        // 干进1、支进3 后奇偶性同时翻转/保持 → 仍为合法干支(同阴阳)
        g.check(ty[0] % 2 == ty[1] % 2, "胎元干支同阴阳");
    }

    // B2 命宫: 全域 2880 组
    for (int yg = 0; yg < 10; yg++)
    for (int mz = 0; mz < 12; mz++)
    for (int hz = 0; hz < 12; hz++)
    for (int az = 0; az < 2; az++) {
        g.ctx = std::string("年干") + Gan[yg] + " 月支" + Zhi[mz] + " 时支" + Zhi[hz] +
                (az ? " 过中气" : " 未过中气");
        auto r = bazi::mingGong(yg, mz, hz, az != 0);
        g.check(r[1] >= 0 && r[1] < 12, "地支在域内");
        g.check(r[0] >= 0 && r[0] < 10, "天干在域内");
        if (r[1] < 0) continue;
        int wantGan = (((yg % 5) * 2 + 2) % 10 + (r[1] - 2 + 12) % 12) % 10;
        g.check(r[0] == wantGan, "天干合五虎遁");
        g.check(r[0] % 2 == r[1] % 2, "命宫干支同阴阳");

        // 独立复现"十四/廿六"口诀
        int m = gongNum(mz);
        if (az) m = m % 12 + 1;
        int sum = m + gongNum(hz);
        int want = ((sum < 14) ? (14 - sum) : (26 - sum));
        g.check(gongNum(r[1]) == want, "本数合口诀");
    }

    // B3 过中气 等价于 月支进一位
    for (int yg = 0; yg < 10; yg++)
    for (int mz = 0; mz < 12; mz++)
    for (int hz = 0; hz < 12; hz++) {
        g.ctx = std::string("年干") + Gan[yg] + " 月支" + Zhi[mz] + " 时支" + Zhi[hz];
        auto a = bazi::mingGong(yg, mz, hz, true);
        auto b = bazi::mingGong(yg, (mz + 1) % 12, hz, false);
        g.check(a[1] == b[1], "过中气 ≡ 次月未过中气(地支)");
    }

    // B4 非法入参
    g.ctx = "非法入参";
    g.check(bazi::taiYuanStr(-1, 0).empty(), "胎元 月干 -1");
    g.check(bazi::taiYuanStr(0, 12).empty(), "胎元 月支 12");
    g.check(bazi::mingGongStr(10, 0, 0, false).empty(), "命宫 年干 10");
    g.check(bazi::mingGongStr(0, -1, 0, false).empty(), "命宫 月支 -1");
    g.check(bazi::mingGongStr(0, 0, 99, false).empty(), "命宫 时支 99");
}

// ═══════════════════════════════════════════════════════════════════════
//  C. 与真实排盘联动
// ═══════════════════════════════════════════════════════════════════════
void testC(int y0, int y1) {
    Group &g = group("C 端到端联动");

    const LiFaType lifas[] = {XianDaiNongLifa_DingQiFa,
                              YuWuWeiZiPingLifa_DingDongZhi,
                              YuWuWeiZiPingLifa_DingXiaZhi};

    for (int Y = y0; Y <= y1; Y += 7)
    for (int M = 1; M <= 12; M += 1)
    for (LiFaType lf : lifas)
    for (int ast = 0; ast < 2; ast++) {
        // 每月取 3 天, 覆盖节前/中气前后/月末
        for (int D : {4, 12, 25}) {
            SBaziInputPara p;
            p.name = "T"; p.gender = false; p.isAst = (ast != 0);
            p.jw = {116.4, 39.9, "t", "t"};
            p.lifa = lf; p.calendar = CalendarSolar;
            p.birthDayTime = Time(Y, (uint8_t)M, (uint8_t)D, 10, 30, 0);

            BaziBase b(p);
            b.calcBaziPaiPan();

            char buf[96];
            std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d %s%s", Y, M, D,
                          lifaName(lf), ast ? " AST" : "");
            g.ctx = buf;

            auto sz = b.getSiZhuIndex();
            long double head = ast ? b.getHeadJieQiAstJd() : b.getHeadJieQiJd();
            long double tail = ast ? b.getTailJieQiAstJd() : b.getTailJieQiJd();
            long double mid  = ast ? b.getMidQiAstJd()     : b.getMidQiJd();

            // C1 中气严格落在两节之间
            g.check(mid > head && mid < tail,
                    "中气应在 (首节,末节) 内: head=" + JD::formatStr(head) +
                    " mid=" + JD::formatStr(mid) + " tail=" + JD::formatStr(tail));

            // C2 胎元 ≡ 月柱干+1 支+3
            g.check(b.getTaiYuan() == bazi::taiYuanStr(sz[2], sz[3]),
                    "胎元 得" + b.getTaiYuan() + " 期望" + bazi::taiYuanStr(sz[2], sz[3]));

            // C3 命宫 ≡ f(年干, 月支, 时支, 过中气)
            std::string wantMg = bazi::mingGongStr(sz[0], sz[3], sz[7], b.isAfterZhongQi());
            g.check(b.getMingGong() == wantMg,
                    "命宫 得" + b.getMingGong() + " 期望" + wantMg);

            // 过中气标志必须与 出生时刻 vs 中气 一致
            long double bd = ast ? b.getBirthAstJd() : b.getBirthJd();
            g.check(b.isAfterZhongQi() == (bd >= mid), "过中气标志自洽");

            // C5 平气法: 24 等分 ⇒ 中气居两节正中。
            // 两点例外, 故留 2 分钟容差而非要求严格相等:
            //  · 气年接缝月的"月末节"取自相邻气年的网格(间隔 k1≠k), 两半差几秒;
            //  · 真太阳时是平时的非线性变换(均时差), 中点不守恒, 半月跨度上可差一分余。
            // 容差仍足以捕获"按固定 15 天估中气"或"用错 k"这类实现错误。
            if (lf != XianDaiNongLifa_DingQiFa) {
                long double half = (head + tail) / 2.0L;
                double devSec = (double)(mid - half) * 86400.0;
                g.check(std::fabs(devSec) < 120.0,
                        "平气中气应近两节正中, 偏差 " + std::to_string(devSec) + "s");
                double frac = (double)((mid - head) / (tail - head));
                g.check(frac > 0.45 && frac < 0.55,
                        "中气应在月内 45%~55% 处, 实为 " + std::to_string(frac));
            }
        }
    }
}

// C4 —— 过中气标志沿时间轴必须单调翻转一次(不得反复横跳)
void testC4(int y0, int y1) {
    Group &g = group("C4 过中气标志单调性");

    for (int Y = y0; Y <= y1; Y += 11)
    for (int M = 1; M <= 12; M += 3) {
        // 以该月 1 日为起点, 逐 6 小时扫描 28 天, 记录标志翻转次数。
        // 期间可能跨一个节(标志由 true 复位为 false), 故允许 0~2 次翻转,
        // 但每次翻转必须与"首节/中气"的位置解释得通。
        int flips = 0;
        bool prev = false, first = true;
        long double prevHead = 0;
        for (int step = 0; step < 28 * 4; step++) {
            int hh = (step % 4) * 6;
            int dd = 1 + step / 4;
            if (dd > 28) break;
            SBaziInputPara p;
            p.name = "T"; p.gender = false; p.isAst = false;
            p.jw = {120.0, 39.9, "t", "t"};
            p.lifa = XianDaiNongLifa_DingQiFa; p.calendar = CalendarSolar;
            p.birthDayTime = Time(Y, (uint8_t)M, (uint8_t)dd, (uint8_t)hh, 0, 0);
            BaziBase b(p);
            b.calcBaziPaiPan();

            bool cur = b.isAfterZhongQi();
            long double head = b.getHeadJieQiJd();
            char buf[64];
            std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d %02d:00", Y, M, dd, hh);
            g.ctx = buf;

            if (!first) {
                bool sameMonth = std::fabs((double)(head - prevHead)) < 1e-6;
                if (cur != prev) {
                    flips++;
                    // 同一个节内只允许 false → true
                    if (sameMonth)
                        g.check(prev == false && cur == true,
                                "同月内标志只能 未过→已过");
                }
                if (sameMonth && prev == true)
                    g.check(cur == true, "同月内已过中气后不得回落");
            }
            prev = cur; prevHead = head; first = false;
        }
        g.check(flips <= 3, "28 天内翻转次数 " + std::to_string(flips) + " 应 ≤3");
    }
}

} // namespace

int main(int argc, char **argv) {
    int y0 = 1905, y1 = 2100;
    for (int i = 1; i < argc; i++) {
        if (!std::strcmp(argv[i], "-v")) gVerbose = true;
        else if (!std::strcmp(argv[i], "-range") && i + 2 < argc) {
            y0 = std::atoi(argv[i + 1]); y1 = std::atoi(argv[i + 2]); i += 2;
        }
    }

    std::printf("胎元 / 命宫 规格符合性测试  年份范围 %d..%d\n", y0, y1);
    std::printf("─────────────────────────────────────────────────────\n");

    testA();
    testB();
    testC(y0, y1);
    testC4(y0, y1);

    int totalPass = 0, totalFail = 0;
    for (auto *g : gGroups) {
        std::printf("%-24s  pass=%-7d fail=%d\n", g->name, g->pass, g->fail);
        for (const auto &m : g->msgs) std::printf("      ✗ %s\n", m.c_str());
        totalPass += g->pass; totalFail += g->fail;
    }
    std::printf("─────────────────────────────────────────────────────\n");
    std::printf("总计 pass=%d fail=%d  → %s\n", totalPass, totalFail,
                totalFail == 0 ? "PASS" : "FAIL");
    return totalFail == 0 ? 0 : 1;
}
