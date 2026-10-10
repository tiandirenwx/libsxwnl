// ─────────────────────────────────────────────────────────────────────────────
// 八字排盘 —— 规格符合性 / 黑盒不变量测试
//
// 本测试不看实现, 只按"外部可陈述的规格"检查 BaziBase 的对外行为。
// 每组用例的判据都来自独立来源(平气法定义 / 天文黄经 / 传统干支口诀 /
// 公开发布的节气时刻), 而不是从 bazi.cpp 里反抄出来的中间量。
//
// ── A. 平气网格符合定义 ───────────────────────────────────────────────
//    平气法(尤氏子平历)的定义: 取一个真至点为起点, 把该至点到下一个同名至点
//    的整段时间 24 等分, 12 个"节"落在奇数分点上。定冬至以真冬至为锚,
//    定夏至以真夏至为锚。
//    参考网格由 qi_accurate() 直接反解太阳视黄经 270°/90° 得到, 不经过
//    SSQ 历谱表, 也不经过 qi_accurate2() 的就近吸附, 与被测路径相互独立。
//    A1 气年锚点 == 真至点
//    A2 12 个节的时刻 == 参考网格的奇数分点
//    A3 节 → 月支 的映射 (小寒→丑 … 大雪→子)
//    A4 立春位于至点后 3/24 (定冬至) 或 15/24 (定夏至) 个气年, 年柱在此刻 +1
//
// ── B. 节气链自洽 ────────────────────────────────────────────────────
//    B1 出生时刻必须落在 [本月首节, 本月末节) 内
//    B2 同一个节, 从它所辖区间内任意时刻查询, 得到的首节/末节必须完全一致
//       (跨公历年、跨气年的节最容易在这里露馅)
//    B3 前一个节的"末节" 必须与后一个节的"首节" 完全相同 —— 节与节之间
//       不允许有缝隙或重叠
//
// ── C. 真太阳时 ──────────────────────────────────────────────────────
//    节气是全球同一瞬间发生的天文事件, 换算真太阳时只是换了一块钟的读数,
//    不改变事件先后。由此:
//    C1 同一物理时刻, 年柱/月柱 不得随 isAst 或出生地经度改变
//    C2 月柱(以及立春处的年柱)切换的物理瞬间不得随 isAst 改变
//    C3 真太阳时口径下: 日柱在真太阳时 23:00 切换, 时柱在真太阳时奇数整点切换
//    C4 真太阳时口径的节气自洽: headAst == calcAST(head), 且
//       headAst <= 真太阳时出生时刻 < tailAst
//    C5 节令交接前后逐秒扫描, 平时口径与真太阳时口径的月柱必须逐秒一致
//
// ── D. SSQ 单例无副作用 ──────────────────────────────────────────────
//    排盘需要临时把 SSQ 的年份窗口挪到别处, 挪完必须还原, 否则其它模块
//    (年历/月历/日详情) 会读到错年的农历表。
//    D1 排盘前后 SSQ 全部对外可见状态逐项相同
//    D2 语义级: 同一天的农历/干支/节气查询, 在中间插入远年排盘后结果不变
//    D3 复用已缓存的 Day 对象 (走缓存短路分支) 时同样不变
//
// ── E. 外部锚点 ──────────────────────────────────────────────────────
//    E1 公开发布的节气时刻 (北京时) 对照 —— 定气模式
//    E2 公认的日干支锚点 2000-01-01 = 戊午日, 按公历日序独立推演
//    E3 1984 = 甲子年 (立春为界) 独立推演年柱
//
// ── F. 干支组合口诀 ──────────────────────────────────────────────────
//    对每一张排出的盘施加传统口诀, 与历法无关, 是最强的横向约束:
//    F1 五虎遁(年上起月): 月干 = (年干 mod 5)*2 + 2 + (月支 - 寅)  mod 10
//    F2 五鼠遁(日上起时): 时干 = (日干 mod 5)*2 + 时支            mod 10
//    F3 月支/时支 与干的奇偶必须同性 (阳干配阳支)
//
// 编译: make bazi_invariants_test
// 运行:
//   ./build/bin/bazi_invariants_test              # 默认采样年份
//   ./build/bin/bazi_invariants_test -v           # 打印明细
//   ./build/bin/bazi_invariants_test -y 2024      # 单年
//   ./build/bin/bazi_invariants_test -range 1900 2100
//   ./build/bin/bazi_invariants_test -g A,C       # 只跑指定组
// ─────────────────────────────────────────────────────────────────────────────

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <sstream>
#include <string>
#include <vector>

#include "bazi.h"
#include "eph.h"
#include "JD.h"
#include "sxtwl.h"
#include "const.h"

namespace {

// ═══════════════════════════════════════════════════════════════════════
//  最小测试框架
// ═══════════════════════════════════════════════════════════════════════

struct Reporter
{
    int total = 0;
    int passed = 0;
    std::vector<std::string> failures;
    std::string group;      // 当前分组
    std::string ctx;        // 当前用例上下文

    void pass() { ++total; ++passed; }
    void fail(const std::string &what)
    {
        ++total;
        if (failures.size() < 120)
            failures.push_back("[" + group + "] " + ctx + " :: " + what);
        else if (failures.size() == 120)
            failures.push_back("... (后续失败已截断)");
    }
    void check(bool cond, const std::string &what)
    {
        if (cond) pass(); else fail(what);
    }
};

Reporter g;

#define CHECK(cond, msg) g.check((cond), (msg))

std::string fmtJd(long double jd)
{
    Time t = JD::JD2DD(jd);
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%d-%02d-%02d %02d:%02d:%06.3f",
                  t.Y, (int)t.M, (int)t.D, (int)t.h, (int)t.m, (double)t.s);
    return buf;
}

// 某个 JD 读数在当日内的秒数 (0..86400)
double secondOfDay(long double jd)
{
    Time t = JD::JD2DD(jd);
    return t.h * 3600.0 + t.m * 60.0 + t.s;
}

std::string fmtSec(long double days)
{
    char buf[48];
    std::snprintf(buf, sizeof(buf), "%.3f 秒", (double)(days * 86400.0L));
    return buf;
}

const char *const kGanStr[10] = {"甲","乙","丙","丁","戊","己","庚","辛","壬","癸"};
const char *const kZhiStr[12] = {"子","丑","寅","卯","辰","巳","午","未","申","酉","戌","亥"};

std::string gz(int g_, int z_)
{
    if (g_ < 0 || g_ >= 10 || z_ < 0 || z_ >= 12) return "??";
    return std::string(kGanStr[g_]) + kZhiStr[z_];
}

std::string lifaName(LiFaType l)
{
    switch (l)
    {
    case YuWuWeiZiPingLifa_DingDongZhi: return "定冬至";
    case YuWuWeiZiPingLifa_DingXiaZhi:  return "定夏至";
    case XianDaiNongLifa_DingQiFa:      return "定气";
    default: return "?";
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  被测系统的黑盒包装
// ═══════════════════════════════════════════════════════════════════════

struct Chart
{
    bool ok = false;
    std::string err;
    std::array<int, 8> sz{};      // 年干 年支 月干 月支 日干 日支 时干 时支
    long double head{}, tail{};          // 平时口径 月首节 / 月末节
    long double headAst{}, tailAst{};    // 真太阳时口径
    long double birth{}, birthAst{};
    int daysAfterJie = -1;
};

Chart cast(long double jd, LiFaType lifa, bool isAst = false, double lon = 120.0)
{
    Chart c;
    try
    {
        SBaziInputPara p;
        p.name = "T";
        p.gender = false;
        p.isAst = isAst;
        p.jw = {lon, 39.9, "test", "test"};
        p.lifa = lifa;
        p.calendar = CalendarSolar;
        p.birthDayTime = JD::JD2DD(jd);
        p.isRun = false;
        p.isSpec = false;
        BaziBase b(p);
        b.calcBaziPaiPan();
        c.sz = b.getSiZhuIndex();
        c.head = b.getHeadJieQiJd();
        c.tail = b.getTailJieQiJd();
        c.headAst = b.getHeadJieQiAstJd();
        c.tailAst = b.getTailJieQiAstJd();
        c.birth = b.getBirthJd();
        c.birthAst = b.getBirthAstJd();
        c.daysAfterJie = b.getDaysAfterJie();
        c.ok = true;
    }
    catch (const std::exception &e) { c.err = e.what(); }
    catch (...)                     { c.err = "unknown"; }
    return c;
}

// 二分定位"某一柱"发生变化的物理瞬间, 返回切换后的第一个 JD
template <typename Pred>
long double bisectSwitch(long double lo, long double hi, Pred sameAsLo,
                         long double epsDays = 1e-9L)
{
    while (hi - lo > epsDays)
    {
        long double mid = (lo + hi) / 2;
        if (sameAsLo(mid)) lo = mid; else hi = mid;
    }
    return hi;
}

// ═══════════════════════════════════════════════════════════════════════
//  独立参考实现: 真至点 与 平气网格
//
//  qi_accurate(W) 直接对"太阳视黄经 == W"求根 (W 为累计弧度, 不取模)。
//  以 1999 春分 为 0 号 15° 分点, 则:
//      Y 年 12 月冬至 (270°) 的分点序号 = 42 + 24*(Y-2000)
//      Y 年  6 月夏至 ( 90°) 的分点序号 = 30 + 24*(Y-2000)
//  这条路径不经 SSQ 历谱整日表, 也不经 qi_accurate2 的 ±5 天就近吸附。
// ═══════════════════════════════════════════════════════════════════════

long double refSolstice(int Y, bool winter)
{
    int n = (winter ? 42 : 30) + 24 * (Y - 2000);
    return qi_accurate(n * (PI / 12.0L), false, 120.0L) + J2000;
}

// 一个气年的 25 个分点 (T[0]=起点至点, T[24]=下一个同名至点)
struct RefFrame
{
    long double T[25];
    long double step;   // 单个气的长度 = 气年/24
};

// 定冬至: 气年 Y = [Y-1 年冬至, Y 年冬至]
// 定夏至: 气年 Y = [Y   年夏至, Y+1 年夏至]
RefFrame refFrame(int Y, LiFaType lifa)
{
    bool dz = (lifa == YuWuWeiZiPingLifa_DingDongZhi);
    long double a0 = dz ? refSolstice(Y - 1, true) : refSolstice(Y, false);
    long double a1 = dz ? refSolstice(Y, true)     : refSolstice(Y + 1, false);
    RefFrame f{};
    f.step = (a1 - a0) / 24.0L;
    for (int i = 0; i <= 24; ++i) f.T[i] = a0 + f.step * i;
    return f;
}

// 分点序号 i (奇数 = 节) 对应的月支
int refZhiOfNode(int i, LiFaType lifa)
{
    bool dz = (lifa == YuWuWeiZiPingLifa_DingDongZhi);
    // 定冬至: T1 小寒→丑(1), T3 立春→寅(2) ... T23 大雪→子(0)
    // 定夏至: T1 小暑→未(7), T3 立秋→申(8) ... T23 芒种→午(6)
    return ((i - 1) / 2 + (dz ? 1 : 7)) % 12;
}

// ═══════════════════════════════════════════════════════════════════════
//  F 组: 干支口诀 —— 对每一张盘都成立, 与历法无关
// ═══════════════════════════════════════════════════════════════════════

void checkGanZhiRules(const Chart &c, const std::string &where)
{
    if (!c.ok) { CHECK(false, where + " 排盘抛异常: " + c.err); return; }
    int yG = c.sz[0], yZ = c.sz[1], mG = c.sz[2], mZ = c.sz[3];
    int dG = c.sz[4], dZ = c.sz[5], hG = c.sz[6], hZ = c.sz[7];

    bool inRange = yG >= 0 && yG < 10 && mG >= 0 && mG < 10 && dG >= 0 && dG < 10 &&
                   hG >= 0 && hG < 10 && yZ >= 0 && yZ < 12 && mZ >= 0 && mZ < 12 &&
                   dZ >= 0 && dZ < 12 && hZ >= 0 && hZ < 12;
    CHECK(inRange, where + " 四柱索引越界");
    if (!inRange) return;

    // F3 阳干只配阳支 (干支序号同奇偶)
    CHECK((yG & 1) == (yZ & 1), where + " 年柱干支阴阳不配: " + gz(yG, yZ));
    CHECK((mG & 1) == (mZ & 1), where + " 月柱干支阴阳不配: " + gz(mG, mZ));
    CHECK((dG & 1) == (dZ & 1), where + " 日柱干支阴阳不配: " + gz(dG, dZ));
    CHECK((hG & 1) == (hZ & 1), where + " 时柱干支阴阳不配: " + gz(hG, hZ));

    // F1 五虎遁: 甲己之年丙作首 —— 寅月天干 = (年干%5)*2+2
    int expMG = ((yG % 5) * 2 + 2 + ((mZ - 2) % 12 + 12) % 12) % 10;
    CHECK(mG == expMG,
          where + " 五虎遁不符: 年=" + gz(yG, yZ) + " 月=" + gz(mG, mZ) +
          " 期望月干=" + kGanStr[expMG]);

    // F2 五鼠遁: 甲己还加甲 —— 子时天干 = (日干%5)*2
    int expHG = ((dG % 5) * 2 + hZ) % 10;
    CHECK(hG == expHG,
          where + " 五鼠遁不符: 日=" + gz(dG, dZ) + " 时=" + gz(hG, hZ) +
          " 期望时干=" + kGanStr[expHG]);
}

// ═══════════════════════════════════════════════════════════════════════
//  A 组 + B 组: 平气网格符合定义 & 节气链自洽
// ═══════════════════════════════════════════════════════════════════════

// 节气时刻判定的容差。参考网格与被测实现走的是同一套 VSOP 太阳黄经模型,
// 差异只应来自浮点舍入, 故取 0.2 秒 —— 任何"算法层面"的偏差都会被抓住。
const long double kNodeTol = 0.2L / 86400.0L;

void groupA_B(int Y, LiFaType lifa, bool verbose)
{
    g.group = "A/B " + lifaName(lifa);
    RefFrame cur = refFrame(Y, lifa);
    RefFrame prv = refFrame(Y - 1, lifa);
    RefFrame nxt = refFrame(Y + 1, lifa);

    // A1: 气年锚点必须是真至点。用气年首节反推: 首节 = 至点 + 气年/24,
    //     这里直接验证参考网格与被测实现在"节"上的一致性即可 (A2), 至点本身
    //     由 refSolstice 独立给出, 若被测实现锚错了至点, A2 必然大面积失败。

    const long double eps = 1.0L / 86400.0L;   // 1 秒

    for (int i = 1; i <= 23; i += 2)
    {
        long double node  = cur.T[i];
        long double nodeN = (i == 23) ? nxt.T[1]  : cur.T[i + 2];   // 下一个节
        long double nodeP = (i == 1)  ? prv.T[23] : cur.T[i - 2];   // 上一个节
        int zhi     = refZhiOfNode(i, lifa);
        int zhiPrev = (zhi + 11) % 12;

        std::ostringstream base;
        base << "Y=" << Y << " 分点#" << i << "(" << kZhiStr[zhi] << "月首节 "
             << fmtJd(node) << ")";

        // ── 节后 1 秒 ──
        {
            g.ctx = base.str() + " 节后1秒";
            Chart c = cast(node + eps, lifa);
            if (!c.ok) { CHECK(false, "排盘失败: " + c.err); continue; }

            // A3 月支
            CHECK(c.sz[3] == zhi,
                  "月支应为 " + std::string(kZhiStr[zhi]) + " 实为 " + kZhiStr[c.sz[3]]);
            // A2 首节 == 参考分点
            CHECK(std::fabs((double)(c.head - node)) < kNodeTol,
                  "首节偏离参考网格 " + fmtSec(c.head - node) + " (实=" + fmtJd(c.head) + ")");
            // A2 末节 == 下一个参考分点
            CHECK(std::fabs((double)(c.tail - nodeN)) < kNodeTol,
                  "末节偏离参考网格 " + fmtSec(c.tail - nodeN) + " (实=" + fmtJd(c.tail) + ")");
            // B1 出生落在区间内
            CHECK(c.head <= c.birth && c.birth < c.tail, "出生时刻未落在 [首节,末节) 内");
            checkGanZhiRules(c, "");
        }

        // ── 节前 1 秒 ──
        {
            g.ctx = base.str() + " 节前1秒";
            Chart c = cast(node - eps, lifa);
            if (!c.ok) { CHECK(false, "排盘失败: " + c.err); continue; }

            CHECK(c.sz[3] == zhiPrev,
                  "月支应为 " + std::string(kZhiStr[zhiPrev]) + " 实为 " + kZhiStr[c.sz[3]]);
            // B3 上一个节的末节必须正是本节
            CHECK(std::fabs((double)(c.tail - node)) < kNodeTol,
                  "末节偏离参考网格 " + fmtSec(c.tail - node) + " (实=" + fmtJd(c.tail) + ")");
            // A2 上一个节 (跨气年时必须用上一气年的间隔)
            CHECK(std::fabs((double)(c.head - nodeP)) < kNodeTol,
                  "首节偏离参考网格 " + fmtSec(c.head - nodeP) + " (实=" + fmtJd(c.head) +
                  " 参考=" + fmtJd(nodeP) + ")");
            CHECK(c.head <= c.birth && c.birth < c.tail, "出生时刻未落在 [首节,末节) 内");
            checkGanZhiRules(c, "");
        }

        // ── B2: 同一节区间内多点取样, 首/末节必须完全一致 ──
        // 这是最贴近真实使用的一条: 同一个节里出生的两个人, 排出来的节令时刻
        // 不能因为出生日期落在公历不同月/不同年而不同。
        {
            g.ctx = base.str() + " 区间内多点取样";
            long double span = nodeN - node;
            Chart ref0 = cast(node + eps, lifa);
            const double fr[5] = {0.05, 0.25, 0.5, 0.75, 0.97};
            for (double r : fr)
            {
                Chart c = cast(node + span * (long double)r, lifa);
                if (!c.ok || !ref0.ok) { CHECK(false, "排盘失败"); break; }
                CHECK(c.sz[3] == ref0.sz[3],
                      "月支随取样点变化: " + std::string(kZhiStr[ref0.sz[3]]) + " → " +
                      kZhiStr[c.sz[3]] + " @r=" + std::to_string(r));
                CHECK(std::fabs((double)(c.head - ref0.head)) < kNodeTol,
                      "首节随取样点变化 " + fmtSec(c.head - ref0.head) +
                      " @r=" + std::to_string(r) +
                      " (" + fmtJd(ref0.head) + " → " + fmtJd(c.head) + ")");
                CHECK(std::fabs((double)(c.tail - ref0.tail)) < kNodeTol,
                      "末节随取样点变化 " + fmtSec(c.tail - ref0.tail) +
                      " @r=" + std::to_string(r));
            }
        }

        if (verbose)
            std::printf("    #%2d %s → %s月\n", i, fmtJd(node).c_str(), kZhiStr[zhi]);
    }

    // ── A4: 立春 ──
    {
        int li = (lifa == YuWuWeiZiPingLifa_DingXiaZhi) ? 15 : 3;
        long double lc = cur.T[li];
        g.ctx = "Y=" + std::to_string(Y) + " 立春 " + fmtJd(lc);
        Chart a = cast(lc - eps, lifa);
        Chart b = cast(lc + eps, lifa);
        if (a.ok && b.ok)
        {
            int dG = ((b.sz[0] - a.sz[0]) % 10 + 10) % 10;
            int dZ = ((b.sz[1] - a.sz[1]) % 12 + 12) % 12;
            CHECK(dG == 1 && dZ == 1,
                  "年柱未在立春递增 1: " + gz(a.sz[0], a.sz[1]) + " → " + gz(b.sz[0], b.sz[1]));
            CHECK(b.sz[3] == 2, "立春后月支应为寅, 实为 " + std::string(kZhiStr[b.sz[3]]));
            CHECK(a.sz[3] == 1, "立春前月支应为丑, 实为 " + std::string(kZhiStr[a.sz[3]]));
        }
        else CHECK(false, "立春处排盘失败");
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  C 组: 真太阳时
// ═══════════════════════════════════════════════════════════════════════

void groupC(int Y, LiFaType lifa, bool verbose)
{
    g.group = "C " + lifaName(lifa);
    const double lons[] = {75.0, 90.0, 105.0, 116.4, 120.0, 135.0};
    const long double eps = 1.0L / 86400.0L;

    // 取本年若干个物理时刻做横向对照, 其中包含节令交接瞬间
    std::vector<long double> probes;
    {
        Time t0{Y, 1, 1, 12, 0, 0};
        long double jd0 = JD::toJD(t0);
        for (int d : {5, 40, 95, 150, 200, 260, 320, 358})
            probes.push_back(jd0 + d + 0.37L);
    }
    RefFrame cur = refFrame(Y, lifa);
    if (lifa != XianDaiNongLifa_DingQiFa)
        for (int i = 1; i <= 23; i += 2)
        {
            probes.push_back(cur.T[i] - eps);
            probes.push_back(cur.T[i] + eps);
        }

    for (long double jd : probes)
    {
        g.ctx = "Y=" + std::to_string(Y) + " t=" + fmtJd(jd);

        Chart basis = cast(jd, lifa, false, 120.0);
        if (!basis.ok) { CHECK(false, "基准排盘失败: " + basis.err); continue; }

        for (double lon : lons)
        {
            Chart plain = cast(jd, lifa, false, lon);
            Chart ast   = cast(jd, lifa, true,  lon);
            if (!plain.ok || !ast.ok) { CHECK(false, "排盘失败 lon=" + std::to_string(lon)); continue; }

            // C1 年柱/月柱 是物理事件的产物, 不得随 isAst / 经度改变
            CHECK(plain.sz[0] == basis.sz[0] && plain.sz[1] == basis.sz[1] &&
                  plain.sz[2] == basis.sz[2] && plain.sz[3] == basis.sz[3],
                  "平时口径年月柱随经度改变 lon=" + std::to_string(lon));
            CHECK(ast.sz[0] == basis.sz[0] && ast.sz[1] == basis.sz[1] &&
                  ast.sz[2] == basis.sz[2] && ast.sz[3] == basis.sz[3],
                  "真太阳时口径年月柱与平时口径不一致 lon=" + std::to_string(lon) +
                  " 平时=" + gz(basis.sz[0], basis.sz[1]) + " " + gz(basis.sz[2], basis.sz[3]) +
                  " 真太阳时=" + gz(ast.sz[0], ast.sz[1]) + " " + gz(ast.sz[2], ast.sz[3]));
            // 平时口径下, 经度不应影响任何一柱
            CHECK(plain.sz == basis.sz, "平时口径四柱随经度改变 lon=" + std::to_string(lon));

            // C4 真太阳时口径的节气自洽
            long double expHead = JD::calcAST(ast.head, lon);
            long double expTail = JD::calcAST(ast.tail, lon);
            CHECK(std::fabs((double)(ast.headAst - expHead)) < 1e-9,
                  "headAst != calcAST(head) lon=" + std::to_string(lon));
            CHECK(std::fabs((double)(ast.tailAst - expTail)) < 1e-9,
                  "tailAst != calcAST(tail) lon=" + std::to_string(lon));
            CHECK(ast.headAst <= ast.birthAst && ast.birthAst < ast.tailAst,
                  "真太阳时口径下出生时刻未落在 [首节,末节) 内 lon=" + std::to_string(lon) +
                  " head=" + fmtJd(ast.headAst) + " birth=" + fmtJd(ast.birthAst) +
                  " tail=" + fmtJd(ast.tailAst));

            // "距节令天数"必须与所选时间口径一致 —— 真太阳时排盘就该用真太阳时读数计天
            CHECK(ast.daysAfterJie == (int)std::floor((double)(ast.birthAst - ast.headAst)),
                  "真太阳时口径 getDaysAfterJie 与 floor(birthAst-headAst) 不符 lon=" +
                  std::to_string(lon) + " 实=" + std::to_string(ast.daysAfterJie));
            CHECK(plain.daysAfterJie == (int)std::floor((double)(plain.birth - plain.head)),
                  "平时口径 getDaysAfterJie 与 floor(birth-head) 不符");

            checkGanZhiRules(ast, "lon=" + std::to_string(lon) + " AST");
        }
    }

    // C2 月柱切换的物理瞬间不得随 isAst 改变
    {
        Time t0{Y, 3, 20, 12, 0, 0};
        long double lo = JD::toJD(t0), hi = lo + 40;
        Chart a0 = cast(lo, lifa, false, 105.0);
        Chart a1 = cast(hi, lifa, false, 105.0);
        if (a0.ok && a1.ok && (a0.sz[2] != a1.sz[2] || a0.sz[3] != a1.sz[3]))
        {
            int g0 = a0.sz[2], z0 = a0.sz[3];
            long double swPlain = bisectSwitch(lo, hi, [&](long double m) {
                Chart c = cast(m, lifa, false, 105.0);
                return c.ok && c.sz[2] == g0 && c.sz[3] == z0;
            });
            long double swAst = bisectSwitch(lo, hi, [&](long double m) {
                Chart c = cast(m, lifa, true, 105.0);
                return c.ok && c.sz[2] == g0 && c.sz[3] == z0;
            });
            g.ctx = "Y=" + std::to_string(Y) + " 月柱切换点";
            CHECK(std::fabs((double)(swPlain - swAst)) < 1e-7,
                  "月柱切换瞬间随 isAst 改变: 平时=" + fmtJd(swPlain) +
                  " 真太阳时=" + fmtJd(swAst) + " 差 " + fmtSec(swAst - swPlain));
        }
    }

    // C3 真太阳时口径: 日柱在真太阳时 23:00 切换, 时柱在真太阳时奇数整点切换
    for (double lon : {87.6, 120.0, 130.0})
    {
        Time t0{Y, 8, 10, 12, 0, 0};
        long double jd0 = JD::toJD(t0);
        g.ctx = "Y=" + std::to_string(Y) + " lon=" + std::to_string(lon) + " 日/时柱边界";

        // 日柱: 在 [jd0, jd0+2) 内找一次切换
        Chart d0 = cast(jd0, lifa, true, lon);
        Chart d1 = cast(jd0 + 1.5L, lifa, true, lon);
        if (d0.ok && d1.ok && (d0.sz[4] != d1.sz[4] || d0.sz[5] != d1.sz[5]))
        {
            int gg = d0.sz[4], zz = d0.sz[5];
            long double sw = bisectSwitch(jd0, jd0 + 1.5L, [&](long double m) {
                Chart c = cast(m, lifa, true, lon);
                return c.ok && c.sz[4] == gg && c.sz[5] == zz;
            });
            // 二分残差使读数可能是 22:59:59.999, 故先就近到秒再判定
            long r = std::lround(secondOfDay(JD::calcAST(sw, lon)));
            CHECK(((r % 86400) + 86400) % 86400 == 23 * 3600,
                  "日柱切换点的真太阳时读数不是 23:00:00, 实为 " +
                  fmtJd(JD::calcAST(sw, lon)));
        }

        // 时柱: 在 [jd0, jd0+0.2) 内找一次切换
        Chart h0 = cast(jd0, lifa, true, lon);
        Chart h1 = cast(jd0 + 0.2L, lifa, true, lon);
        if (h0.ok && h1.ok && (h0.sz[6] != h1.sz[6] || h0.sz[7] != h1.sz[7]))
        {
            int gg = h0.sz[6], zz = h0.sz[7];
            long double sw = bisectSwitch(jd0, jd0 + 0.2L, [&](long double m) {
                Chart c = cast(m, lifa, true, lon);
                return c.ok && c.sz[6] == gg && c.sz[7] == zz;
            });
            long r = std::lround(secondOfDay(JD::calcAST(sw, lon)));
            long hour = ((r / 3600) % 24 + 24) % 24;
            CHECK(r % 3600 == 0 && hour % 2 == 1,
                  "时柱切换点的真太阳时读数不是奇数整点, 实为 " +
                  fmtJd(JD::calcAST(sw, lon)));
        }
    }

    // C5 节令交接窗口逐秒扫描: 平时/真太阳时 月柱必须逐秒一致
    if (lifa != XianDaiNongLifa_DingQiFa)
    {
        long double node = cur.T[3];   // 立春附近, 同时覆盖年柱切换
        g.ctx = "Y=" + std::to_string(Y) + " 立春 ±60 秒逐秒扫描";
        for (int s = -60; s <= 60; ++s)
        {
            long double jd = node + (long double)s / 86400.0L;
            Chart p = cast(jd, lifa, false, 100.0);
            Chart a = cast(jd, lifa, true,  100.0);
            if (!p.ok || !a.ok) { CHECK(false, "排盘失败 s=" + std::to_string(s)); break; }
            if (p.sz[0] != a.sz[0] || p.sz[1] != a.sz[1] ||
                p.sz[2] != a.sz[2] || p.sz[3] != a.sz[3])
            {
                CHECK(false, "s=" + std::to_string(s) + " 秒 平时年月柱=" +
                      gz(p.sz[0], p.sz[1]) + " " + gz(p.sz[2], p.sz[3]) +
                      " 真太阳时=" + gz(a.sz[0], a.sz[1]) + " " + gz(a.sz[2], a.sz[3]));
                break;
            }
        }
        g.pass();
    }
}

// ═══════════════════════════════════════════════════════════════════════
//  D 组: SSQ 单例无副作用
// ═══════════════════════════════════════════════════════════════════════

// 同一天的全部对外可见历法信息
struct DayFacts
{
    int lunarYearLichun = 0, lunarYearJun = 0;
    int lunarMonth = 0, lunarDay = 0, lunarMonthDays = 0;
    bool leap = false, nextIsZheng = false;
    int yGanLc = 0, yZhiLc = 0, mGan = 0, mZhi = 0, dGan = 0, dZhi = 0;
    int constellation = 0;
    int8_t jqIdx = 0;

    bool operator==(const DayFacts &o) const
    {
        return lunarYearLichun == o.lunarYearLichun && lunarYearJun == o.lunarYearJun &&
               lunarMonth == o.lunarMonth && lunarDay == o.lunarDay &&
               lunarMonthDays == o.lunarMonthDays && leap == o.leap &&
               nextIsZheng == o.nextIsZheng && yGanLc == o.yGanLc && yZhiLc == o.yZhiLc &&
               mGan == o.mGan && mZhi == o.mZhi && dGan == o.dGan && dZhi == o.dZhi &&
               constellation == o.constellation && jqIdx == o.jqIdx;
    }
    std::string str() const
    {
        std::ostringstream os;
        os << "农历年(立春界)=" << lunarYearLichun << " 农历年(春节界)=" << lunarYearJun
           << " 月=" << lunarMonth << (leap ? "闰" : "") << " 日=" << lunarDay
           << " 月大小=" << lunarMonthDays << " 下月正=" << nextIsZheng
           << " 年柱=" << gz(yGanLc, yZhiLc) << " 月柱=" << gz(mGan, mZhi)
           << " 日柱=" << gz(dGan, dZhi) << " 星座=" << constellation
           << " 节气idx=" << (int)jqIdx;
        return os.str();
    }
};

// reuse=false: 每次新建 Day 对象; reuse=true: 复用传入的对象(走内部缓存短路)
DayFacts readDayFacts(Day *d)
{
    DayFacts f;
    f.lunarYearLichun = d->getLunarYear(false);
    f.lunarYearJun    = d->getLunarYear(true);
    f.lunarMonth      = d->getLunarMonth();
    f.lunarDay        = d->getLunarDay();
    f.lunarMonthDays  = d->getLunarMonthDays();
    f.leap            = d->isLunarLeap();
    f.nextIsZheng     = d->isNextLunarMonthZheng();
    GZ y = d->getYearGZ(false); f.yGanLc = y.tg; f.yZhiLc = y.dz;
    GZ m = d->getMonthGZ();     f.mGan = m.tg;   f.mZhi = m.dz;
    GZ dd = d->getDayGZ();      f.dGan = dd.tg;  f.dZhi = dd.dz;
    f.constellation   = d->getConstellation();
    f.jqIdx           = d->hasLiPuJieQi() ? (int8_t)d->getLiPuJieQi() : (int8_t)-1;
    return f;
}

// ═══════════════════════════════════════════════════════════════════════
//  G 组: 报告文本自洽性
//  同一份排盘报告里出现的每个"干支年", 必须指向同一个立春。农历干支年若沿用
//  solar2Lunar2 的定气整日立春, 而年柱走所选历法的立春时刻, 报告就会出现
//  "农历庚子年 … 年柱己亥"这种自相矛盾。立春前后逐小时扫描来钉死这条。
// ═══════════════════════════════════════════════════════════════════════
void groupG(int Y, LiFaType lf, bool verbose)
{
    g.group = "G 报告自洽";

    for (int D = 2; D <= 7; ++D)
    for (int H = 0; H < 24; H += 3)
    {
        SBaziInputPara p;
        p.name = "T"; p.gender = false; p.isAst = false;
        p.jw = {120.0, 39.9, "test", "test"};
        p.lifa = lf; p.calendar = CalendarSolar;
        p.birthDayTime = Time(Y, 2, D, H, 0, 0);
        p.isRun = false; p.isSpec = false;

        std::array<int, 8> sz{};
        std::string lunarStr, sxStr, sxRef;
        try
        {
            BaziBase b(p);
            b.calcBaziPaiPan();
            sz = b.getSiZhuIndex();
            auto info = b.getLunarInfo();
            sxStr    = std::get<1>(info);
            lunarStr = std::get<2>(info);
            sxRef    = b.getShenXiao();
        }
        catch (...) { continue; }

        char buf[64];
        std::snprintf(buf, sizeof buf, "%d-02-%02d %02d:00 %s", Y, D, H, lifaName(lf).c_str());
        g.ctx = buf;

        // G1: 农历生日里的干支年 == 四柱年柱
        std::string want = "农历" + gz(sz[0], sz[1]) + "年";
        CHECK(lunarStr.rfind(want, 0) == 0,
              "农历干支年与年柱不一致: 年柱=" + gz(sz[0], sz[1]) + ", 农历串=\"" + lunarStr + "\"");

        // G2: getLunarInfo 的生肖 == getShenXiao()
        CHECK(sxStr == sxRef,
              "getLunarInfo 生肖(" + sxStr + ") 与 getShenXiao(" + sxRef + ") 不一致");
    }
}

// 排一张盘, 返回四柱字符串; 供"排盘结果与 SSQ 入态无关"的检查用
static std::string paipanSiZhu(int y, int m, int d, LiFaType lf, bool ast)
{
    SBaziInputPara p;
    p.name = "T"; p.gender = false; p.isAst = ast;
    p.jw = {116.4, 39.9, "t", "t"};
    p.lifa = lf; p.calendar = CalendarSolar;
    p.birthDayTime = Time(y, m, d, 10, 30, 0);
    p.isRun = false; p.isSpec = false;
    BaziBase b(p);
    b.calcBaziPaiPan();
    const auto &a = b.getSiZhuIndex();
    std::ostringstream os;
    for (int i = 0; i < 8; ++i) os << (int)a[i] << ",";
    return os.str();
}

// SSQ 单例是 thread_local 的共享缓存, 八字排盘会把它的"当前年份窗口"挪到出生年。
// 本组不断言单例的内部字节状态不变 —— 那是实现细节; 只断言用户可见的契约:
// **插入任意排盘, 不改变其它模块对同一天的查询结果**, 以及反过来
// **排盘结果不受 SSQ 入态影响**。读方(Day::checkSSQ 等)负责窗口不覆盖时重算。
void groupD(int Y, bool verbose)
{
    g.group = "D SSQ单例";

    // 让 SSQ 处于"Y 年窗口"这一确定初始状态
    Day *warm = sxtwl::fromSolar(Y, 6, 15);
    DayFacts baseline = readDayFacts(warm);

    // sxtwl 层(不经过 Day 缓存)的基线
    const uint8_t baseRun   = sxtwl::getRunMonth(Y);
    const uint8_t baseMDays = sxtwl::getLunarMonthNum(Y, 5, false);
    Day *baseLd = sxtwl::fromLunar(Y, 5, 12);
    auto ymd = [](Day *x) { return x->getSolarYear() * 10000 + (int)x->getSolarMonth() * 100 + x->getSolarDay(); };
    const int baseLdD0 = ymd(baseLd);
    delete baseLd;

    struct Case { int y, m, d; LiFaType lf; bool ast; };
    // 刻意挑选与 Y 相距很远的年份, 放大污染
    const std::vector<Case> cases = {
        {Y,      6, 15, YuWuWeiZiPingLifa_DingDongZhi, false},
        {Y,      6, 15, YuWuWeiZiPingLifa_DingXiaZhi,  false},
        {Y,      6, 15, XianDaiNongLifa_DingQiFa,      false},
        {Y - 137, 1,  3, YuWuWeiZiPingLifa_DingDongZhi, true},
        {Y + 211, 7,  3, YuWuWeiZiPingLifa_DingXiaZhi,  true},
        {Y + 400, 12, 28, XianDaiNongLifa_DingQiFa,     false},
    };

    for (const auto &cs : cases)
    {
        g.ctx = "Y=" + std::to_string(Y) + " 插入排盘 " + std::to_string(cs.y) + "-" +
                std::to_string(cs.m) + "-" + std::to_string(cs.d) + " " + lifaName(cs.lf);

        SBaziInputPara p;
        p.name = "T"; p.gender = false; p.isAst = cs.ast;
        p.jw = {116.4, 39.9, "t", "t"};
        p.lifa = cs.lf; p.calendar = CalendarSolar;
        p.birthDayTime = Time(cs.y, cs.m, cs.d, 10, 30, 0);
        p.isRun = false; p.isSpec = false;
        bool threw = false;
        try
        {
            BaziBase b(p);
            b.calcBaziPaiPan();
            // 这些只读接口内部也会用 SSQ, 同样不许留下痕迹
            (void)b.getAge();
            (void)b.getJieLing();
            (void)b.getQiYun();
        }
        catch (...) { threw = true; }
        CHECK(!threw, "排盘抛异常");

        // ── D1: 复用已缓存的 Day 对象 (走缓存短路分支) ──
        // 必须排在新建 Day 之前: 新建会触发 checkSSQ 把窗口拉回 Y 年, 先建就掩盖了污染。
        // 这条专盯"内部缓存已填 -> 提前 return -> 却又直接读 SSQ"的写法。
        DayFacts f2 = readDayFacts(warm);
        CHECK(f2 == baseline,
              "复用 Day 对象查询结果被污染\n        基线: " + baseline.str() +
              "\n        实得: " + f2.str());

        // ── D2: 新建 Day 查同一天, 结果必须与基线一致 ──
        Day *probe = sxtwl::fromSolar(Y, 6, 15);
        DayFacts f = readDayFacts(probe);
        CHECK(f == baseline,
              "新建 Day 查询结果被污染\n        基线: " + baseline.str() +
              "\n        实得: " + f.str());
        delete probe;

        // ── D3: 不经过 Day 的 sxtwl 层接口同样不受影响 ──
        CHECK(sxtwl::getRunMonth(Y) == baseRun, "getRunMonth 被污染");
        CHECK(sxtwl::getLunarMonthNum(Y, 5, false) == baseMDays, "getLunarMonthNum 被污染");
        Day *ld = sxtwl::fromLunar(Y, 5, 12);
        CHECK(ymd(ld) == baseLdD0, "fromLunar 被污染");
        delete ld;
    }

    // ── D4: 反向 —— 排盘结果不受 SSQ 入态影响 ──
    // 先把窗口分别置于三个相距很远的年份, 再排同一张盘, 四柱必须完全一致
    {
        g.ctx = "Y=" + std::to_string(Y) + " 排盘对 SSQ 入态免疫";
        for (LiFaType lf : {YuWuWeiZiPingLifa_DingDongZhi, YuWuWeiZiPingLifa_DingXiaZhi,
                            XianDaiNongLifa_DingQiFa})
        {
            std::string ref, got;
            for (int off : {0, -311, 523})
            {
                Day *w = sxtwl::fromSolar(Y + off, 3, 9);
                (void)w->getYearGZ(false);   // 强制把窗口挪到 Y+off
                delete w;
                got = paipanSiZhu(Y, 2, 4, lf, true);
                if (off == 0) ref = got;
                else CHECK(got == ref,
                           lifaName(lf) + " 排盘结果随 SSQ 入态变化 (窗口置于 " +
                           std::to_string(Y + off) + " 年): 期望 " + ref + " 实得 " + got);
            }
        }
    }
    delete warm;
}

// ═══════════════════════════════════════════════════════════════════════
//  E 组: 外部锚点
// ═══════════════════════════════════════════════════════════════════════

// 公历(格里高利)日序, 用于独立推算日柱; 仅用于 1583 年以后的日期
long long gregorianDayNumber(int y, int m, int d)
{
    long long a = (14 - m) / 12;
    long long yy = y + 4800 - a;
    long long mm = m + 12 * a - 3;
    return d + (153 * mm + 2) / 5 + 365 * yy + yy / 4 - yy / 100 + yy / 400 - 32045;
}

void groupE(bool verbose)
{
    g.group = "E 外部锚点";

    // ── E1: 公开发布的节气时刻 (北京时) ──
    //  取自公开历书。定气模式下, 排盘给出的"月首节/月末节"应与之吻合。
    //  容差 3 分钟 —— 足以区分算法错误, 又不苛求与某本历书末位秒一致。
    struct JieQi { int y, mo, d, h, mi; const char *name; };
    const std::vector<JieQi> known = {
        {2023, 12,  7, 17, 33, "2023 大雪"},
        {2024,  1,  6,  4, 49, "2024 小寒"},
        {2024,  2,  4, 16, 27, "2024 立春"},
        {2024,  3,  5, 10, 23, "2024 惊蛰"},
        {2024,  8,  7,  8,  9, "2024 立秋"},
        {2025,  2,  3, 22, 10, "2025 立春"},
    };
    const long double tol = 3.0L / 1440.0L;
    for (const auto &q : known)
    {
        g.ctx = q.name;
        Time t{q.y, (uint8_t)q.mo, (uint8_t)q.d, (double)q.h, (double)q.mi, 0};
        long double want = JD::toJD(t);
        // 该时刻之后 1 小时排盘, 本月首节就应当是这个节
        Chart c = cast(want + 1.0L / 24.0L, XianDaiNongLifa_DingQiFa);
        if (!c.ok) { CHECK(false, "排盘失败"); continue; }
        CHECK(std::fabs((double)(c.head - want)) < tol,
              "定气首节 " + fmtJd(c.head) + " 与历书 " + fmtJd(want) +
              " 相差 " + fmtSec(c.head - want));
        // 该时刻之前 1 小时排盘, 本月末节就应当是这个节
        Chart c2 = cast(want - 1.0L / 24.0L, XianDaiNongLifa_DingQiFa);
        if (!c2.ok) { CHECK(false, "排盘失败"); continue; }
        CHECK(std::fabs((double)(c2.tail - want)) < tol,
              "定气末节 " + fmtJd(c2.tail) + " 与历书 " + fmtJd(want) + " 不符");
    }

    // ── E2: 日柱锚点 2000-01-01 = 戊午日 ──
    //  戊=4, 午=6; 日柱每日 +1 (mod 60), 子时(23:00)起算下一日。
    const long long anchor = gregorianDayNumber(2000, 1, 1);
    struct DayCase { int y, m, d; };
    const std::vector<DayCase> dayCases = {
        {2000, 1, 1}, {1984, 2, 4}, {1949, 10, 1}, {2024, 2, 4},
        {2026, 9, 20}, {1900, 1, 1}, {2100, 6, 30},
    };
    for (const auto &dc : dayCases)
    {
        long long off = gregorianDayNumber(dc.y, dc.m, dc.d) - anchor;
        int expG = (int)(((4 + off) % 10 + 10) % 10);
        int expZ = (int)(((6 + off) % 12 + 12) % 12);
        for (LiFaType lf : {YuWuWeiZiPingLifa_DingDongZhi, YuWuWeiZiPingLifa_DingXiaZhi,
                            XianDaiNongLifa_DingQiFa})
        {
            g.ctx = std::to_string(dc.y) + "-" + std::to_string(dc.m) + "-" +
                    std::to_string(dc.d) + " 日柱 " + lifaName(lf);
            Time t{dc.y, (uint8_t)dc.m, (uint8_t)dc.d, 10, 0, 0};   // 巳时, 远离子时边界
            Chart c = cast(JD::toJD(t), lf);
            if (!c.ok) { CHECK(false, "排盘失败"); continue; }
            CHECK(c.sz[4] == expG && c.sz[5] == expZ,
                  "日柱应为 " + gz(expG, expZ) + " 实为 " + gz(c.sz[4], c.sz[5]));
        }
    }

    // ── E3: 年柱锚点 1984 = 甲子年 (立春为界) ──
    struct YearCase { int y, mo, d; int g_, z_; const char *why; };
    const std::vector<YearCase> yearCases = {
        {1984,  6, 1, 0, 0, "1984 立春后 = 甲子"},
        {1984,  1, 1, 9, 11, "1984 立春前 = 癸亥(1983年柱)"},
        {2024,  6, 1, 0, 4, "2024 立春后 = 甲辰"},
        {2000,  6, 1, 6, 4, "2000 立春后 = 庚辰"},
        {2026,  6, 1, 2, 6, "2026 立春后 = 丙午"},
    };
    for (const auto &yc : yearCases)
    {
        for (LiFaType lf : {YuWuWeiZiPingLifa_DingDongZhi, YuWuWeiZiPingLifa_DingXiaZhi,
                            XianDaiNongLifa_DingQiFa})
        {
            g.ctx = std::string(yc.why) + " " + lifaName(lf);
            Time t{yc.y, (uint8_t)yc.mo, (uint8_t)yc.d, 12, 0, 0};
            Chart c = cast(JD::toJD(t), lf);
            if (!c.ok) { CHECK(false, "排盘失败"); continue; }
            CHECK(c.sz[0] == yc.g_ && c.sz[1] == yc.z_,
                  "年柱应为 " + gz(yc.g_, yc.z_) + " 实为 " + gz(c.sz[0], c.sz[1]));
        }
    }

    // ── E5: 农历输入端到端 —— 用农历日期排盘, 必须与其对应公历日期排盘完全一致 ──
    //  使用者可以按农历报生日; 库先把农历折成公历再排盘, 两条入口不允许有分歧。
    //  覆盖普通月、闰月, 以及历史上月建/重月特殊的年份。
    for (int y : {1900, 1984, 2000, 2004, 2023, 2024, 2033, 762, 239, 1500})
    {
        int runMonth = sxtwl::getRunMonth(y);
        struct LunarCase { int m; bool leap; };
        std::vector<LunarCase> lcs = {{1, false}, {6, false}, {12, false}};
        if (runMonth > 0) lcs.push_back({runMonth, true});

        for (const auto &lc : lcs)
        {
            for (int day : {1, 15, 28})
            {
                Day *dp = nullptr;
                try { dp = sxtwl::fromLunar(y, (uint8_t)lc.m, day, lc.leap, false); }
                catch (...) { continue; }
                if (!dp) continue;
                int sy = dp->getSolarYear(), sm = dp->getSolarMonth(), sd = dp->getSolarDay();
                delete dp;

                for (LiFaType lf : {YuWuWeiZiPingLifa_DingDongZhi,
                                    YuWuWeiZiPingLifa_DingXiaZhi,
                                    XianDaiNongLifa_DingQiFa})
                {
                    g.ctx = "农历 " + std::to_string(y) + "-" + (lc.leap ? "闰" : "") +
                            std::to_string(lc.m) + "-" + std::to_string(day) +
                            " → 公历 " + std::to_string(sy) + "-" + std::to_string(sm) +
                            "-" + std::to_string(sd) + " " + lifaName(lf);

                    SBaziInputPara pl;
                    pl.name = "T"; pl.gender = false; pl.isAst = false;
                    pl.jw = {120, 39.9, "t", "t"};
                    pl.lifa = lf; pl.calendar = CalendarLunar;
                    pl.birthDayTime = Time(y, (uint8_t)lc.m, day, 9, 20, 0);
                    pl.isRun = lc.leap; pl.isSpec = false;

                    std::array<int, 8> viaLunar{};
                    bool ok = true;
                    try { BaziBase b(pl); b.calcBaziPaiPan(); viaLunar = b.getSiZhuIndex(); }
                    catch (...) { ok = false; }
                    if (!ok) { CHECK(false, "农历入口排盘抛异常"); continue; }

                    Time st{sy, (uint8_t)sm, (uint8_t)sd, 9, 20, 0};
                    Chart viaSolar = cast(JD::toJD(st), lf);
                    if (!viaSolar.ok) { CHECK(false, "公历入口排盘失败"); continue; }

                    CHECK(viaLunar == viaSolar.sz,
                          "农历入口与公历入口四柱不一致: " +
                          gz(viaLunar[0], viaLunar[1]) + " " + gz(viaLunar[2], viaLunar[3]) + " " +
                          gz(viaLunar[4], viaLunar[5]) + " " + gz(viaLunar[6], viaLunar[7]) + " vs " +
                          gz(viaSolar.sz[0], viaSolar.sz[1]) + " " + gz(viaSolar.sz[2], viaSolar.sz[3]) + " " +
                          gz(viaSolar.sz[4], viaSolar.sz[5]) + " " + gz(viaSolar.sz[6], viaSolar.sz[7]));
                }
            }
        }
    }

    // ── E4: 三种历法在"远离节令"的日子必须给出相同的年月柱 ──
    //  平气与定气对同一个节的时刻最多差两天出头, 故在节中段(月中)应当一致。
    //  这是使用者最直观的期待: 只有卡在交节那一两天才会出现分歧。
    for (int y : {1900, 1984, 2000, 2024, 2100})
    {
        for (int mo = 1; mo <= 12; ++mo)
        {
            g.ctx = std::to_string(y) + "-" + std::to_string(mo) + "-18 三历法比对";
            Time t{y, (uint8_t)mo, 18, 12, 0, 0};
            long double jd = JD::toJD(t);
            Chart a = cast(jd, YuWuWeiZiPingLifa_DingDongZhi);
            Chart b = cast(jd, YuWuWeiZiPingLifa_DingXiaZhi);
            Chart c = cast(jd, XianDaiNongLifa_DingQiFa);
            if (!a.ok || !b.ok || !c.ok) { CHECK(false, "排盘失败"); continue; }
            CHECK(a.sz[3] == c.sz[3] && b.sz[3] == c.sz[3],
                  "月支分歧: 定冬至=" + std::string(kZhiStr[a.sz[3]]) +
                  " 定夏至=" + kZhiStr[b.sz[3]] + " 定气=" + kZhiStr[c.sz[3]]);
            CHECK(a.sz[0] == c.sz[0] && a.sz[1] == c.sz[1] &&
                  b.sz[0] == c.sz[0] && b.sz[1] == c.sz[1],
                  "年柱分歧: 定冬至=" + gz(a.sz[0], a.sz[1]) +
                  " 定夏至=" + gz(b.sz[0], b.sz[1]) + " 定气=" + gz(c.sz[0], c.sz[1]));
        }
    }
}

}  // namespace

// ═══════════════════════════════════════════════════════════════════════
int main(int argc, char **argv)
{
    bool verbose = false;
    int onlyY = 0, rs = 0, re = 0;
    std::string groups = "ABCDEG";

    for (int i = 1; i < argc; ++i)
    {
        std::string s = argv[i];
        if (s == "-v") verbose = true;
        else if (s == "-y" && i + 1 < argc) onlyY = std::atoi(argv[++i]);
        else if (s == "-range" && i + 2 < argc) { rs = std::atoi(argv[++i]); re = std::atoi(argv[++i]); }
        else if (s == "-g" && i + 1 < argc) groups = argv[++i];
    }
    auto want = [&](char c) { return groups.find(c) != std::string::npos; };

    std::vector<int> years;
    if (onlyY) years = {onlyY};
    else if (rs || re) { if (rs > re) std::swap(rs, re); for (int y = rs; y <= re; ++y) years.push_back(y); }
    else
        years = {
            // 上古 / 儒略历时代 —— 至点在公历上的位置漂移最大, 跨年/跨月分档最易出错
            -700, -221, -104, 1, 200, 500, 900, 1200, 1400, 1500,
            // 格里高利改历前后
            1581, 1582, 1583,
            // 历法段切换 / 现代
            1645, 1700, 1800, 1900, 1959, 1960, 1984, 2000, 2024, 2025, 2026,
            // 未来 —— 线性外推误差累积
            2100, 2500, 3000, 5000, 9000,
        };

    std::time_t t0 = std::time(nullptr);

    if (want('A') || want('B') || want('C'))
    {
        for (LiFaType lf : {YuWuWeiZiPingLifa_DingDongZhi, YuWuWeiZiPingLifa_DingXiaZhi})
        {
            std::printf("\n==== %s ====\n", lifaName(lf).c_str());
            for (int Y : years)
            {
                if (verbose) std::printf("  ---- Y=%d ----\n", Y);
                if (want('A') || want('B')) groupA_B(Y, lf, verbose);
                if (want('C')) groupC(Y, lf, verbose);
            }
        }
        if (want('C'))
        {
            std::printf("\n==== 定气 (真太阳时组) ====\n");
            for (int Y : years) groupC(Y, XianDaiNongLifa_DingQiFa, verbose);
        }
    }

    if (want('D'))
    {
        std::printf("\n==== SSQ 单例无副作用 ====\n");
        for (int Y : years) groupD(Y, verbose);
    }

    if (want('E'))
    {
        std::printf("\n==== 外部锚点 ====\n");
        groupE(verbose);
    }

    if (want('G'))
    {
        std::printf("\n==== 报告文本自洽 ====\n");
        for (LiFaType lf : {YuWuWeiZiPingLifa_DingDongZhi, YuWuWeiZiPingLifa_DingXiaZhi,
                            XianDaiNongLifa_DingQiFa})
            for (int Y : years) groupG(Y, lf, verbose);
    }

    std::time_t t1 = std::time(nullptr);
    std::printf("\n───────── 汇总 ─────────\n");
    std::printf("检查 %d 项, 通过 %d, 失败 %d, 用时 %ld 秒\n",
                g.total, g.passed, (int)g.failures.size(), (long)(t1 - t0));
    if (!g.failures.empty())
    {
        std::printf("\n失败明细:\n");
        for (const auto &f : g.failures) std::printf("  ✗ %s\n", f.c_str());
    }
    return g.failures.empty() ? 0 : 1;
}
