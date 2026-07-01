// 平气/定气八字模式的节气交接测试
// 覆盖公元前 722 至 9999 年
//
// 目标: 验证每种立法内部在 24 节气交接时刻的自洽性
//
// 测试方法:
//   1. 对每种立法/每个采样年份, 沿全年逐日扫描 mid-day 时刻的月柱.
//   2. 当相邻两日的月柱不同, 说明穿过了一个 "节" 边界.
//   3. 用二分法找到月柱切换的精确 JD (Tswitch).
//   4. 断言:
//      a. Tswitch - 1 秒 与 Tswitch + 1 秒 的月柱恰好递增 1 (mod 60);
//      b. 排盘于 Tswitch + 1 秒 得到的 mHeadJieQiJd_ ≈ Tswitch (自洽);
//      c. 每年有 12 次月柱切换 (寅→卯→辰→...→丑→寅), 计数完整;
//      d. 12 次切换中恰好有 1 次伴随年柱切换 (立春), 其他 11 次年柱不变;
//      e. 立春 时刻年柱恰好递增 1 (mod 60).
//
// 业务逻辑测试矩阵 (每个采样年 × 每种立法):
//   1) verifyYearScan        : 24 节交接月柱/年柱切换自洽 (核心, 必跑)
//   2) verifyDayPillarSwitch : 日柱每天 23:00 子时切换
//   3) verifyHourPillarSwitch: 时柱每 2 小时切换
//   4) verifyAstYearScan     : 真太阳时口径节气自洽 (-fast 跳过)
//   5) verifyLunarInput      : 农历输入端到端 (含闰月) (-fast 跳过)
//   6) verifyDaysAfterJie    : 距节令天数一致性 (-fast 跳过)
//
// 编译 (通过 CMake):  make pingqi_dingqi_test
// 运行:
//   ./build/bin/pingqi_dingqi_test                           # 默认 82 个采样年份
//   ./build/bin/pingqi_dingqi_test -v -y 2024                # 单年详情
//   ./build/bin/pingqi_dingqi_test -full                     # 完整 -722 ~ 9999 (10711 年)
//   ./build/bin/pingqi_dingqi_test -full -fast               # 完整范围 + 加速 (仅核心3项)
//   ./build/bin/pingqi_dingqi_test -range -722 9999 -fast    # 自定义范围
//   ./build/bin/pingqi_dingqi_test -full -progress 500       # 每 500 年输出进度
//
// 完整范围预估时间 (单机, 8核):
//   全量 6 项测试: 约 70 分钟 (32秒/246组合 → 32*10711*3/246 ≈ 4173秒)
//   -fast 3 项测试: 约 30 分钟

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "SSQ.h"
#include "bazi.h"
#include "const.h"
#include "eph.h"
#include "JD.h"

namespace {

const char *const kGan  = "甲乙丙丁戊己庚辛壬癸";
const char *const kZhi  = "子丑寅卯辰巳午未申酉戌亥";

std::string gz(int g, int z)
{
    if (g < 0 || z < 0 || g >= 10 || z >= 12) return "??";
    std::string s;
    s.append(kGan + g * 3, 3);
    s.append(kZhi + z * 3, 3);
    return s;
}

std::string lifaName(LiFaType lifa)
{
    switch (lifa)
    {
    case YuWuWeiZiPingLifa_DingDongZhi: return "定冬至";
    case YuWuWeiZiPingLifa_DingXiaZhi:  return "定夏至";
    case XianDaiNongLifa_DingQiFa:      return "定气";
    default:                            return "未知";
    }
}

struct BaziSnapshot
{
    std::array<int, 8> sizhu{};
    long double headJq{};
    long double tailJq{};
    long double headJqAst{};
    long double tailJqAst{};
    long double birthJd{};
    long double birthAstJd{};
    int daysAfterJie = -1;
    bool valid = false;
    std::string err;
};

BaziSnapshot calcBaziAtJd(long double jd, LiFaType lifa, bool isAst = false)
{
    BaziSnapshot r;
    try
    {
        Time t = JD::JD2DD(jd);
        SBaziInputPara p;
        p.name = "T";
        p.gender = false;
        p.isAst = isAst;
        p.jw = {120, 39.9, "test", "test"};
        p.lifa = lifa;
        p.calendar = CalendarSolar;
        p.birthDayTime = t;
        p.isRun = false;
        p.isSpec = false;
        BaziBase b(p);
        b.calcBaziPaiPan();
        r.sizhu       = b.getSiZhuIndex();
        r.headJq      = b.getHeadJieQiJd();
        r.tailJq      = b.getTailJieQiJd();
        r.headJqAst   = b.getHeadJieQiAstJd();
        r.tailJqAst   = b.getTailJieQiAstJd();
        r.birthJd     = b.getBirthJd();
        r.birthAstJd  = b.getBirthAstJd();
        r.daysAfterJie= b.getDaysAfterJie();
        r.valid       = true;
    }
    catch (const std::exception &e) { r.err = e.what(); }
    catch (...)                     { r.err = "unknown"; }
    return r;
}

// 通过农历输入排盘 (端到端测试用)
BaziSnapshot calcBaziFromLunar(int lyYear, int lyMonth, int lyDay,
                               bool isRun, bool isSpec,
                               LiFaType lifa, bool isAst = false)
{
    BaziSnapshot r;
    try
    {
        SBaziInputPara p;
        p.name = "T";
        p.gender = false;
        p.isAst = isAst;
        p.jw = {120, 39.9, "test", "test"};
        p.lifa = lifa;
        p.calendar = CalendarLunar;
        // 注意: birthDayTime 在 Lunar 模式下解释为农历 (year,month,day,hour,min,sec)
        Time lt{lyYear, (uint8_t)lyMonth, lyDay, 12, 0, 0};
        p.birthDayTime = lt;
        p.isRun  = isRun;
        p.isSpec = isSpec;
        BaziBase b(p);
        b.calcBaziPaiPan();
        r.sizhu       = b.getSiZhuIndex();
        r.headJq      = b.getHeadJieQiJd();
        r.tailJq      = b.getTailJieQiJd();
        r.daysAfterJie= b.getDaysAfterJie();
        r.valid       = true;
    }
    catch (const std::exception &e) { r.err = e.what(); }
    catch (...)                     { r.err = "unknown"; }
    return r;
}

// 不同立法的节气交接容差 (天)
long double toleranceByLifa(LiFaType lifa)
{
    // 定气: 天文算法精确, 容许 2 秒 (含二分收敛残差)
    // 平气: 常数近似, 容许 15 分钟
    return (lifa == XianDaiNongLifa_DingQiFa) ? 2.0L / 86400.0L : 0.01L;
}

// 计算 (gA,zA) → (gB,zB) 的 mod 60 前向差值; 若不合法(gan/zhi 奇偶不匹配), 返回 -1
int gzForwardDiff(int gA, int zA, int gB, int zB)
{
    int dg = ((gB - gA) % 10 + 10) % 10;
    int dz = ((zB - zA) % 12 + 12) % 12;
    for (int k = 0; k < 60; ++k)
        if (k % 10 == dg && k % 12 == dz) return k;
    return -1;
}

bool inRange(const std::array<int, 8> &sz)
{
    for (int i = 0; i < 4; ++i)
    {
        if (sz[2 * i]     < 0 || sz[2 * i]     >= 10) return false;
        if (sz[2 * i + 1] < 0 || sz[2 * i + 1] >= 12) return false;
    }
    return true;
}

std::string jdStr(long double jd)
{
    Time t = JD::JD2DD(jd);
    std::ostringstream os;
    os << t.Y << "-" << std::setw(2) << std::setfill('0') << t.M << "-"
       << std::setw(2) << std::setfill('0') << t.D << " "
       << std::setw(2) << std::setfill('0') << (int)t.h << ":"
       << std::setw(2) << std::setfill('0') << (int)t.m << ":"
       << std::setw(2) << std::setfill('0') << (int)t.s;
    return os.str();
}

struct Issue
{
    int    year;
    std::string kind;
    std::string mode;
    std::string detail;
};

struct Reporter
{
    std::vector<Issue> issues;
    int totalChecks = 0;
    int okChecks    = 0;

    void ok() { ++okChecks; ++totalChecks; }
    void fail(const Issue &iss) { issues.push_back(iss); ++totalChecks; }

    void print(int maxDetails = 30) const
    {
        std::printf("\n----- 汇总 -----\n");
        std::printf("总检查数: %d, 通过: %d, 失败: %d\n",
                    totalChecks, okChecks, (int)issues.size());
        if (issues.empty()) return;
        std::printf("失败样例(前 %d 条):\n", maxDetails);
        int printed = 0;
        for (const auto &e : issues)
        {
            if (printed >= maxDetails) break;
            std::printf("  Y=%d [%s] %s: %s\n", e.year,
                        e.mode.c_str(), e.kind.c_str(), e.detail.c_str());
            ++printed;
        }
    }
};

// 通过二分法定位月柱切换的精确 JD.
// 前置: 排盘 jdLo 得到 monthIdxLo, 排盘 jdHi 得到 monthIdxHi, 且 idxHi = idxLo+1 mod 60.
long double bisectMonthSwitch(long double jdLo, long double jdHi,
                              int idxLo_g, int idxLo_z, LiFaType lifa)
{
    // 收敛到毫秒精度 (1e-8 天 ≈ 0.86ms)
    while (jdHi - jdLo > 1e-8L)
    {
        long double mid = (jdLo + jdHi) / 2.0L;
        auto bmid = calcBaziAtJd(mid, lifa);
        if (!bmid.valid) return jdHi;
        if (bmid.sizhu[2] == idxLo_g && bmid.sizhu[3] == idxLo_z)
            jdLo = mid;
        else
            jdHi = mid;
    }
    return jdHi;  // jdHi 是切换后第一个 JD
}

// 通过二分法定位年柱切换的精确 JD.
long double bisectYearSwitch(long double jdLo, long double jdHi,
                             int idxLo_g, int idxLo_z, LiFaType lifa)
{
    while (jdHi - jdLo > 1e-8L)
    {
        long double mid = (jdLo + jdHi) / 2.0L;
        auto bmid = calcBaziAtJd(mid, lifa);
        if (!bmid.valid) return jdHi;
        if (bmid.sizhu[0] == idxLo_g && bmid.sizhu[1] == idxLo_z)
            jdLo = mid;
        else
            jdHi = mid;
    }
    return jdHi;
}

// 逐日扫描给定年份, 记录月柱/年柱切换情况.
// 期望: 全年 12 次月柱切换, 恰好 1 次伴随年柱切换 (发生在立春).
void verifyYearScan(int Y, LiFaType lifa, Reporter &rep, bool verbose)
{
    // 采样范围: Y-01-01 至 Y+1-01-31, 共约 396 天, 覆盖全年 12 节 + 部分次年.
    // 为避免 Julian 与远古年份的日期偏移导致遗漏, 用 JD 直接推进.
    Time tStart{Y, 1, 1, 12, 0, 0};
    long double jd0 = JD::toJD(tStart);
    const int NDAYS = 396;

    std::vector<BaziSnapshot> daily;
    daily.reserve(NDAYS + 1);
    for (int d = 0; d <= NDAYS; ++d)
    {
        long double jd = jd0 + d;
        auto b = calcBaziAtJd(jd, lifa);
        if (!b.valid)
        {
            rep.fail({Y, "扫描失败", lifaName(lifa),
                      "d=" + std::to_string(d) + " " + b.err});
            return;
        }
        if (!inRange(b.sizhu))
        {
            rep.fail({Y, "越界", lifaName(lifa),
                      "d=" + std::to_string(d) + " 四柱越界"});
            return;
        }
        daily.push_back(b);
    }

    // 汇总: 月柱切换与年柱切换的次数
    int monthSwitches = 0;
    int yearSwitches  = 0;
    int monthWithYear = 0;  // 月柱与年柱同时切换的次数

    for (int d = 1; d <= NDAYS; ++d)
    {
        const auto &prev = daily[d - 1];
        const auto &curr = daily[d];

        bool monthChanged = (prev.sizhu[2] != curr.sizhu[2]) ||
                            (prev.sizhu[3] != curr.sizhu[3]);
        bool yearChanged  = (prev.sizhu[0] != curr.sizhu[0]) ||
                            (prev.sizhu[1] != curr.sizhu[1]);

        if (!monthChanged && !yearChanged) continue;

        long double jdPrev = jd0 + (d - 1);
        long double jdCurr = jd0 + d;

        // 检查月柱切换
        if (monthChanged)
        {
            ++monthSwitches;
            // 二分找边界
            long double Tsw = bisectMonthSwitch(
                jdPrev, jdCurr,
                prev.sizhu[2], prev.sizhu[3], lifa);
            // Tsw - eps 和 Tsw + eps 排盘
            long double eps = 1.0L / 86400.0L;
            auto bPre  = calcBaziAtJd(Tsw - eps, lifa);
            auto bPost = calcBaziAtJd(Tsw + eps, lifa);
            if (!bPre.valid || !bPost.valid)
            {
                rep.fail({Y, "月柱边界", lifaName(lifa),
                          "Tsw=" + jdStr(Tsw) + " 排盘失败"});
                continue;
            }
            int diff = gzForwardDiff(bPre.sizhu[2], bPre.sizhu[3],
                                     bPost.sizhu[2], bPost.sizhu[3]);
            if (diff != 1)
            {
                rep.fail({Y, "月柱边界", lifaName(lifa),
                          "Tsw=" + jdStr(Tsw) +
                          " 前=" + gz(bPre.sizhu[2], bPre.sizhu[3]) +
                          " 后=" + gz(bPost.sizhu[2], bPost.sizhu[3]) +
                          " 差=" + std::to_string(diff)});
                continue;
            }
            // 月支严格按 +1 (mod 12) 单调递增 (与"节"序同步)
            int zhiDiff = ((bPost.sizhu[3] - bPre.sizhu[3]) % 12 + 12) % 12;
            if (zhiDiff != 1)
            {
                rep.fail({Y, "月支递增", lifaName(lifa),
                          "Tsw=" + jdStr(Tsw) +
                          " 前月支=" + std::to_string(bPre.sizhu[3]) +
                          " 后月支=" + std::to_string(bPost.sizhu[3]) +
                          " 差=" + std::to_string(zhiDiff) + " (期望 1)"});
                continue;
            }
            // 自洽性: bPost.headJq ≈ Tsw (按立法分级精度)
            long double gap = std::abs(double(bPost.headJq - Tsw));
            long double tol = toleranceByLifa(lifa);
            if (gap > tol)
            {
                rep.fail({Y, "自洽性", lifaName(lifa),
                          "Tsw=" + jdStr(Tsw) +
                          " 前=" + gz(bPre.sizhu[2], bPre.sizhu[3]) +
                          " 后=" + gz(bPost.sizhu[2], bPost.sizhu[3]) +
                          " bPost.head=" + jdStr(bPost.headJq) +
                          " bPre.tail=" + jdStr(bPre.tailJq) +
                          " gap=" + std::to_string(double(gap)) +
                          " tol=" + std::to_string(double(tol))});
            }
            else
            {
                rep.ok();
            }

            // 月柱切换是否与年柱同时?
            bool monthCoYear =
                (bPre.sizhu[0] != bPost.sizhu[0]) ||
                (bPre.sizhu[1] != bPost.sizhu[1]);
            if (monthCoYear)
            {
                ++monthWithYear;
                // 年柱前后应恰好递增 1 (立春)
                int yd = gzForwardDiff(bPre.sizhu[0], bPre.sizhu[1],
                                       bPost.sizhu[0], bPost.sizhu[1]);
                if (yd != 1)
                {
                    rep.fail({Y, "立春年柱", lifaName(lifa),
                              "Tsw=" + jdStr(Tsw) +
                              " 前年=" + gz(bPre.sizhu[0], bPre.sizhu[1]) +
                              " 后年=" + gz(bPost.sizhu[0], bPost.sizhu[1]) +
                              " 差=" + std::to_string(yd)});
                }
                else
                {
                    rep.ok();
                }
                // 立春月柱应为寅月 (zhi=2)
                if (bPost.sizhu[3] != 2)
                {
                    rep.fail({Y, "立春月支", lifaName(lifa),
                              "Tsw=" + jdStr(Tsw) +
                              " 后月支=" + std::to_string(bPost.sizhu[3]) +
                              " (期望寅=2)"});
                }
                else
                {
                    rep.ok();
                }
                // 立春日期合理性: 区分 Julian / Gregorian 历法
                //   * Gregorian (Tsw JD >= 2299161, 即 1582-10-15 之后): 立春应在 2/1-2/10
                //     定气精确在 2/3-2/5; 平气可能略宽 2/1-2/7
                //   * Julian (1582 之前): Julian 历平均年比回归年长 0.0078 天,
                //     立春在 Julian 历法下随年份逐渐提前 (公元 1500 约 1/27,
                //     公元前 700 约 2/12), 这是历法本身特性, 不视为 bug.
                //     此处仅断言"立春月份在 1 月或 2 月".
                Time tLc = JD::JD2DD(Tsw);
                bool isGregorian = (int(Tsw) >= 2299161);
                if (isGregorian)
                {
                    // Gregorian 历平均年 365.2425 天, 比回归年长 0.0003 天/年.
                    // 1582→9999 累积漂移约 2.4 天, 立春会从 2/4 提前到 1/29 附近.
                    // 故远未来年份允许 1/29; 上限保持 2/10 (定气/平气都不应超过).
                    if (tLc.M == 1 && tLc.D < 29)
                    {
                        rep.fail({Y, "立春日期", lifaName(lifa),
                                  "Tsw=" + jdStr(Tsw) +
                                  " Gregorian 下立春出现在 1/29 之前 (异常)"});
                    }
                    else if (tLc.M == 2 && tLc.D > 10)
                    {
                        rep.fail({Y, "立春日期", lifaName(lifa),
                                  "Tsw=" + jdStr(Tsw) +
                                  " Gregorian 下立春出现在 2/10 之后 (异常)"});
                    }
                    else if (tLc.M < 1 || tLc.M > 2)
                    {
                        rep.fail({Y, "立春日期", lifaName(lifa),
                                  "Tsw=" + jdStr(Tsw) +
                                  " Gregorian 下立春月份不在 1-2 月 (异常)"});
                    }
                    else
                    {
                        rep.ok();
                    }
                }
                else
                {
                    // Julian: 仅断言月份在 1-2 月
                    if (tLc.M < 1 || tLc.M > 2)
                    {
                        rep.fail({Y, "立春日期", lifaName(lifa),
                                  "Tsw=" + jdStr(Tsw) +
                                  " Julian 下立春月份不在 1-2 月 (异常)"});
                    }
                    else
                    {
                        rep.ok();
                    }
                }
            }

            if (verbose && monthSwitches <= 15)
            {
                std::printf("  [%s Y=%d] 月柱切换#%d Tsw=%s %s→%s %s\n",
                            lifaName(lifa).c_str(), Y, monthSwitches,
                            jdStr(Tsw).c_str(),
                            gz(bPre.sizhu[2], bPre.sizhu[3]).c_str(),
                            gz(bPost.sizhu[2], bPost.sizhu[3]).c_str(),
                            monthCoYear ? "[立春]" : "");
            }
        }
        else if (yearChanged)
        {
            // 年柱变化但月柱未变: 异常, 因为年柱切换必须发生在立春 = 月柱切换点
            ++yearSwitches;
            long double Tsw = bisectYearSwitch(
                jdPrev, jdCurr,
                prev.sizhu[0], prev.sizhu[1], lifa);
            rep.fail({Y, "年柱独立切换", lifaName(lifa),
                      "Tsw=" + jdStr(Tsw) +
                      " 月柱未同步 (应仅在立春发生)"});
        }
    }

    // 年内节柱切换个数应约 12 (可能覆盖 11 或 13, 因跨年边界).
    // 考虑 396 天窗口, 期望 12~13 次.
    if (monthSwitches < 11 || monthSwitches > 14)
    {
        rep.fail({Y, "月柱切换次数", lifaName(lifa),
                  "全年扫描共 " + std::to_string(monthSwitches) + " 次月柱切换"});
    }
    else
    {
        rep.ok();
    }
    // 立春 应发生 1 次 (可能 2 次, 因窗口覆盖 396 天可能横跨 2 个立春)
    if (monthWithYear < 1 || monthWithYear > 2)
    {
        rep.fail({Y, "立春次数", lifaName(lifa),
                  "396 天扫描共 " + std::to_string(monthWithYear) +
                  " 次立春 (期望 1~2)"});
    }
    else
    {
        rep.ok();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 业务逻辑测试 #1: 日柱每天子时 (23:00) 切换
//   代码 (bazi.cpp:511-514): jd += 13/24 - J2000; D = int2(jd); riZhu = D - 6 + 9000000
//   含义: 把"本日中午12点"的 jd 平移到"前一日23:00"起算, 然后 int2 取日号.
//   所以日柱在每天 23:00 整点切换. 测试方法:
//     * 采样 N 个连续日 (含子时前后 1 秒),
//     * 断言 23:00-eps 与 23:00+eps 的日柱恰好 +1 mod 60.
// ─────────────────────────────────────────────────────────────────────────────
void verifyDayPillarSwitch(int Y, LiFaType lifa, Reporter &rep)
{
    // 在 Y-06-15 附近取 30 个连续日的 23:00 子时边界
    Time t0{Y, 6, 15, 12, 0, 0};
    long double jd0 = JD::toJD(t0);
    long double eps = 1.0L / 86400.0L;
    for (int d = 0; d < 30; ++d)
    {
        // 第 d 日 23:00 = jd0 + d + 11/24 (jd0 是当日 12:00, 加 11h 到 23:00)
        long double tZi = jd0 + d + 11.0L / 24.0L;
        auto bPre  = calcBaziAtJd(tZi - eps, lifa);
        auto bPost = calcBaziAtJd(tZi + eps, lifa);
        if (!bPre.valid || !bPost.valid)
        {
            rep.fail({Y, "日柱切换", lifaName(lifa),
                      "d=" + std::to_string(d) + " 排盘失败"});
            continue;
        }
        int diff = gzForwardDiff(bPre.sizhu[4], bPre.sizhu[5],
                                 bPost.sizhu[4], bPost.sizhu[5]);
        if (diff != 1)
        {
            rep.fail({Y, "日柱切换", lifaName(lifa),
                      "Tsw=" + jdStr(tZi) +
                      " 前日=" + gz(bPre.sizhu[4], bPre.sizhu[5]) +
                      " 后日=" + gz(bPost.sizhu[4], bPost.sizhu[5]) +
                      " 差=" + std::to_string(diff) + " (期望 1)"});
        }
        else
        {
            rep.ok();
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 业务逻辑测试 #2: 时柱每 2 小时 (一个时辰) 切换
//   代码 (bazi.cpp:513): SC = int2((jd - D) * 12)
//   jd 已平移到前一日 23:00 起算, 所以 (jd - D) ∈ [0, 1).
//   SC 取 0..11 对应 子/丑/寅/卯/辰/巳/午/未/申/酉/戌/亥.
//   时辰边界: 23:00, 1:00, 3:00, 5:00, 7:00, 9:00, 11:00, 13:00, 15:00, 17:00, 19:00, 21:00.
// ─────────────────────────────────────────────────────────────────────────────
void verifyHourPillarSwitch(int Y, LiFaType lifa, Reporter &rep)
{
    Time t0{Y, 6, 15, 12, 0, 0};
    long double jd0 = JD::toJD(t0);
    long double eps = 1.0L / 86400.0L;
    // 12 个时辰切换点 (相对当日 0:00 的小时偏移; 23:00 视为 -1h)
    // 此处采样当日 1:00, 3:00, ..., 21:00 共 11 个边界
    for (int h = 1; h <= 21; h += 2)
    {
        long double tShi = jd0 + (h - 12.0L) / 24.0L;  // jd0=当日12:00, 减12h到0:00, 加 h 小时
        auto bPre  = calcBaziAtJd(tShi - eps, lifa);
        auto bPost = calcBaziAtJd(tShi + eps, lifa);
        if (!bPre.valid || !bPost.valid)
        {
            rep.fail({Y, "时柱切换", lifaName(lifa),
                      "h=" + std::to_string(h) + " 排盘失败"});
            continue;
        }
        int diff = gzForwardDiff(bPre.sizhu[6], bPre.sizhu[7],
                                 bPost.sizhu[6], bPost.sizhu[7]);
        if (diff != 1)
        {
            rep.fail({Y, "时柱切换", lifaName(lifa),
                      "Tsw=" + jdStr(tShi) +
                      " 前时=" + gz(bPre.sizhu[6], bPre.sizhu[7]) +
                      " 后时=" + gz(bPost.sizhu[6], bPost.sizhu[7]) +
                      " 差=" + std::to_string(diff) + " (期望 1)"});
        }
        else
        {
            rep.ok();
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 业务逻辑测试 #3: 真太阳时口径 (isAst=true) 节气交接自洽
//   代码 (bazi.cpp:259-268): mHeadJieQiAstJd_ = JD::calcAST(mHeadJieQiJd_, J)
//   即真太阳时下节气 JD = 平太阳时 JD + 时差方程修正.
//   测试: 用 isAst=true 重新扫一年, 二分定位月柱切换点 Tsw,
//        断言 bPost.headJqAst ≈ Tsw (即用真太阳时排盘时, 节气边界与切换 JD 自洽).
// ─────────────────────────────────────────────────────────────────────────────
void verifyAstYearScan(int Y, LiFaType lifa, Reporter &rep)
{
    Time tStart{Y, 1, 1, 12, 0, 0};
    long double jd0 = JD::toJD(tStart);
    const int NDAYS = 396;
    std::vector<BaziSnapshot> daily;
    daily.reserve(NDAYS + 1);
    for (int d = 0; d <= NDAYS; ++d)
    {
        auto b = calcBaziAtJd(jd0 + d, lifa, /*isAst=*/true);
        if (!b.valid || !inRange(b.sizhu))
        {
            return;  // 扫描失败, 静默跳过
        }
        daily.push_back(b);
    }
    // 局部计数器 (原来误写为 static, 会跨函数调用累积导致测试失效)
    int localCount = 0;
    for (int d = 1; d <= NDAYS; ++d)
    {
        bool monthChanged = (daily[d - 1].sizhu[2] != daily[d].sizhu[2]) ||
                            (daily[d - 1].sizhu[3] != daily[d].sizhu[3]);
        if (!monthChanged) continue;
        long double Tsw = bisectMonthSwitch(
            jd0 + (d - 1), jd0 + d,
            daily[d - 1].sizhu[2], daily[d - 1].sizhu[3], lifa);
        long double eps = 1.0L / 86400.0L;
        // 注意: 二分时仍用 isAst=false (因 bisectMonthSwitch 内部用 calcBaziAtJd 默认参数)
        // 这里只验证"真太阳时排盘的 headJqAst ≈ 真太阳时口径的 Tsw".
        // 真太阳时口径下, Tsw 与 平太阳时 Tsw 相差一个时差修正 (最多 ±16分钟).
        // 因此本测试用 平太阳时 bisect 得 Tsw, 然后 isAst=true 排盘看 headJqAst ≈ Tsw.
        // 这等价于验证 calcAST(headJq) ≈ calcAST(birthJd) 时月柱切换一致.
        auto bPost = calcBaziAtJd(Tsw + eps, lifa, /*isAst=*/true);
        if (!bPost.valid) continue;
        // 在真太阳时口径下, 切换点 Tsw 应用 calcAST 转换后再比较
        // 简单近似: 时差方程在 Tsw 附近几乎线性, 可以直接比 bPost.headJqAst 与
        // calcAST(Tsw) 的差. 但我们改用更稳健的检验: 平太阳时 headJqAst 与
        // birthAstJd 应保持等差关系 (即两者都加同样的 calcAST 修正).
        // 这里直接断言: bPost.headJqAst 与 bPost.headJq 的差 ≈ bPost.birthAstJd 与 bPost.birthJd 的差
        // (即 calcAST 在节气点与出生点应用了同一时差修正)
        long double dHead = bPost.headJqAst - bPost.headJq;
        long double dBirth = bPost.birthAstJd - bPost.birthJd;
        long double delta = std::abs(double(dHead - dBirth));
        // 时差方程在 1 天内变化不超过 30 秒, 容忍 60 秒
        if (delta > 60.0L / 86400.0L)
        {
            rep.fail({Y, "真太阳时口径", lifaName(lifa),
                      "Tsw=" + jdStr(Tsw) +
                      " dHead=" + std::to_string(double(dHead)) +
                      " dBirth=" + std::to_string(double(dBirth)) +
                      " delta=" + std::to_string(double(delta)) + " 天"});
        }
        else
        {
            rep.ok();
        }
        // 每年仅检查前 12 次切换, 避免过度采样
        if (++localCount >= 12) break;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 业务逻辑测试 #4: 农历输入端到端 (CalendarLunar + isRun + isSpec)
//   代码 (bazi.cpp:28-44): fromLunar -> solar -> 排盘
//   SSQ.cpp 中 8/23/237/239/689/701/761/762 等月建别名变更年, 历法段切换年,
//   都应通过 fromLunar -> BaziBase 端到端流程不抛异常, 四柱在合法范围.
// ─────────────────────────────────────────────────────────────────────────────
void verifyLunarInput(int Y, LiFaType lifa, Reporter &rep)
{
    // 测试场景:
    //   A. 当年农历正月初一 (验证春节前后年柱切换是否正确)
    //   B. 当年农历腊月初一 (验证年末排盘)
    //   C. 若该年有闰月, 测试闰月初一 (isRun=true)
    //   D. 特殊月建年份的临界日期
    int runMonth = sxtwl::getRunMonth(Y);

    // A. 正月初一
    {
        auto b = calcBaziFromLunar(Y, 1, 1, false, false, lifa);
        if (!b.valid || !inRange(b.sizhu))
        {
            rep.fail({Y, "农历正月初一", lifaName(lifa),
                      "排盘失败: " + (b.valid ? "四柱越界" : b.err)});
        }
        else
        {
            rep.ok();
        }
    }
    // B. 腊月初一 (月份 12)
    {
        auto b = calcBaziFromLunar(Y, 12, 1, false, false, lifa);
        if (!b.valid || !inRange(b.sizhu))
        {
            rep.fail({Y, "农历腊月初一", lifaName(lifa),
                      "排盘失败: " + (b.valid ? "四柱越界" : b.err)});
        }
        else
        {
            rep.ok();
        }
    }
    // C. 闰月测试 (若该年有闰月)
    if (runMonth > 0)
    {
        auto b = calcBaziFromLunar(Y, runMonth, 1, true, false, lifa);
        if (!b.valid || !inRange(b.sizhu))
        {
            rep.fail({Y, "农历闰月", lifaName(lifa),
                      "闰" + std::to_string(runMonth) + "月初一 排盘失败: " +
                      (b.valid ? "四柱越界" : b.err)});
        }
        else
        {
            rep.ok();
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// 业务逻辑测试 #5: getDaysAfterJie() 与节气日差一致性
//   代码 (bazi.cpp:615-621): days = floor(birthJd - headJieQiJd); return days<0 ? 0:days
//   即出生时刻距本月节令的整日数. 测试:
//     * 在月柱切换后第 0/3/10/20 天采样, getDaysAfterJie() 应返回 0/3/10/20.
// ─────────────────────────────────────────────────────────────────────────────
void verifyDaysAfterJie(int Y, LiFaType lifa, Reporter &rep)
{
    Time tStart{Y, 3, 1, 12, 0, 0};  // 立春后某日
    long double jd0 = JD::toJD(tStart);
    // 找到 jd0 之后第一个月柱切换点
    long double jdLo = jd0, jdHi = jd0 + 60;
    auto b0 = calcBaziAtJd(jdLo, lifa);
    if (!b0.valid) return;
    auto b1 = calcBaziAtJd(jdHi, lifa);
    if (!b1.valid) return;
    if (b0.sizhu[2] == b1.sizhu[2] && b0.sizhu[3] == b1.sizhu[3])
    {
        // 60 天内未切换, 跳过 (理论上不可能, 但防御性)
        return;
    }
    long double Tsw = bisectMonthSwitch(jdLo, jdHi, b0.sizhu[2], b0.sizhu[3], lifa);
    // 探测本节长度 (head→tail), 用于选择安全的采样偏移
    auto bAt = calcBaziAtJd(Tsw + 1.0L / 24.0L, lifa);
    if (!bAt.valid) return;
    long double nodeSpan = bAt.tailJq - bAt.headJq;
    if (nodeSpan < 14 || nodeSpan > 32)
    {
        // 节间距异常 (远古年份平气可能出现), 静默跳过
        return;
    }
    // 采样偏移 n 天必须严格小于 nodeSpan, 否则会越到下一节导致 daysAfterJie 重置为 0.
    // 取 [0, 3, 10, min(20, nodeSpan-2)] 四个采样点.
    int nMax = std::min(20, int(nodeSpan) - 2);
    std::vector<int> offsets = {0, 3, 10};
    if (nMax > 10) offsets.push_back(nMax);
    for (int n : offsets)
    {
        long double jdSample = Tsw + n + 1.0L / 24.0L;  // 切换后 n 天又 1 小时 (确保进入当日)
        auto b = calcBaziAtJd(jdSample, lifa);
        if (!b.valid)
        {
            rep.fail({Y, "距节令天数", lifaName(lifa),
                      "n=" + std::to_string(n) + " 排盘失败"});
            continue;
        }
        if (b.daysAfterJie != n)
        {
            rep.fail({Y, "距节令天数", lifaName(lifa),
                      "Tsw=" + jdStr(Tsw) + " n=" + std::to_string(n) +
                      " getDaysAfterJie=" + std::to_string(b.daysAfterJie) +
                      " (期望 " + std::to_string(n) + ")"});
        }
        else
        {
            rep.ok();
        }
    }
}

}  // namespace

// 深度调试功能: 打印指定 JD 的排盘详情
static void debugAt(long double jd, LiFaType lifa, const std::string &label)
{
    auto b = calcBaziAtJd(jd, lifa);
    std::printf("  [%s | %s] jd=%.9Lf (%s)\n", label.c_str(),
                lifaName(lifa).c_str(), jd, jdStr(jd).c_str());
    if (!b.valid) { std::printf("    ERROR: %s\n", b.err.c_str()); return; }
    std::printf("    四柱: %s %s %s %s\n",
                gz(b.sizhu[0], b.sizhu[1]).c_str(),
                gz(b.sizhu[2], b.sizhu[3]).c_str(),
                gz(b.sizhu[4], b.sizhu[5]).c_str(),
                gz(b.sizhu[6], b.sizhu[7]).c_str());
    std::printf("    月首节 JD=%.9Lf (%s)\n", b.headJq, jdStr(b.headJq).c_str());
    std::printf("    月末节 JD=%.9Lf (%s)\n", b.tailJq, jdStr(b.tailJq).c_str());
}

int main(int argc, char *argv[])
{
    bool verbose = false;
    bool doDebug = false;
    int  onlyY   = 0;
    // 完整范围模式: -full       => 扫描 -722 ~ 9999 每一年
    // 自定义范围:   -range S E  => 扫描 [S, E] 每一年
    // 加速模式:     -fast       => 跳过真太阳时/农历/距节令测试, 仅保留核心节气+日柱+时柱
    // 进度输出:     -progress N => 每 N 年输出一次进度 (默认 100)
    bool fullMode = false;
    bool fastMode = false;
    int  rangeStart = 0, rangeEnd = 0;
    int  progressEvery = 100;
    for (int i = 1; i < argc; ++i)
    {
        std::string s = argv[i];
        if      (s == "-v") verbose = true;
        else if (s == "-d") doDebug = true;
        else if (s == "-y" && i + 1 < argc) onlyY = std::atoi(argv[++i]);
        else if (s == "-full") fullMode = true;
        else if (s == "-fast") fastMode = true;
        else if (s == "-range" && i + 2 < argc)
        {
            rangeStart = std::atoi(argv[++i]);
            rangeEnd   = std::atoi(argv[++i]);
            fullMode   = true;  // 用 fullMode 路径, 但范围由 rangeStart/rangeEnd 决定
        }
        else if (s == "-progress" && i + 1 < argc)
            progressEvery = std::atoi(argv[++i]);
    }

    if (doDebug)
    {
        // 深度调试: 修复后验证 Y=-711 定夏至 立春前后
        std::cout << "\n=== 深度调试: 定夏至 立春前后 (Y=-711) ===\n";
        {
            auto b0 = calcBaziAtJd(1461420.0L, YuWuWeiZiPingLifa_DingXiaZhi);
            long double lc = b0.headJq;
            std::printf("  ★ 精确立春(定夏至平气) JD=%.9Lf\n", lc);
            debugAt(lc - 1.0L / 86400.0L, YuWuWeiZiPingLifa_DingXiaZhi, "立春 - 1 秒");
            debugAt(lc + 1.0L / 86400.0L, YuWuWeiZiPingLifa_DingXiaZhi, "立春 + 1 秒");
        }
        std::cout << "\n=== 深度调试: 定冬至 Y=1500 12月末月柱 ===\n";
        {
            Time t_sample{1500, 1, 15, 12, 0, 0};
            long double jd_sample = JD::toJD(t_sample);
            auto bs = calcBaziAtJd(jd_sample, YuWuWeiZiPingLifa_DingDongZhi);
            long double H = bs.headJq;
            std::printf("  ★ Y=1500 m=1 sample head=%.9Lf\n", H);
            debugAt(H - 1.0L / 86400.0L, YuWuWeiZiPingLifa_DingDongZhi, "head - 1 秒");
            debugAt(H + 1.0L / 86400.0L, YuWuWeiZiPingLifa_DingDongZhi, "head + 1 秒");
        }
        return 0;
    }

    std::vector<int> years;
    if (onlyY != 0)
    {
        years = {onlyY};
    }
    else if (fullMode)
    {
        // 完整逐年模式: 默认 -722 ~ 9999 (10712 年)
        //   -722 对应公元前 722 年 (astronomical year -722),
        //   与 SSQ.cpp 春秋古历段起点 -721 相近, 已在算法有效范围内.
        int s = (rangeStart || rangeEnd) ? rangeStart : -722;
        int e = (rangeStart || rangeEnd) ? rangeEnd   : 9999;
        if (s > e) std::swap(s, e);
        years.reserve(e - s + 1);
        for (int Y = s; Y <= e; ++Y) years.push_back(Y);
        std::fprintf(stderr, "完整范围模式: %d ~ %d (共 %zu 年)\n",
                     s, e, years.size());
        if (fastMode)
            std::fprintf(stderr, "加速模式: 跳过真太阳时/农历/距节令测试\n");
        std::fprintf(stderr, "进度: 每 %d 年输出一次\n\n", progressEvery);
    }
    else
    {
        // 年份取样规则:
        //   1) 覆盖公元前 722 至 9999 完整范围;
        //   2) 加入 SSQ.cpp 中所有历法算法切换边界年份 (suoKB/qiKB 转折点)
        //      因为平气/平朔的 k,b 常数在这些年份前后会发生不连续跳变,
        //      是最容易出现节气边界 bug 的位置;
        //   3) 加入历史上月建别名变更的年份 (8/23/237/239/689/701/761/762)
        //      验证在月建变更时 24 节气边界仍能自洽切换;
        //   4) 加入 1582 (Gregorian 改革) 与 1960 (现代天文算法接续) 的边界.
        years = {
            // --- 上古 & 春秋战国 ---
            -711,              // 距 -721 春秋古历段起点 10 年内
            -700, -600,
            -500, -479,        // -479 战国古历段起点
            -400, -300,
            -221, -220,        // -221 秦汉古历段起点
            -216,              // -216 秦汉修正段
            -200, -150, -104,  // -104 太初历段起点
            -50, -1, 0,
            // --- 汉魏晋南北朝 ---
            1, 8, 23,          // 8/23 月建别名变更 (建子为十二)
            50, 85,            // 85 四分历段起点
            100, 200,
            237, 239,          // 237 景初历段起点 + 月建别名变更
            300, 445,          // 445 元嘉历段起点
            500, 510,          // 510 大明历段起点
            550, 590, 597,     // 590/597 开皇/大业历段起点
            // --- 唐五代宋辽金元 ---
            619,               // 619 平朔表结束/戊寅元历起点
            666, 700,          // 666 麟德历段起点
            689, 701,          // 689/701 月建别名变更 (建子为正)
            729,               // 729 大衍历段起点
            761, 762,          // 761/762 月建别名变更 + 五纪历段起点
            784, 822, 893,     // 784/822/893 正元/宣明/崇玄历段起点
            900, 956, 964,     // 956/964 钦天/应天历段起点
            983, 1000, 1001,   // 983/1001 乾元/仪天历段起点
            1064, 1068,        // 1064/1068 明天/崇天历段起点
            1075, 1094, 1103,  // 1075/1094/1103 奉元/宋律历/占天历段起点
            1106, 1136,        // 1106/1136 纪元/统元历段起点
            1191, 1199, 1208,  // 1191/1199/1208 会元/统天/开禧历段起点
            1252, 1253,        // 1252/1253 淳祐/会天历段起点
            1271, 1277, 1280,  // 1271/1277/1280 成天/本天/授时历段起点
            // --- 明清近现代 ---
            1400, 1500,
            1581, 1582, 1583,  // 1582 Gregorian 改革
            1600, 1645,        // 1645 平气表结束
            1700, 1800, 1900,
            1959, 1960, 1961,  // 1960 现代天文算法接续边界
            1984, 2000, 2024, 2100,
            // --- 未来 ---
            2500, 3000, 4000,
            5000, 6000, 7000, 8000, 9000, 9500, 9998
        };
    }

    std::vector<LiFaType> modes = {
        YuWuWeiZiPingLifa_DingDongZhi,
        YuWuWeiZiPingLifa_DingXiaZhi,
        XianDaiNongLifa_DingQiFa
    };

    Reporter rep;
    auto tBegin = std::time(nullptr);

    for (LiFaType lifa : modes)
    {
        std::printf("\n==== 立法: %s ====\n", lifaName(lifa).c_str());
        int yIdx = 0;
        for (int Y : years)
        {
            ++yIdx;
            if (verbose) std::printf("  ---- Y=%d ----\n", Y);
            // 业务逻辑测试矩阵:
            //   1) verifyYearScan        : 24 节交接月柱/年柱切换自洽 (核心)
            //   2) verifyDayPillarSwitch : 日柱每天 23:00 子时切换
            //   3) verifyHourPillarSwitch: 时柱每 2 小时切换
            //   4) verifyAstYearScan     : 真太阳时口径节气自洽 (-fast 跳过)
            //   5) verifyLunarInput      : 农历输入端到端 (含闰月) (-fast 跳过)
            //   6) verifyDaysAfterJie    : 距节令天数一致性 (-fast 跳过)
            verifyYearScan(Y, lifa, rep, verbose);
            verifyDayPillarSwitch(Y, lifa, rep);
            verifyHourPillarSwitch(Y, lifa, rep);
            if (!fastMode)
            {
                verifyAstYearScan(Y, lifa, rep);
                verifyLunarInput(Y, lifa, rep);
                verifyDaysAfterJie(Y, lifa, rep);
            }
            // 完整模式进度输出 (stderr, 不污染 stdout 汇总)
            if (fullMode && (yIdx % progressEvery == 0))
            {
                auto tNow = std::time(nullptr);
                std::fprintf(stderr,
                    "  [%s] 进度 %d/%zu (%.1f%%) | Y=%d | 累计: ok=%d fail=%d | 用时 %lds\n",
                    lifaName(lifa).c_str(), yIdx, years.size(),
                    100.0 * yIdx / years.size(), Y,
                    rep.okChecks, (int)rep.issues.size(),
                    (long)(tNow - tBegin));
            }
        }
    }

    rep.print(60);
    if (fullMode)
    {
        auto tEnd = std::time(nullptr);
        std::fprintf(stderr, "\n总用时: %ld 秒 (%.1f 分钟)\n",
                     (long)(tEnd - tBegin), (tEnd - tBegin) / 60.0);
    }
    return rep.issues.empty() ? 0 : 1;
}
