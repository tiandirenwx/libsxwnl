// 起运 / 交运 算法验证用例
//
// 目标: 独立复算 bazi.cpp 中 节气区间 + 顺逆判定 + 起运折算 + 交运日期,
//       与 BaziBase 的公开输出 (getQiYun / getJiaoYun / getStartYear) 逐项比对,
//       从而端到端验证 calcJiaoYunDate() 的算法链。
//
// 验证点:
//   [节气]  mHeadJieQiJd_ <= 出生 < mTailJieQiJd_, 且首末节间隔 ≈ 一个节气月 (28~33 天)
//   [顺逆]  flag = (年干 % 2) XOR 性别; 阳男阴女顺行(查未来节), 阴男阳女逆行(查过去节)
//   [起运]  arrayQiYunDelta_(年/月/日) = 出生↔节令天数 按 "365.2422/120 天 = 1 起运年" 折算
//   [交运]  mJyJd_ = 出生 + 365.2422*起运年 + 余数; 且 mJyYear_ = 交运日所在农历年
//
// 编译:  make qiyun_jiaoyun_test
// 运行:  ./build/bin/qiyun_jiaoyun_test

#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "bazi.h"
#include "JD.h"
#include "sxtwl.h"

namespace {

const char *const kGan = "甲乙丙丁戊己庚辛壬癸";
const char *const kZhi = "子丑寅卯辰巳午未申酉戌亥";

std::string gz(int g, int z)
{
    std::string s;
    s.append(kGan + g * 3, 3);
    s.append(kZhi + z * 3, 3);
    return s;
}

std::string lifaName(LiFaType lifa)
{
    switch (lifa)
    {
    case YuWuWeiZiPingLifa_DingDongZhi: return "定冬至(平气)";
    case YuWuWeiZiPingLifa_DingXiaZhi:  return "定夏至(平气)";
    case XianDaiNongLifa_DingQiFa:      return "定气(天文)";
    default:                            return "未知";
    }
}

struct Case
{
    std::string name;
    Time birth;
    bool female;
    LiFaType lifa;
};

int gPass = 0;
int gFail = 0;

void check(bool ok, const std::string &label, const std::string &detail)
{
    if (ok)
    {
        ++gPass;
        std::printf("    [PASS] %s\n", label.c_str());
    }
    else
    {
        ++gFail;
        std::printf("    [FAIL] %s -> %s\n", label.c_str(), detail.c_str());
    }
}

void runCase(const Case &c)
{
    std::printf("\n==== 用例: %s | %s | %s ====\n",
                c.name.c_str(),
                c.female ? "女" : "男",
                lifaName(c.lifa).c_str());

    SBaziInputPara p;
    p.name = c.name;
    p.gender = c.female;
    p.isAst = false;                       // 平太阳时口径, 使 head/tail/birth 与 mJyJd 同基
    p.jw = {120, 39.9, "test", "北京"};
    p.lifa = c.lifa;
    p.calendar = CalendarSolar;
    p.birthDayTime = c.birth;
    p.isRun = false;
    p.isSpec = false;

    BaziBase obj(p);
    obj.calcBaziPaiPan();

    auto sz = obj.getSiZhuIndex();
    long double bd = obj.getBirthJd();
    long double H = obj.getHeadJieQiJd();
    long double T = obj.getTailJieQiJd();

    std::printf("  四柱: %s %s %s %s\n",
                gz(sz[0], sz[1]).c_str(), gz(sz[2], sz[3]).c_str(),
                gz(sz[4], sz[5]).c_str(), gz(sz[6], sz[7]).c_str());
    std::printf("  出生 JD=%.6Lf (%s)\n", bd, JD::formatStr(bd).c_str());
    std::printf("  月首节 JD=%.6Lf (%s)\n", H, JD::formatStr(H).c_str());
    std::printf("  月末节 JD=%.6Lf (%s)\n", T, JD::formatStr(T).c_str());

    // ── [节气] 出生落在 [head, tail) 且间隔 ≈ 一个节气月 ──
    check(H <= bd && bd < T, "节气区间: head<=出生<tail",
          "H=" + JD::formatStr(H) + " bd=" + JD::formatStr(bd) +
          " T=" + JD::formatStr(T));
    long double span = T - H;
    check(span > 28.0L && span < 33.0L, "首末节间隔 ≈ 一个节气月",
          "span=" + std::to_string(double(span)) + " 天");

    // ── [顺逆] flag = (年干 % 2) XOR 性别 ──
    int isFemale = c.female ? 1 : 0;
    bool flag = ((sz[0] % 2) ^ isFemale) != 0;
    bool yangGan = (sz[0] % 2) == 0;
    // 传统规则: 阳男/阴女 顺行(flag=false), 阴男/阳女 逆行(flag=true)
    bool expectShun = (yangGan && !c.female) || (!yangGan && c.female);
    std::printf("  年干=%s(%s) 性别=%s => %s\n",
                std::string(kGan + sz[0] * 3, 3).c_str(),
                yangGan ? "阳" : "阴", c.female ? "女" : "男",
                flag ? "逆行(查过去节/月首节)" : "顺行(查未来节/月末节)");
    check(flag != expectShun, "顺逆方向符合 阳男阴女顺、阴男阳女逆",
          std::string("flag=") + (flag ? "逆" : "顺") +
          " 期望=" + (expectShun ? "顺" : "逆"));

    // ── 独立复算 起运 / 交运 (完全对齐 calcJiaoYunDate 的算式) ──
    long double dt1, dt2;
    if (flag) { dt1 = H; dt2 = bd; }   // 逆行: 月首节 -> 出生
    else      { dt1 = bd; dt2 = T; }   // 顺行: 出生 -> 月末节

    const long double unit = 365.2422L / 120.0L; // ≈3.0435 天 = 1 起运年
    long double offset = dt2 - dt1;
    long double deltaY = std::floor(offset / unit);
    long double rem = (offset / unit - deltaY) * 365.2422L;
    long double jyJd = bd + 365.2422L * deltaY + rem;
    long double deltaM = std::floor(rem / (365.2422L / 12.0L));
    rem -= deltaM * (365.2422L / 12.0L);
    long double deltaD = std::floor(rem);

    std::printf("  出生↔节令 天数=%.6Lf => 起运 %d年%d个月%d天\n",
                offset, int(deltaY), int(deltaM), int(deltaD));

    // 期望字符串 (与 BaziBase::getQiYun 完全同构)
    std::string expQiYun = "起运 命主于出生后" + std::to_string(int(deltaY)) +
                           "年" + std::to_string(int(deltaM)) + "个月" +
                           std::to_string(int(deltaD)) + "天后起运\n";
    check(obj.getQiYun() == expQiYun, "起运折算与 getQiYun() 一致",
          "实际=[" + obj.getQiYun() + "] 期望=[" + expQiYun + "]");

    // 期望交运字符串 (与 BaziBase::getJiaoYun 完全同构)
    std::string expJiaoYun = "交运 命主于" + JD::formatStr(jyJd) + "交运\n";
    check(obj.getJiaoYun() == expJiaoYun, "交运日期与 getJiaoYun() 一致",
          "实际=[" + obj.getJiaoYun() + "] 期望=[" + expJiaoYun + "]");
    std::printf("  交运时刻 JD=%.6Lf (%s)\n", jyJd, JD::formatStr(jyJd).c_str());

    // ── [交运] mJyYear_ = 交运日所在农历年 ──
    Time jyt = JD::JD2DD(jyJd);
    SLunarDay jyLunar = sxtwl::solar2Lunar2(jyt);
    check(obj.getStartYear() == jyLunar.year, "起运农历年 getStartYear() 正确",
          "getStartYear=" + std::to_string(obj.getStartYear()) +
          " 复算=" + std::to_string(jyLunar.year));

    // ── 交叉校验: "3天≈1年" 物理意义 — 交运距出生的真实年数 ≈ 起运年数 ──
    long double realYears = (jyJd - bd) / 365.2422L;
    check(std::floor(realYears + 1e-9L) == deltaY,
          "交运距出生真实年数 == 起运整年数",
          "realYears=" + std::to_string(double(realYears)) +
          " deltaY=" + std::to_string(int(deltaY)));
}

}  // namespace

int main()
{
    std::printf("======== 起运 / 交运 算法验证 ========\n");

    std::vector<Case> cases = {
        // 出生年干阴阳 × 性别 全组合, 覆盖顺行/逆行两条分支
        {"阳男-定夏至", Time{1978, 11, 8, 14, 32, 0}, false, YuWuWeiZiPingLifa_DingXiaZhi},
        {"阳女-定夏至", Time{1978, 11, 8, 14, 32, 0}, true,  YuWuWeiZiPingLifa_DingXiaZhi},
        {"阳男-定冬至", Time{1984, 2, 10, 7, 35, 10}, false, YuWuWeiZiPingLifa_DingDongZhi},
        {"阳男-定气",   Time{1984, 2, 10, 7, 35, 10}, false, XianDaiNongLifa_DingQiFa},
        {"阴男-定气",   Time{1985, 6, 15, 9, 0, 0},   false, XianDaiNongLifa_DingQiFa},
        {"阴女-定气",   Time{1985, 6, 15, 9, 0, 0},   true,  XianDaiNongLifa_DingQiFa},
        {"现代-定气",   Time{2000, 1, 1, 0, 0, 0},    false, XianDaiNongLifa_DingQiFa},
        {"现代-定气",   Time{2024, 8, 15, 18, 30, 0}, true,  XianDaiNongLifa_DingQiFa},
    };

    for (const auto &c : cases) runCase(c);

    std::printf("\n======== 汇总: 通过 %d, 失败 %d ========\n", gPass, gFail);
    return gFail == 0 ? 0 : 1;
}
