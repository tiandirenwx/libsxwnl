// 历谱口径节气 (Day::getLiPuJieQiName) 回归用例
//
// 验证 Day 的两套节气口径:
//   - 天文口径 getJieQiName()   : qi_accurate 精确交气时刻所在日
//   - 历谱口径 getLiPuJieQiName(): SSQ 整日表 + QB 定气修正 (对应权威 sxwnl ob.Ljq)
//
// 关键断言 (对齐权威 sxwnl 网页版):
//   公元 900 年(儒略, 1645 年以前, QB 修正生效): 两套口径相差 1~2 天
//     小寒: 历谱 1/2, 天文 1/1
//     立春: 历谱 2/1, 天文 1/30
//   现代年份(2026): 两套口径完全重合
//
// 编译:  make lipu_jieqi_test
// 运行:  ./build/bin/lipu_jieqi_test

#include <cstdio>
#include <memory>
#include <string>

#include "day.h"
#include "sxtwl.h"

namespace {

int gPass = 0, gFail = 0;

void check(bool ok, const std::string &label, const std::string &detail = "")
{
    if (ok) { ++gPass; std::printf("  [PASS] %s\n", label.c_str()); }
    else    { ++gFail; std::printf("  [FAIL] %s -> %s\n", label.c_str(), detail.c_str()); }
}

// 取指定公历(古代按儒略历)日期的两套节气名
struct JQ { std::string astro, lipu; };
JQ jqAt(int y, int m, int d)
{
    std::unique_ptr<Day> day(sxtwl::fromSolar(y, (uint8_t)m, d));
    return { day->getJieQiName(), day->getLiPuJieQiName() };
}

}  // namespace

int main()
{
    std::printf("======== 历谱口径节气 (ob.Ljq) 回归 ========\n");

    std::printf("\n[公元 900 年 (儒略, QB 修正生效, 两套口径应分离)]\n");
    // 小寒: 天文 1/1, 历谱 1/2
    check(jqAt(900, 1, 1).astro == "小寒", "900-01-01 天文口径 = 小寒", jqAt(900,1,1).astro);
    check(jqAt(900, 1, 1).lipu.empty(), "900-01-01 历谱口径 = 空(非节气日)", "['"+jqAt(900,1,1).lipu+"']");
    check(jqAt(900, 1, 2).lipu == "小寒", "900-01-02 历谱口径 = 小寒", jqAt(900,1,2).lipu);
    check(jqAt(900, 1, 2).astro.empty(), "900-01-02 天文口径 = 空", "['"+jqAt(900,1,2).astro+"']");

    // 立春: 天文 1/30, 历谱 2/1
    check(jqAt(900, 1, 30).astro == "立春", "900-01-30 天文口径 = 立春", jqAt(900,1,30).astro);
    check(jqAt(900, 1, 30).lipu.empty(), "900-01-30 历谱口径 = 空", "['"+jqAt(900,1,30).lipu+"']");
    check(jqAt(900, 2, 1).lipu == "立春", "900-02-01 历谱口径 = 立春", jqAt(900,2,1).lipu);
    check(jqAt(900, 2, 1).astro.empty(), "900-02-01 天文口径 = 空", "['"+jqAt(900,2,1).astro+"']");

    std::printf("\n[现代 2026 年 (1645 年以后, 两套口径应重合)]\n");
    {
        auto mz = jqAt(2026, 6, 5);   // 芒种
        auto xz = jqAt(2026, 6, 21);  // 夏至
        check(mz.astro == "芒种" && mz.lipu == "芒种", "2026-06-05 天文=历谱=芒种",
              "astro='"+mz.astro+"' lipu='"+mz.lipu+"'");
        check(xz.astro == "夏至" && xz.lipu == "夏至", "2026-06-21 天文=历谱=夏至",
              "astro='"+xz.astro+"' lipu='"+xz.lipu+"'");
    }

    std::printf("\n======== 汇总: 通过 %d, 失败 %d ========\n", gPass, gFail);
    return gFail == 0 ? 0 : 1;
}
