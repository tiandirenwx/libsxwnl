// 简易验证: 对若干样本日, 输出神煞/宜忌/吉时/用事
//
// 用法 (在 repo 根目录):
//   c++ -std=c++17 -Isrc -Lbuild/lib scripts/verify_almanac.cpp -lsxwnl -o /tmp/verify_alm
//   /tmp/verify_alm
#include "sxwnl_capi.h"
#include <cstdio>
#include <cstring>
#include <vector>

static const char *kZhi[12] = {"子","丑","寅","卯","辰","巳","午","未","申","酉","戌","亥"};

static void show(int y, int m, int d, const char* label) {
    SxwnlAlmanac a{};
    SxwnlDayInfo info{};
    if (sxwnl_get_day_info(y, m, d, &info) != 0 ||
        sxwnl_get_almanac (y, m, d, &a)    != 0) {
        printf("[%s %04d-%02d-%02d] fetch failed\n", label, y, m, d);
        return;
    }
    printf("════════════════════════════════════\n");
    printf("[%s] %04d-%02d-%02d  %s%s %s日柱 (建除:%s)\n",
        label, y, m, d, info.lunar_month_name, info.lunar_day_name,
        info.day_gz, info.jian12_name);
    printf("二十八宿: %s | 十二神: %s | 冲: %s | 煞: %s\n",
        a.xiu, a.twelve_god, a.chong_sheng_xiao, a.sha);

    printf("\n[神煞] 共 %d:\n", a.shen_sha_count);
    for (int i = 0; i < a.shen_sha_count; ++i)
        printf("  %s %s (w%d)\n",
            a.shen_sha[i].is_lucky ? "[吉]" : "[凶]",
            a.shen_sha[i].name, a.shen_sha[i].weight);

    printf("[宜]");
    for (int i = 0; i < a.yi_count; ++i) printf(" %s", a.yi[i]);
    printf("\n[忌]");
    for (int i = 0; i < a.ji_count; ++i) printf(" %s", a.ji[i]);
    printf("\n");

    printf("[吉时]");
    for (int i = 0; i < a.lucky_hour_count; ++i) {
        printf(" %s:%s时", a.lucky_hours[i].name, kZhi[a.lucky_hours[i].zhi]);
    }
    printf("\n");

    printf("[用事择吉]\n");
    for (int i = 0; i < a.event_count; ++i)
        printf("  %s:%s (%s)\n",
            a.events[i].event,
            a.events[i].suitable ? "宜" : "忌",
            a.events[i].reason);

    printf("[提示]\n");
    for (int i = 0; i < a.note_count; ++i) printf("  - %s\n", a.notes[i]);

    if (a.quote_count > 0) {
        printf("[董公诀] %s | %s\n",
            a.quotes[0].title, a.quotes[0].luck);
    }
    printf("\n");
}

static void show_topics() {
    SxwnlAlmanacTopic items[32];
    int n = sxwnl_get_almanac_topics(items, 32);
    printf("════ topics: %d ════\n", n);
    for (int i = 0; i < n; ++i) {
        printf("  [%s] %s\n    body_len=%zu\n    %s\n",
            items[i].category, items[i].title,
            std::strlen(items[i].body), items[i].body);
    }
}

int main() {
    show_topics();
    // 关键样本:
    // 1. 2025-02-03 立春 (节气切换)
    // 2. 2025-02-12 杨公忌(正月十三, 农历)
    // 3. 2025-03-22 春分前一日 (四离)
    // 4. 2025-04-04 清明
    // 5. 2025-04-21 谷雨, 农历三月廿四
    // 6. 2024-09-22 秋分前一日 (四离)
    // 7. 2025-01-29 春节 (乙巳年大年初一)
    show(2025, 2,  3,  "立春");
    show(2025, 2, 10,  "正月十三 杨公忌");
    show(2025, 3, 19,  "春分前一日 四离");
    show(2025, 4, 30,  "三月初三");
    show(2024, 9, 22,  "秋分前一日");
    show(2025, 1, 29,  "乙巳大年初一");
    show(2025, 7,  7,  "甲午天赦?");
    show(2025, 6, 20,  "夏至前一日 四离");
    return 0;
}
