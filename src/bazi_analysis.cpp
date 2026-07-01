#include "bazi_analysis.h"
#include "sx_lang_zh.h" // 复用 Gan/Zhi/ShiShen/NaYinWuXing/gCangGan 等定义
#include <algorithm>

namespace bazi {
namespace {

// 天干五行(木火土金水): 甲乙木 丙丁火 戊己土 庚辛金 壬癸水
const int kGanWuXing[10] = {0, 0, 1, 1, 2, 2, 3, 3, 4, 4};
// 地支本气五行: 子水 丑土 寅木 卯木 辰土 巳火 午火 未土 申金 酉金 戌土 亥水
const int kZhiWuXing[12] = {4, 2, 0, 0, 2, 1, 1, 2, 3, 3, 2, 4};

// 十神简称(顺序对应 sx_lang_zh.h 中 ShiShen 全称)
const char *const kShiShenShort[10] = {"比", "劫", "食", "伤", "才",
                                       "财", "杀", "官", "枭", "印"};

// 十二长生名
const char *const kChangSheng[12] = {"长生", "沐浴", "冠带", "临官", "帝旺", "衰",
                                     "病", "死", "墓", "绝", "胎", "养"};
// 各天干长生所在地支
const int kChangShengStart[10] = {11, 6, 2, 9, 2, 9, 5, 0, 8, 3};

// 三合局分组: 0=申子辰(水) 1=寅午戌(火) 2=巳酉丑(金) 3=亥卯未(木)
int sanHeGroup(int zhi) {
    switch (zhi) {
    case 8: case 0: case 4:  return 0;
    case 2: case 6: case 10: return 1;
    case 5: case 9: case 1:  return 2;
    case 11: case 3: case 7: return 3;
    }
    return -1;
}

// 人元司令分野表: 每月(地支)的当令天干及其天数(顺序即先后)
struct SiLingSeg { int gan; int days; };
const std::vector<SiLingSeg> kSiLing[12] = {
    /*子*/ {{8, 10}, {9, 20}},
    /*丑*/ {{9, 9}, {7, 3}, {5, 18}},
    /*寅*/ {{4, 7}, {2, 7}, {0, 16}},
    /*卯*/ {{0, 10}, {1, 20}},
    /*辰*/ {{1, 9}, {9, 3}, {4, 18}},
    /*巳*/ {{4, 5}, {6, 9}, {2, 16}},
    /*午*/ {{2, 10}, {5, 9}, {3, 11}},
    /*未*/ {{3, 9}, {1, 3}, {5, 18}},
    /*申*/ {{4, 10}, {8, 3}, {6, 17}},
    /*酉*/ {{6, 10}, {7, 20}},
    /*戌*/ {{7, 9}, {3, 3}, {4, 18}},
    /*亥*/ {{4, 7}, {0, 5}, {8, 18}},
};

// 天乙贵人(日干→两地支)
const int kTianYi[10][2] = {
    {1, 7}, {0, 8}, {11, 9}, {11, 9}, {1, 7},
    {0, 8}, {1, 7}, {6, 2}, {3, 5}, {3, 5}};
// 文昌贵人(日干→地支)
const int kWenChang[10] = {5, 6, 8, 9, 8, 9, 11, 0, 2, 3};
// 禄神(日干→地支)
const int kLuShen[10] = {2, 3, 5, 6, 5, 6, 8, 9, 11, 0};
// 羊刃(日干→地支)
const int kYangRen[10] = {3, 4, 6, 7, 6, 7, 9, 10, 0, 1};
// 金舆(日干→地支)
const int kJinYu[10] = {4, 5, 7, 8, 7, 8, 10, 11, 1, 2};

// 三合局→地支(顺序对应分组 0水 1火 2金 3木)
const int kTaoHua[4]   = {9, 3, 6, 0};   // 桃花(咸池)
const int kYiMa[4]     = {2, 8, 11, 5};  // 驿马
const int kHuaGai[4]   = {4, 10, 1, 7};  // 华盖
const int kJiangXing[4]= {0, 6, 9, 3};   // 将星
const int kJieSha[4]   = {5, 11, 2, 8};  // 劫煞
const int kZaiSha[4]   = {6, 0, 3, 9};   // 灾煞
const int kWangShen[4] = {11, 5, 8, 2};  // 亡神
const int kPanAn[4]    = {1, 7, 10, 4};  // 攀鞍

// 孤辰/寡宿(年支四孟分组→孤辰、寡宿地支)
void guChenGuaSu(int yearZhi, int &gu, int &gua) {
    // 亥子丑→寅戌 寅卯辰→巳丑 巳午未→申辰 申酉戌→亥未
    if (yearZhi == 11 || yearZhi == 0 || yearZhi == 1) { gu = 2;  gua = 10; }
    else if (yearZhi == 2 || yearZhi == 3 || yearZhi == 4) { gu = 5;  gua = 1; }
    else if (yearZhi == 5 || yearZhi == 6 || yearZhi == 7) { gu = 8;  gua = 4; }
    else { gu = 11; gua = 7; }
}

// 月德贵人(月支三合局→天干)
int yueDeGan(int monthZhi) {
    switch (sanHeGroup(monthZhi)) {
    case 1: return 2; // 寅午戌→丙
    case 0: return 8; // 申子辰→壬
    case 3: return 0; // 亥卯未→甲
    case 2: return 6; // 巳酉丑→庚
    }
    return -1;
}

// 天德贵人(月支→目标). type:0天干 1地支; idx<0 表示无(子月为巽,八字略)
// 《三命通会》: 寅丁卯申辰壬巳辛午亥未甲申癸酉寅戌丙亥乙子(巽)丑庚
struct DeTarget { int type; int idx; };
const DeTarget kTianDe[12] = {
    {-1, -1}, // 子(巽,略)
    {0, 6},   // 丑→庚
    {0, 3},   // 寅→丁
    {1, 8},   // 卯→申
    {0, 8},   // 辰→壬
    {0, 7},   // 巳→辛
    {1, 11},  // 午→亥
    {0, 0},   // 未→甲
    {0, 9},   // 申→癸
    {1, 2},   // 酉→寅
    {0, 2},   // 戌→丙
    {0, 1},   // 亥→乙
};

// 天干五合: 合(g) = (g+5)%10  (甲己 乙庚 丙辛 丁壬 戊癸)
inline int ganHe(int g) { return (g + 5) % 10; }

void appendUnique(std::vector<std::string> &v, const std::string &s) {
    if (std::find(v.begin(), v.end(), s) == v.end()) v.push_back(s);
}

} // namespace

// ─── 基础属性 ───
int ganWuXing(int gan) { return (gan >= 0 && gan < 10) ? kGanWuXing[gan] : -1; }
int zhiWuXing(int zhi) { return (zhi >= 0 && zhi < 12) ? kZhiWuXing[zhi] : -1; }
bool ganIsYang(int gan) { return (gan % 2) == 0; }

int jiaZiIndex(int gan, int zhi) {
    if (gan < 0 || gan > 9 || zhi < 0 || zhi > 11) return -1;
    for (int n = 0; n < 60; n++)
        if (n % 10 == gan && n % 12 == zhi) return n;
    return -1;
}

// ─── 纳音 ───
std::string naYin(int gan, int zhi) {
    int n = jiaZiIndex(gan, zhi);
    if (n < 0) return "";
    return NaYinWuXing[n / 2];
}

void buildShiShenMap(int dayGan, int out[10]) {
    for (int i = 0; i < 10; i++) out[i] = -1;
    if (dayGan < 0 || dayGan > 9) return;
    int j = dayGan;
    if (j % 2 == 0) {
        for (int i = 0; i < 10; i++) {
            out[j] = i;
            j = (j + 1) % 10;
        }
    } else {
        for (int k = 0; k < 9; k += 2) {
            out[j] = k;
            out[j - 1] = k + 1;
            j = (j + 2) % 10;
        }
    }
}

int shiShenIndex(int dayGan, int targetGan) {
    if (dayGan < 0 || dayGan > 9 || targetGan < 0 || targetGan > 9) return -1;
    int map[10];
    buildShiShenMap(dayGan, map);
    return map[targetGan];
}

// ─── 十神 ───
std::string shiShen(int dayGan, int targetGan) {
    int idx = shiShenIndex(dayGan, targetGan);
    if (idx < 0 || idx >= 10) return "";
    return ShiShen[idx];
}
std::string shiShenShort(int dayGan, int targetGan) {
    int idx = shiShenIndex(dayGan, targetGan);
    if (idx < 0 || idx >= 10) return "";
    return kShiShenShort[idx];
}

// ─── 地支藏干(本气在前, 顺序同 gCangGan) ───
std::vector<int> cangGan(int zhi) {
    std::vector<int> out;
    if (zhi < 0 || zhi > 11) return out;
    for (int i = 0; i < 3; i++)
        if (gCangGan[zhi][i] < 10) out.push_back(gCangGan[zhi][i]);
    return out;
}

// 本气天干: gCangGan 已将本气置于首位
int benQiGan(int zhi) {
    return (zhi >= 0 && zhi <= 11) ? gCangGan[zhi][0] : -1;
}

// ─── 十二长生 ───
std::string changSheng(int gan, int zhi) {
    if (gan < 0 || gan > 9 || zhi < 0 || zhi > 11) return "";
    int start = kChangShengStart[gan];
    int idx = ganIsYang(gan) ? ((zhi - start + 12) % 12)
                             : ((start - zhi + 12) % 12);
    return kChangSheng[idx];
}

// ─── 旬空(空亡) ───
std::array<int, 2> kongWang(int gan, int zhi) {
    if (gan < 0 || zhi < 0) return {-1, -1};
    int a = (zhi - gan + 10) % 12;
    int b = (zhi - gan + 11) % 12;
    return {a, b};
}
std::string kongWangStr(int gan, int zhi) {
    auto kw = kongWang(gan, zhi);
    if (kw[0] < 0) return "";
    return std::string(Zhi[kw[0]]) + Zhi[kw[1]];
}

// ─── 五行统计 ───
std::array<int, 5> wuXingCount(const std::array<int, 8> &sz, bool includeCangGan) {
    std::array<int, 5> cnt{0, 0, 0, 0, 0};
    for (int i = 0; i < 8; i += 2) { // 四天干
        int wx = ganWuXing(sz[i]);
        if (wx >= 0) cnt[wx]++;
    }
    for (int i = 1; i < 8; i += 2) { // 四地支
        int zhi = sz[i];
        if (includeCangGan) {
            for (int g : cangGan(zhi)) {
                int wx = ganWuXing(g);
                if (wx >= 0) cnt[wx]++;
            }
        } else {
            int wx = zhiWuXing(zhi);
            if (wx >= 0) cnt[wx]++;
        }
    }
    return cnt;
}

// ─── 五行旺衰 ───
std::array<std::string, 5> wuXingStatus(int monthZhi) {
    // 令: 旺; 令生者: 相; 生令者: 休; 令克者: 死; 克令者: 囚
    std::array<std::string, 5> st;
    int ling = zhiWuXing(monthZhi); // 月支本气五行即当令五行
    for (int e = 0; e < 5; e++) {
        if (e == ling) st[e] = "旺";
        else if (e == (ling + 1) % 5) st[e] = "相";  // 令生者相
        else if (e == (ling + 4) % 5) st[e] = "休";  // 生令者休
        else if (e == (ling + 2) % 5) st[e] = "死";  // 令克者死
        else st[e] = "囚";                            // 克令者囚
    }
    return st;
}

// ─── 人元司令 ───
int siLing(int monthZhi, int daysAfterJie) {
    if (monthZhi < 0 || monthZhi > 11) return -1;
    if (daysAfterJie < 0) daysAfterJie = 0;
    const auto &segs = kSiLing[monthZhi];
    int acc = 0;
    for (const auto &s : segs) {
        acc += s.days;
        if (daysAfterJie < acc) return s.gan;
    }
    return segs.empty() ? -1 : segs.back().gan;
}

// ─── 神煞 ───
std::vector<std::string> shenSha(const std::array<int, 8> &sz,
                                 int targetGan, int targetZhi) {
    std::vector<std::string> out;
    int yearGan = sz[0], yearZhi = sz[1];
    int monthZhi = sz[3], dayGan = sz[4], dayZhi = sz[5];

    // —— 以日干为基(命中目标地支) ——
    if (targetZhi >= 0 && dayGan >= 0) {
        if (targetZhi == kTianYi[dayGan][0] || targetZhi == kTianYi[dayGan][1])
            appendUnique(out, "天乙贵人");
        if (targetZhi == kWenChang[dayGan]) appendUnique(out, "文昌贵人");
        if (targetZhi == kLuShen[dayGan])   appendUnique(out, "禄神");
        if (targetZhi == kYangRen[dayGan])  appendUnique(out, "羊刃");
        if (targetZhi == kJinYu[dayGan])    appendUnique(out, "金舆");
    }

    // —— 以年支、日支三合局为基 ——
    auto checkSanHe = [&](int refZhi) {
        int g = sanHeGroup(refZhi);
        if (g < 0 || targetZhi < 0) return;
        if (targetZhi == kTaoHua[g])    appendUnique(out, "桃花");
        if (targetZhi == kYiMa[g])      appendUnique(out, "驿马");
        if (targetZhi == kHuaGai[g])    appendUnique(out, "华盖");
        if (targetZhi == kJiangXing[g]) appendUnique(out, "将星");
        if (targetZhi == kJieSha[g])    appendUnique(out, "劫煞");
        if (targetZhi == kZaiSha[g])    appendUnique(out, "灾煞");
        if (targetZhi == kWangShen[g])  appendUnique(out, "亡神");
        if (targetZhi == kPanAn[g])     appendUnique(out, "攀鞍");
    };
    checkSanHe(yearZhi);
    checkSanHe(dayZhi);

    // —— 孤辰寡宿(以年支为基) ——
    if (targetZhi >= 0 && yearZhi >= 0) {
        int gu, gua;
        guChenGuaSu(yearZhi, gu, gua);
        if (targetZhi == gu)  appendUnique(out, "孤辰");
        if (targetZhi == gua) appendUnique(out, "寡宿");
    }

    // —— 月德贵人 / 月德合(目标天干) ——
    int md = yueDeGan(monthZhi);
    if (md >= 0) {
        if (targetGan == md) appendUnique(out, "月德贵人");
        if (targetGan == ganHe(md)) appendUnique(out, "月德合");
    }

    // —— 天德贵人 / 天德合 ——
    DeTarget td = kTianDe[monthZhi];
    if (td.idx >= 0) {
        if (td.type == 0) { // 天德为天干
            if (targetGan == td.idx) appendUnique(out, "天德贵人");
            if (targetGan == ganHe(td.idx)) appendUnique(out, "天德合");
        } else { // 天德为地支
            if (targetZhi == td.idx) appendUnique(out, "天德贵人");
        }
    }

    (void)yearGan;
    return out;
}

} // namespace bazi
