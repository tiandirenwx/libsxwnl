#pragma once
#include <array>
#include <string>
#include <vector>

// ═══════════════════════════════════════════════════════════════
//  bazi_analysis —— 八字命理派生信息计算（纯函数, 无外部依赖）
//
//  设计目标:
//   - 仅以整数索引为输入输出, 不依赖天文历法/JD 等重模块,
//     便于在 iOS / Android / Web 等平台直接复用。
//   - 所有规则均来自《渊海子平》《三命通会》等公开命理典籍,
//     与 lunar-java(6tail) 等成熟开源项目一致, 不做主观臆造。
//
//  索引约定:
//   - 天干 gan ∈ [0,9]   甲乙丙丁戊己庚辛壬癸
//   - 地支 zhi ∈ [0,11]  子丑寅卯辰巳午未申酉戌亥
//   - 五行 wx  ∈ [0,4]   木火土金水
//
//  四柱数组 sz[8] 约定:
//   [0]年干 [1]年支 [2]月干 [3]月支 [4]日干 [5]日支 [6]时干 [7]时支
// ═══════════════════════════════════════════════════════════════

// 名称表(天干Gan/地支Zhi/十神ShiShen/纳音NaYinWuXing/藏干gCangGan等)
// 统一复用 sx_lang_zh.h, 本模块不再重复定义。

namespace bazi {

// ─── 基础属性 ───
int  ganWuXing(int gan);        // 天干 → 五行
int  zhiWuXing(int zhi);        // 地支(本气) → 五行
bool ganIsYang(int gan);        // 天干阴阳, true=阳
int  jiaZiIndex(int gan, int zhi); // 干支 → 六十甲子序号 [0,59], 非法返回 -1

// ─── 派生信息 ───
std::string naYin(int gan, int zhi);                 // 纳音, 例 "海中金"
// 十神索引表: out[gan] = ShiShen 下标, 与 BaziBase::vecShiShen_ 建表规则一致
void buildShiShenMap(int dayGan, int out[10]);
int  shiShenIndex(int dayGan, int targetGan);         // 查表, 非法返回 -1

std::string shiShen(int dayGan, int targetGan);      // 十神全称
std::string shiShenShort(int dayGan, int targetGan); // 十神简称
std::vector<int> cangGan(int zhi);                   // 地支藏干(顺序同 gCangGan)
int benQiGan(int zhi);                               // 地支本气天干索引
std::string changSheng(int gan, int zhi);            // 十二长生, 例 "帝旺"
std::array<int, 2> kongWang(int gan, int zhi);       // 旬空两地支索引
std::string kongWangStr(int gan, int zhi);           // 旬空, 例 "寅卯"

// 五行统计(木火土金水). includeCangGan=true 含藏干(含藏气)
std::array<int, 5> wuXingCount(const std::array<int, 8> &sz, bool includeCangGan);

// 五行旺衰: 以月支定令, 返回木火土金水各自状态("旺""相""休""囚""死")
std::array<std::string, 5> wuXingStatus(int monthZhi);

// 人元司令分野: 月支 + 距本月节令天数 → 当令天干索引
int siLing(int monthZhi, int daysAfterJie);

// 神煞: 基于命局 sz 上下文, 判断目标柱(targetGan,targetZhi)所带神煞
std::vector<std::string> shenSha(const std::array<int, 8> &sz,
                                 int targetGan, int targetZhi);

} // namespace bazi
