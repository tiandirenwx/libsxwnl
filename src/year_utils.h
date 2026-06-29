#pragma once

#include <string>

// ═══════════════════════════════════════════════════════════════
//  year_utils —— 纪年与时间字符串工具
//
//  对应上游 sxwnl/tools.js 中的 year2Ayear / Ayear2year / timeStr2hour
// ═══════════════════════════════════════════════════════════════

namespace year_utils
{

// ── 历史纪年 ↔ 天文学纪年 ─────────────────────────────────────
//
//  天文学纪年: 没有"公元0年", 但用 0 代表公元前 1 年, -1 代表公元前 2 年 …
//  历史纪年:   "公元X年" / "公元前X年", 没有 0
//
//  规则:
//    历史 X (>0)   <->  天文 X
//    公元前 X      <->  天文 (1 - X)

// 将历史纪年字符串解析为天文学纪年 (整数)
// 接受形式: "2024" / "2024年" / "公元2024" / "公元前221" / "前221" / "-221"
// 解析失败返回 INT32_MIN 作为错误标记
int year2Ayear(const std::string &s);

// 将天文学纪年转回字符串
//  fullStyle = true (默认): "公元2024年" / "公元前221年"
//  fullStyle = false       : "2024" / "前221"
std::string Ayear2year(int aYear, bool fullStyle = true);

// ── 时间字符串 → 小时数(浮点) ─────────────────────────────────
//
//  支持:   "HH"           => h
//          "HH:MM"        => h + m/60
//          "HH:MM:SS"     => h + m/60 + s/3600
//          "HH:MM:SS.fff" => 同上, 秒可带小数
//          纯数字 (按小数小时直接解析)
//
//  失败返回 NaN (std::nan)
double timeStr2hour(const std::string &s);

} // namespace year_utils
