// ─────────────────────────────────────────────────────────────
//  公元/天文纪年转换 (对应 sxwnl/tools.js: year2Ayear / Ayear2year)
//
//  ▷ 普通纪年 ↔ 天文纪年 约定:
//    天文年 0  = 公元前 1 年       (B1)
//    天文年 1  = 公元 1 年
//    天文年 -211 = 公元前 212 年    (B212)
//
//  ▷ 输入支持的格式:
//    "2026"   → 公元 2026 年       → 天文年 2026
//    "-211"   → 直接当天文年        → -211       (=公元前212)
//    "B212"   → 公元前 212 年       → 天文年 1-212 = -211
//    "b212"   → 同上 (大小写不敏感)
//    "*212"   → 同上
//
//  ▷ astroYearToStr 显示约定:
//    y > 0      → "2026"
//    y <= 0     → "B" + (1-y)   (即  0 → "B1",  -211 → "B212")
// ─────────────────────────────────────────────────────────────
export class YearUtil {
  // 字符串 → 天文年; 解析失败返回 NaN
  static yearStrToAstro(input: string): number {
    if (input === null || input === undefined) return NaN
    // 保留 0-9, B, b, *, -
    const cleaned = ('' + input).replace(/[^0-9Bb\*\-]/g, '')
    if (cleaned.length === 0) return NaN

    const head = cleaned.substr(0, 1)
    if (head === 'B' || head === 'b' || head === '*') {
      const n = parseInt(cleaned.substr(1), 10)
      // 公元前从 B.C.1 开始, 没有公元 0 年
      if (isNaN(n) || n < 1) return NaN
      return 1 - n
    }
    const v = parseInt(cleaned, 10)
    return isNaN(v) ? NaN : v
  }

  // 天文年 → 显示字符串 (公元后直接数字, 公元前用 B 前缀)
  static astroYearToStr(y: number): string {
    if (isNaN(y)) return ''
    if (y <= 0) return 'B' + (1 - y)
    return '' + y
  }

  // 校验天文年是否在本库可计算范围内 (跟 JS 约定一致: -4712 ~ 9999)
  static isAstroYearValid(y: number): boolean {
    return !isNaN(y) && y >= -4712 && y <= 9999
  }
}
