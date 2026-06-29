package com.sxwnl.calendar.util

import com.sxwnl.calendar.bridge.SxwnlNative

/**
 * 与鸿蒙端 common/YearUtil.ets 对齐的天文纪年工具。
 *
 *   天文学纪年 (ISO 8601):
 *     公元 1 年 = 1, 无 0 年; 公元前 1 年 = 0, 公元前 211 年 = -210
 *
 *   支持的输入格式 (用于年份输入框):
 *     "2026"   → 2026
 *     "-211"   → -211      (公元前 212 年)
 *     "B212"   → -211      (公元前 212 年, 'B' = Before)
 *     "b212"   → -211
 *     "*212"   → -211
 *     "公元前 221 年" → -220
 *
 *   显示文本:
 *     astroYearToStr(2026)       → "2026"
 *     astroYearToStr(-211)       → "B212"      (简写, 与鸿蒙月历输入一致)
 *     astroYearToStr(-211, true) → "公元前212年"
 */
object YearUtil {

    const val MIN_ASTRO_YEAR = -4712
    const val MAX_ASTRO_YEAR = 9999

    fun isAstroYearValid(y: Int): Boolean = y in MIN_ASTRO_YEAR..MAX_ASTRO_YEAR

    fun yearStrToAstro(s: String): Int {
        val t = s.trim()
        if (t.isEmpty()) return Int.MIN_VALUE
        return SxwnlNative.yearStrToAstro(t)
    }

    fun astroYearToStr(year: Int, full: Boolean = false): String {
        val r = SxwnlNative.astroYearToStr(year, full)
        if (r.isNotEmpty()) return r
        // 回退到本地实现
        return if (year > 0) {
            if (full) "公元${year}年" else "$year"
        } else {
            val before = 1 - year
            if (full) "公元前${before}年" else "B$before"
        }
    }

    fun astroYearToFull(year: Int): String = astroYearToStr(year, true)
}
