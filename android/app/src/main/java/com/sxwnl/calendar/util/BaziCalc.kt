package com.sxwnl.calendar.util

import com.sxwnl.calendar.data.BaziParams
import com.sxwnl.calendar.data.LunarMonth

/**
 * 与鸿蒙端 model/BaziInput.ets 对齐的八字纯计算工具(无 UI 依赖)。
 */
object BaziCalc {

    val GAN: Array<String> = arrayOf(
        "甲", "乙", "丙", "丁", "戊", "己", "庚", "辛", "壬", "癸"
    )
    val ZHI: Array<String> = arrayOf(
        "子", "丑", "寅", "卯", "辰", "巳",
        "午", "未", "申", "酉", "戌", "亥"
    )

    /** 输入模式 */
    const val MODE_SOLAR = 0
    const val MODE_LUNAR = 1
    const val MODE_REVERSE = 2

    /** 输入模式 — 枚举形式 (与鸿蒙端 BaziInput.ets 的整型常量等价) */
    enum class BirthInputMode(val rawValue: Int) {
        SOLAR(MODE_SOLAR), LUNAR(MODE_LUNAR), REVERSE(MODE_REVERSE);

        companion object {
            fun of(v: Int): BirthInputMode = when (v) {
                MODE_LUNAR -> LUNAR
                MODE_REVERSE -> REVERSE
                else -> SOLAR
            }
        }
    }

    /** 出生时间选择结果 */
    data class DateSelection(
        val inputMode: BirthInputMode = BirthInputMode.SOLAR,
        val year: Int = 1990,
        val month: Int = 1,
        val day: Int = 1,
        val hour: Int = 12,
        val minute: Int = 0,
        val isLeap: Boolean = false,
        val isSpec: Boolean = false
    )

    /** 历法 — 仅 3 项 (与鸿蒙 LIFA_OPTIONS 一致) */
    const val LIFA_DING_QI = 11
    const val LIFA_PING_DONG_ZHI = 12
    const val LIFA_PING_XIA_ZHI = 13

    data class LiFaOption(val label: String, val value: Int)

    val LIFA_OPTIONS: List<LiFaOption> = listOf(
        LiFaOption("定气", LIFA_DING_QI),
        LiFaOption("平气定夏至", LIFA_PING_XIA_ZHI),
        LiFaOption("平气定冬至", LIFA_PING_DONG_ZHI)
    )

    fun jiaZiName(idx: Int): String {
        val i = ((idx % 60) + 60) % 60
        return GAN[i % 10] + ZHI[i % 12]
    }

    fun gz60(gan: Int, zhi: Int): Int {
        for (n in 0 until 60) {
            if (n % 10 == gan && n % 12 == zhi) return n
        }
        return 0
    }

    /** 五虎遁: 由年干得 12 个合法月柱(寅月起)的六十甲子序号 */
    fun validMonths(yearIdx: Int): IntArray {
        val yg = ((yearIdx % 60) + 60) % 60 % 10
        val base = ((yg % 5) * 2 + 2) % 10
        return IntArray(12) { k -> gz60((base + k) % 10, (2 + k) % 12) }
    }

    /** 五鼠遁: 由日干得 12 个合法时柱(子时起)的六十甲子序号 */
    fun validHours(dayIdx: Int): IntArray {
        val dg = ((dayIdx % 60) + 60) % 60 % 10
        val base = ((dg % 5) * 2) % 10
        return IntArray(12) { h -> gz60((base + h) % 10, h) }
    }

    fun remapByZhi(validList: IntArray, oldIdx: Int): Int {
        val z = ((oldIdx % 60) + 60) % 60 % 12
        for (v in validList) {
            if (v % 12 == z) return v
        }
        return validList.firstOrNull() ?: 0
    }

    fun lnGan(year: Int): Int = (((year - 4) % 10) + 10) % 10
    fun lnZhi(year: Int): Int = (((year - 4) % 12) + 12) % 12

    fun parseYear(s: String): Int? {
        val v = YearUtil.yearStrToAstro(s)
        return if (v == Int.MIN_VALUE) null else v
    }

    fun formatYear(y: Int): String = if (y <= 0) "公元前${1 - y}年" else "${y}年"

    fun solarDaysInMonth(year: Int, month: Int): Int = when (month) {
        2 -> {
            val leap = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)
            if (leap) 29 else 28
        }
        4, 6, 9, 11 -> 30
        else -> 31
    }

    fun pad2(n: Int): String = if (n < 10) "0$n" else "$n"

    /**
     * 农历日名(初一../三十). 统一走底层 C++(复用 sx_lang_zh.h 的 Rmc[] 表),
     * 不在各平台重复维护字符串表; native 异常时回退 "${d}日".
     */
    fun lunarDayName(d: Int): String {
        val s = com.sxwnl.calendar.bridge.SxwnlNative.getLunarDayName(d)
        return s.ifEmpty { "${d}日" }
    }

    /** 出生时间一条记录文本 */
    fun formatRecord(p: BaziParams, lunarMonths: List<LunarMonth>): String {
        val ymd = "${formatYear(p.year)}${pad2(p.month)}月${pad2(p.day)}日"
        val hm = "${p.hour}点${pad2(p.minute)}分"
        if (p.inputMode == MODE_LUNAR) {
            val mName = lunarMonths.firstOrNull {
                it.month == p.month && it.isLeap == p.isLeap
            }?.name ?: "${p.month}月"
            return "农历 ${formatYear(p.year)}${mName}${lunarDayName(p.day)} $hm"
        }
        return "公历 $ymd $hm"
    }
}
