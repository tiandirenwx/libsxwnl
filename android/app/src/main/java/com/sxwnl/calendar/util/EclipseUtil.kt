package com.sxwnl.calendar.util

import androidx.compose.ui.graphics.Color
import kotlin.math.abs
import kotlin.math.cos
import kotlin.math.floor
import kotlin.math.sin
import kotlin.math.sqrt

// ════════════════════════════════════════════════════════════════
//  EclipseUtil — 日月食通用工具
//
//  ① 类型徽章配色 (日食: T/A/H/P; 月食: T/P/B)
//  ② JD ↔ 公历, ΔT (Espenak/Meeus 五千年模型), TD↔UT↔本地时
//  ③ 等距圆柱地图投影 (lon/lat 弧度 → [-1,1])
//  ④ 食分进度、时长格式化、弧度度分秒
//  ⑤ 帧线性插值
// ════════════════════════════════════════════════════════════════
object EclipseUtil {

    // ── ① 徽章配色 ─────────────────────────────────────────────

    /** 日食类型徽章颜色; T=全, A=环, H=全环, P=偏 */
    fun solarBadgeColor(type: String): Color = when {
        type.startsWith("T") -> Color(0xFFD32F2F)
        type.startsWith("A") -> Color(0xFFF57C00)
        type.startsWith("H") -> Color(0xFF7B1FA2)
        type == "P"           -> Color(0xFF616161)
        else                  -> Color(0xFF9E9E9E)
    }

    /** 月食类型徽章颜色; T=全, P=偏, B=半影 */
    fun lunarBadgeColor(type: String): Color = when (type) {
        "T"  -> Color(0xFFD32F2F)
        "P"  -> Color(0xFF616161)
        "B"  -> Color(0xFF1976D2)
        else -> Color(0xFF9E9E9E)
    }

    /** 日食类型 emoji (与 iOS / 鸿蒙端一致) */
    fun solarEmoji(type: String): String = when {
        type.startsWith("T") -> "🌑"
        type.startsWith("A") -> "🌒"
        type.startsWith("H") -> "🌓"
        else                  -> "🌗"
    }

    /** 月食类型 emoji (与 iOS / 鸿蒙端一致) */
    fun lunarEmoji(type: String): String = when (type) {
        "T"  -> "🌕"
        "P"  -> "🌖"
        "B"  -> "🌗"
        else -> "🌘"
    }

    // ── ② JD / ΔT / 时间转换 ───────────────────────────────────

    /**
     * ΔT = TD − UT, 单位: 秒. Espenak & Meeus 五千年多项式 (适用 -1999..+3000).
     * 参考: NASA Eclipse Web Site, "Polynomial Expressions for Delta T".
     */
    fun deltaT(year: Int, month: Int = 6): Double {
        val y = year + (month - 0.5) / 12.0
        return when {
            y < -500 -> {
                val u = (y - 1820) / 100.0
                -20.0 + 32.0 * u * u
            }
            y < 500 -> {
                val u = y / 100.0
                10583.6 + u * (-1014.41 + u * (33.78311 + u * (-5.952053 +
                    u * (-0.1798452 + u * (0.022174192 + u * 0.0090316521)))))
            }
            y < 1600 -> {
                val u = (y - 1000) / 100.0
                1574.2 + u * (-556.01 + u * (71.23472 + u * (0.319781 +
                    u * (-0.8503463 + u * (-0.005050998 + u * 0.0083572073)))))
            }
            y < 1700 -> {
                val t = y - 1600
                120.0 - 0.9808 * t - 0.01532 * t * t + t * t * t / 7129.0
            }
            y < 1800 -> {
                val t = y - 1700
                8.83 + t * (0.1603 + t * (-0.0059285 + t * (0.00013336 - t / 1174000.0)))
            }
            y < 1860 -> {
                val t = y - 1800
                13.72 + t * (-0.332447 + t * (0.0068612 + t * (0.0041116 + t * (-0.00037436 +
                    t * (0.0000121272 + t * (-0.0000001699 + t * 0.000000000875))))))
            }
            y < 1900 -> {
                val t = y - 1860
                7.62 + t * (0.5737 + t * (-0.251754 + t * (0.01680668 +
                    t * (-0.0004473624 + t / 233174.0))))
            }
            y < 1920 -> {
                val t = y - 1900
                -2.79 + t * (1.494119 + t * (-0.0598939 + t * (0.0061966 - t * 0.000197)))
            }
            y < 1941 -> {
                val t = y - 1920
                21.20 + t * (0.84493 + t * (-0.076100 + t * 0.0020936))
            }
            y < 1961 -> {
                val t = y - 1950
                29.07 + 0.407 * t - t * t / 233.0 + t * t * t / 2547.0
            }
            y < 1986 -> {
                val t = y - 1975
                45.45 + 1.067 * t - t * t / 260.0 - t * t * t / 718.0
            }
            y < 2005 -> {
                val t = y - 2000
                63.86 + t * (0.3345 + t * (-0.060374 + t * (0.0017275 +
                    t * (0.000651814 + t * 0.00002373599))))
            }
            y < 2050 -> {
                val t = y - 2000
                62.92 + t * (0.32217 + t * 0.005589)
            }
            y < 2150 -> {
                val u = (y - 1820) / 100.0
                -20 + 32 * u * u - 0.5628 * (2150 - y)
            }
            else -> {
                val u = (y - 1820) / 100.0
                -20 + 32 * u * u
            }
        }
    }

    /** Meeus 算法: JD → 公历日期 (返回 Y/M/D 整数 + 小时小数) */
    private data class Ymd(val y: Int, val m: Int, val d: Int, val hourFrac: Double)
    private fun jdToYmd(jd: Double): Ymd {
        val jdAdj = jd + 0.5
        val Z = floor(jdAdj).toLong()
        val F = jdAdj - Z
        var A = Z
        if (Z >= 2299161L) {
            val a = ((Z - 1867216L) - 0.25) / 36524.25
            A = Z + 1 + a.toLong() - (a / 4).toLong()
        }
        val B = A + 1524
        val C = ((B - 122.1) / 365.25).toLong()
        val D = (365.25 * C).toLong()
        val E = ((B - D) / 30.6001).toLong()
        val day = (B - D - (30.6001 * E).toLong()) + F
        val mon = if (E < 14) (E - 1).toInt() else (E - 13).toInt()
        val yr  = if (mon > 2) (C - 4716).toInt() else (C - 4715).toInt()
        val dInt = floor(day).toInt()
        val dFrac = day - dInt
        return Ymd(yr, mon, dInt, dFrac * 24.0)
    }

    private fun hmsFromFrac(hourFrac: Double): Triple<Int, Int, Int> {
        var hh = hourFrac.toInt()
        var rest = hourFrac - hh
        var mm = (rest * 60).toInt()
        rest = rest * 60 - mm
        var ss = (rest * 60 + 0.5).toInt()
        if (ss >= 60) { ss = 0; mm += 1 }
        if (mm >= 60) { mm = 0; hh += 1 }
        return Triple(hh, mm, ss)
    }

    /** 力学时 (TD) JD → "YYYY-MM-DD HH:MM:SS TD" */
    fun jdToDateTime(jd: Double): String {
        if (jd <= 0.0) return "--"
        val (y, m, d, hf) = jdToYmd(jd)
        val (hh, mm, ss) = hmsFromFrac(hf)
        return "%04d-%02d-%02d %02d:%02d:%02d".format(y, m, d, hh, mm, ss)
    }

    /** 仅 HH:MM:SS */
    fun jdToTime(jd: Double): String {
        if (jd <= 0.0) return "--"
        return jdToDateTime(jd).substring(11)
    }

    /** TD → 本地民用时间字符串.
     *  tzHours = 时区偏移 (东 +), 默认 +8 (北京)
     *  withDate=true 包含日期
     */
    fun jdTdToLocal(jd: Double, tzHours: Double = 8.0, withDate: Boolean = false): String {
        if (jd <= 0.0) return "--"
        val (y0, m0, _, _) = jdToYmd(jd)
        val dT = deltaT(y0, m0) / 86400.0
        val localJd = jd - dT + tzHours / 24.0
        val (y, m, d, hf) = jdToYmd(localJd)
        val (hh, mm, ss) = hmsFromFrac(hf)
        return if (withDate)
            "%04d-%02d-%02d %02d:%02d:%02d".format(y, m, d, hh, mm, ss)
        else "%02d:%02d:%02d".format(hh, mm, ss)
    }

    /** "+0800" 这种短标签 */
    fun tzLabel(tzHours: Double): String {
        val sign = if (tzHours >= 0) "+" else "-"
        val abs = abs(tzHours)
        val h = abs.toInt()
        val m = ((abs - h) * 60).toInt()
        return "UTC%s%02d:%02d".format(sign, h, m)
    }

    // ── ③ 等距圆柱投影 (世界地图绘制用) ────────────────────────

    /**
     * 经度 (弧度) → 屏幕 X 比例 [0, 1]
     * 输入 lon 是 [-π, π]; 输出 0 = 左边 (-180°), 1 = 右边 (180°)
     */
    fun projectLonX(lonRad: Double): Float {
        val norm = ((lonRad + Math.PI).rem(2 * Math.PI) + 2 * Math.PI).rem(2 * Math.PI)
        return (norm / (2 * Math.PI)).toFloat()
    }

    /** 纬度 (弧度) → 屏幕 Y 比例 [0, 1]; 0=北极, 1=南极 (屏幕 Y 向下) */
    fun projectLatY(latRad: Double): Float {
        val clamped = latRad.coerceIn(-Math.PI / 2, Math.PI / 2)
        return (0.5 - clamped / Math.PI).toFloat()
    }

    // ── ④ 格式化辅助 ───────────────────────────────────────────

    /** 食分 → 0..1 进度 (>=1 视为满) */
    fun magnitudeProgress(mag: Double): Float =
        (mag.coerceAtLeast(0.0).coerceAtMost(1.5) / 1.5).toFloat()

    /** 秒 → "Hh MMm SSs" / "MMm SSs" */
    fun formatDuration(sec: Double): String {
        if (sec <= 0) return "—"
        val s = sec.toInt()
        val h = s / 3600
        val m = (s % 3600) / 60
        val ss = s % 60
        return if (h > 0) "%dh%02dm%02ds".format(h, m, ss)
        else "%dm%02ds".format(m, ss)
    }

    /** 经/纬度漂亮格式: "116.40°E", "39.90°N" */
    fun lonLabel(lonDeg: Double): String {
        val v = ((lonDeg + 540).rem(360)) - 180
        val hemi = if (v >= 0) "E" else "W"
        return "%.2f°%s".format(abs(v), hemi)
    }
    fun latLabel(latDeg: Double): String {
        val hemi = if (latDeg >= 0) "N" else "S"
        return "%.2f°%s".format(abs(latDeg), hemi)
    }

    // ── ⑤ 帧线性插值 (动画 scrubber 用) ───────────────────────

    /**
     * 在已排序的 frames 中按 jd 找到合适帧 (返回索引 + 插值 t∈[0,1]).
     */
    fun frameLerp(frameJds: DoubleArray, targetJd: Double): Pair<Int, Float> {
        if (frameJds.isEmpty()) return 0 to 0f
        if (targetJd <= frameJds.first()) return 0 to 0f
        if (targetJd >= frameJds.last()) return (frameJds.size - 1) to 0f
        // 二分查找
        var lo = 0
        var hi = frameJds.size - 1
        while (lo + 1 < hi) {
            val mid = (lo + hi) ushr 1
            if (frameJds[mid] <= targetJd) lo = mid else hi = mid
        }
        val span = frameJds[hi] - frameJds[lo]
        val t = if (span > 0) ((targetJd - frameJds[lo]) / span).toFloat() else 0f
        return lo to t
    }

    fun lerp(a: Double, b: Double, t: Float): Double = a + (b - a) * t
    fun lerp(a: Float, b: Float, t: Float): Float = a + (b - a) * t
}
