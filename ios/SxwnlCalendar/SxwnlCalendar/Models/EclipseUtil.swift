import Foundation
import SwiftUI

// ════════════════════════════════════════════════════════════════
//  EclipseUtil — 日月食类型分类 + JD 格式化
//
//  与 Android 端 EclipseUtil.kt / 鸿蒙端 EclipseUtil.ets 行为对齐
// ════════════════════════════════════════════════════════════════

enum EclipseUtil {

    /// 日食类型徽章颜色
    static func solarBadgeColor(_ type: String) -> Color {
        switch type.first {
        case "T": return Color(hex: 0xD32F2F)  // 全食 — 红
        case "A": return Color(hex: 0xF57C00)  // 环食 — 橙
        case "H": return Color(hex: 0x7B1FA2)  // 全环混合 — 紫
        case "P": return Color(hex: 0x616161)  // 偏食 — 灰
        default:  return Color(hex: 0x9E9E9E)
        }
    }

    static func lunarBadgeColor(_ type: String) -> Color {
        switch type {
        case "T": return Color(hex: 0xD32F2F)
        case "P": return Color(hex: 0x616161)
        case "B": return Color(hex: 0x1976D2)
        default:  return Color(hex: 0x9E9E9E)
        }
    }

    /// 力学时 JD → "YYYY-MM-DD HH:MM:SS"
    static func jdToDateTime(_ jd: Double) -> String {
        guard jd > 0 else { return "--" }
        let jdAdj = jd + 0.5
        let Z = Int(floor(jdAdj))
        let F = jdAdj - Double(Z)
        var A = Z
        if Z >= 2_299_161 {
            let alpha = Int((Double(Z) - 1_867_216.25) / 36_524.25)
            A = Z + 1 + alpha - alpha / 4
        }
        let B = A + 1524
        let C = Int((Double(B) - 122.1) / 365.25)
        let D = Int(365.25 * Double(C))
        let E = Int((Double(B - D)) / 30.6001)
        let day = Double(B - D - Int(30.6001 * Double(E))) + F
        let mon = E < 14 ? E - 1 : E - 13
        let yr  = mon > 2 ? C - 4716 : C - 4715
        let dInt = Int(floor(day))
        var dFrac = day - Double(dInt)
        var hh = Int(dFrac * 24); dFrac = dFrac * 24 - Double(hh)
        var mm = Int(dFrac * 60); dFrac = dFrac * 60 - Double(mm)
        var ss = Int(dFrac * 60 + 0.5)
        if ss >= 60 { ss = 0; mm += 1 }
        if mm >= 60 { mm = 0; hh += 1 }
        return String(format: "%04d-%02d-%02d %02d:%02d:%02d",
                      yr, mon, dInt, hh, mm, ss)
    }

    static func jdToTime(_ jd: Double) -> String {
        guard jd > 0 else { return "--" }
        let s = jdToDateTime(jd)
        return s.count >= 19 ? String(s.suffix(8)) : "--"
    }

    /// 食分 → 0..1 进度
    static func magnitudeProgress(_ mag: Double) -> Double {
        max(0, min(1.5, mag)) / 1.5
    }

    /// 秒 → "Hh MMm SSs" / "MMm SSs"
    static func formatDuration(_ sec: Double) -> String {
        if sec <= 0 { return "—" }
        let s = Int(sec)
        let h = s / 3600
        let m = (s % 3600) / 60
        let ss = s % 60
        return h > 0
            ? String(format: "%dh%02dm%02ds", h, m, ss)
            : String(format: "%dm%02ds", m, ss)
    }

    // ─── ΔT (TD − UT1) Espenak/Meeus 5000-year polynomial ──────
    static func deltaT(year: Int) -> Double {
        let y = Double(year)
        if y < -500 || y > 2150 {
            let u = (y - 1820) / 100
            return -20 + 32 * u * u
        }
        switch year {
        case ..<(-500): let u = (y - 1820)/100; return -20 + 32*u*u
        case (-500)..<500:
            let u = y/100
            return 10583.6 + u*(-1014.41 + u*(33.78311 + u*(-5.952053 + u*(-0.1798452 + u*(0.022174192 + u*0.0090316521)))))
        case 500..<1600:
            let u = (y-1000)/100
            return 1574.2 + u*(-556.01 + u*(71.23472 + u*(0.319781 + u*(-0.8503463 + u*(-0.005050998 + u*0.0083572073)))))
        case 1600..<1700:
            let t = y-1600
            return 120 + t*(-0.9808 + t*(-0.01532 + t/7129))
        case 1700..<1800:
            let t = y-1700
            return 8.83 + t*(0.1603 + t*(-0.0059285 + t*(0.00013336 - t/1174000)))
        case 1800..<1860:
            let t = y-1800
            return 13.72 + t*(-0.332447 + t*(0.0068612 + t*(0.0041116 + t*(-0.00037436 + t*(0.0000121272 + t*(-0.0000001699 + t*0.000000000875))))))
        case 1860..<1900:
            let t = y-1860
            return 7.62 + t*(0.5737 + t*(-0.251754 + t*(0.01680668 + t*(-0.0004473624 + t/233174))))
        case 1900..<1920:
            let t = y-1900
            return -2.79 + t*(1.494119 + t*(-0.0598939 + t*(0.0061966 - t*0.000197)))
        case 1920..<1941:
            let t = y-1920
            return 21.20 + t*(0.84493 + t*(-0.076100 + t*0.0020936))
        case 1941..<1961:
            let t = y-1950
            return 29.07 + t*(0.407 + t*(-1/233 + t/2547))
        case 1961..<1986:
            let t = y-1975
            return 45.45 + t*(1.067 + t*(-1/260 - t/718))
        case 1986..<2005:
            let t = y-2000
            return 63.86 + t*(0.3345 + t*(-0.060374 + t*(0.0017275 + t*(0.000651814 + t*0.00002373599))))
        case 2005..<2050:
            let t = y-2000
            return 62.92 + t*(0.32217 + t*0.005589)
        case 2050..<2150:
            return -20 + 32*pow((y-1820)/100, 2) - 0.5628*(2150-y)
        default:
            let u = (y-1820)/100; return -20 + 32*u*u
        }
    }

    /// 力学时 JD → 本地时 "HH:MM:SS" (含 ΔT, 时区)
    static func jdTdToLocal(_ jd: Double, tzHours: Double) -> String {
        guard jd > 0 else { return "--" }
        // 估算年份用于 deltaT
        let s = jdToDateTime(jd)
        let yr = Int(s.prefix(4)) ?? 2000
        let dt = deltaT(year: yr) / 86400.0
        let ut = jd - dt + tzHours / 24.0
        return jdToTime(ut)
    }

    static func tzLabel(_ tzHours: Double) -> String {
        let sign = tzHours >= 0 ? "+" : "-"
        let abs = Swift.abs(tzHours)
        let h = Int(abs)
        let m = Int((abs - Double(h)) * 60)
        return String(format: "UTC%@%02d:%02d", sign, h, m)
    }

    // ─── 等距圆柱投影 (Plate Carrée) ─────────────────────────
    /// 经度(弧度) → 屏幕 X 比例 0..1
    static func projectLonX(_ lonRad: Double) -> Double {
        (lonRad + .pi) / (2 * .pi)
    }
    /// 纬度(弧度) → 屏幕 Y 比例 0..1 (北上)
    static func projectLatY(_ latRad: Double) -> Double {
        0.5 - latRad / .pi
    }

    /// 帧线性插值
    static func lerp(_ a: Double, _ b: Double, _ t: Double) -> Double { a + (b - a) * t }

    // ─── 类型 emoji (与 Android / 鸿蒙端一致) ─────────────────
    static func solarEmoji(_ type: String) -> String {
        switch type.first {
        case "T": return "🌑"
        case "A": return "🌒"
        case "H": return "🌓"
        default:  return "🌗"
        }
    }

    static func lunarEmoji(_ type: String) -> String {
        switch type {
        case "T": return "🌕"
        case "P": return "🌖"
        case "B": return "🌗"
        default:  return "🌘"
        }
    }
}
