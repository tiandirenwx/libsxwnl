import Foundation

// ════════════════════════════════════════════════════════════════
//  YearUtil — 与鸿蒙 common/YearUtil.ets 对齐的天文纪年工具
//
//  天文学纪年 (astronomical year, ISO 8601):
//    公元 1 年 = 1, 无 0 年; 公元前 1 年 = 0, 公元前 211 年 = -210
//
//  支持的输入格式 (用于年份输入框):
//    "2026"   → 2026
//    "-211"   → -211      (公元前 212 年)
//    "B212"   → -211      (公元前 212 年, 'B' = Before)
//    "b212"   → -211
//    "*212"   → -211      (与 sxwnl JS 风格兼容)
//    "公元前 221 年" → -220
//    "前 221"          → -220
//
//  显示文本:
//    astroYearToStr(2026)  → "2026"
//    astroYearToStr(-211)  → "B212"   (简写, 与鸿蒙月历输入一致)
//    astroYearToFull(-211) → "公元前212年"
// ════════════════════════════════════════════════════════════════

enum YearUtil {

    static let minAstroYear: Int = -4712
    static let maxAstroYear: Int = 9999

    static func isAstroYearValid(_ y: Int) -> Bool {
        y >= minAstroYear && y <= maxAstroYear
    }

    /// 字符串解析为天文学纪年
    static func yearStrToAstro(_ s: String) -> Int {
        let t = s.trimmingCharacters(in: .whitespaces)
        if t.isEmpty { return Int.min }
        // 优先调用 C API (与原 libsxwnl 共用解析逻辑)
        let result = t.withCString { sxwnl_year_str_to_astro($0) }
        return Int(result)
    }

    /// 天文学纪年 → 字符串
    ///   full=false: "2026" / "B212"  (输入框/简写)
    ///   full=true:  "公元2026年" / "公元前212年"
    static func astroYearToStr(_ y: Int, full: Bool = false) -> String {
        var buf = [CChar](repeating: 0, count: 64)
        let ok = buf.withUnsafeMutableBufferPointer {
            sxwnl_astro_year_to_str(Int32(y), full, $0.baseAddress, 64)
        }
        if ok == 0 {
            return String(cString: buf)
        }
        // C API 失败时退回到本地实现
        if y > 0 { return full ? "公元\(y)年" : "\(y)" }
        let before = 1 - y
        return full ? "公元前\(before)年" : "B\(before)"
    }

    /// "公元前212年" 完整版
    static func astroYearToFull(_ y: Int) -> String {
        astroYearToStr(y, full: true)
    }
}
