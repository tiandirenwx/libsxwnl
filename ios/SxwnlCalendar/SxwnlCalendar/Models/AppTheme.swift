import Foundation
import SwiftUI

// ════════════════════════════════════════════════════════════════
//  AppTheme — 与鸿蒙 common/Constants.ets AppColors/AppDimens 对齐
// ════════════════════════════════════════════════════════════════

enum AppColors {
    // ── 主色板 ─────────────────────────────────────────────
    static let primary       = Color(hex: 0x1A237E)
    static let primaryLight  = Color(hex: 0x534BAE)
    static let primaryDark   = Color(hex: 0x000051)
    static let accent        = Color(hex: 0xFFB300)
    static let accentLight   = Color(hex: 0xFFE54C)

    // ── 背景/表面 ──────────────────────────────────────────
    static let background    = Color(hex: 0xF5F5F5)
    static let surface       = Color(hex: 0xFFFFFF)
    static let onPrimary     = Color(hex: 0xFFFFFF)
    static let onSurface     = Color(hex: 0x212121)

    // ── 文字 ──────────────────────────────────────────────
    static let textSecondary = Color(hex: 0x757575)
    static let textTertiary  = Color(hex: 0x9E9E9E)

    // ── 装饰 ──────────────────────────────────────────────
    static let divider       = Color(hex: 0xE0E0E0)
    static let todayHighlight = Color(hex: 0xFFB300)
    static let jieQi         = Color(hex: 0xE65100)
    static let weekend       = Color(hex: 0xC62828)
    static let lunarText     = Color(hex: 0x9E9E9E)
    static let holiday       = Color(hex: 0xD32F2F)

    static let gradientStart = Color(hex: 0x1A237E)
    static let gradientEnd   = Color(hex: 0x283593)
    static let goldStart     = Color(hex: 0xFFB300)
    static let goldEnd       = Color(hex: 0xFF8F00)

    // ── 日月食专属色板 ───────────────────────────────────────
    // 与 Android Color.kt 完全对应
    static let skyDeepNight = Color(hex: 0x070B22)
    static let skyMidNight  = Color(hex: 0x121A36)
    static let skyDawn      = Color(hex: 0x2B2347)
    static let sunCore      = Color(hex: 0xFFE082)
    static let sunGlow      = Color(hex: 0xFFB300)
    static let sunHalo      = Color(hex: 0x4C2B00)
    static let moonDark     = Color(hex: 0x12131A)
    static let moonRim      = Color(hex: 0x2A2F45)
    static let moonBright   = Color(hex: 0xE6E6FA)
    static let moonBlood    = Color(hex: 0xB23A2E)
    static let earthUmbra   = Color(hex: 0x000814)
    static let earthPenumbra = Color(hex: 0x1E2840)
    static let mapLand      = Color(hex: 0x223453)
    static let mapOcean     = Color(hex: 0x0B1124)
    static let mapGrid      = Color(hex: 0x172040)
    static let mapBorder    = Color(hex: 0x4C5A86)
    static let pathCenter   = Color(hex: 0xFFC107)
    static let pathUmbra    = Color(hex: 0xEF5350)
    static let pathPenumbra = Color(hex: 0x90CAF9)
    static let pathMarker   = Color(hex: 0xFFFFFF)
}

enum AppDimens {
    static let radiusSm: CGFloat = 8
    static let radiusMd: CGFloat = 12
    static let radiusLg: CGFloat = 16
    static let radiusXl: CGFloat = 20

    static let paddingXs: CGFloat = 4
    static let paddingSm: CGFloat = 8
    static let paddingMd: CGFloat = 12
    static let paddingLg: CGFloat = 16
    static let paddingXl: CGFloat = 20
    static let paddingXxl: CGFloat = 24

    static let fontTitle: CGFloat = 20
    static let fontSubtitle: CGFloat = 16
    static let fontBody: CGFloat = 14
    static let fontCaption: CGFloat = 12
    static let fontSmall: CGFloat = 10

    static let calendarCellSize: CGFloat = 64
    static let tabBarHeight: CGFloat = 56
}

enum AppText {
    static let weekNames: [String] = ["日", "一", "二", "三", "四", "五", "六"]
    static let tianGan: [String] = ["甲", "乙", "丙", "丁", "戊",
                                    "己", "庚", "辛", "壬", "癸"]
    static let diZhi: [String] = ["子", "丑", "寅", "卯", "辰", "巳",
                                  "午", "未", "申", "酉", "戌", "亥"]
    static let shengXiao: [String] = ["鼠", "牛", "虎", "兔", "龙", "蛇",
                                      "马", "羊", "猴", "鸡", "狗", "猪"]
}

// MARK: - Color 扩展: 0xRRGGBB 整数构造

extension Color {
    init(hex: UInt32, alpha: Double = 1.0) {
        let r = Double((hex >> 16) & 0xFF) / 255.0
        let g = Double((hex >>  8) & 0xFF) / 255.0
        let b = Double( hex        & 0xFF) / 255.0
        self.init(.sRGB, red: r, green: g, blue: b, opacity: alpha)
    }
}
