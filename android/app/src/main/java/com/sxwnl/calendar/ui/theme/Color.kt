package com.sxwnl.calendar.ui.theme

import androidx.compose.ui.graphics.Color

// ════════════════════════════════════════════════════════════════
//  与鸿蒙端 common/Constants.ets AppColors 对齐
// ════════════════════════════════════════════════════════════════

// ── 主色板 ─────────────────────────────────────────────────
val Primary = Color(0xFF1A237E)
val PrimaryLight = Color(0xFF534BAE)
val PrimaryDark = Color(0xFF000051)
val Accent = Color(0xFFFFB300)
val AccentLight = Color(0xFFFFE54C)

// ── 背景/表面 ──────────────────────────────────────────────
val Background = Color(0xFFF5F5F5)
val Surface = Color(0xFFFFFFFF)
val OnPrimary = Color(0xFFFFFFFF)
val OnSurface = Color(0xFF212121)

// ── 文字 ──────────────────────────────────────────────────
val TextPrimary = Color(0xFF1C1B1F)
val TextSecondary = Color(0xFF757575)
val TextTertiary = Color(0xFF9E9E9E)

// ── 分隔/装饰 ─────────────────────────────────────────────
val DividerColor = Color(0xFFE0E0E0)
val TodayHighlight = Color(0xFFFFB300)
val JieQiColor = Color(0xFFE65100)
val WeekendColor = Color(0xFFC62828)
val LunarText = Color(0xFF9E9E9E)
val HolidayRed = Color(0xFFD32F2F)

val GradientStart = Color(0xFF1A237E)
val GradientEnd = Color(0xFF283593)
val GoldGradientStart = Color(0xFFFFB300)
val GoldGradientEnd = Color(0xFFFF8F00)
val CardShadow = Color(0x1A000000)

// ── 四柱卡片 ──────────────────────────────────────────────
val PillarBgYear = Color(0xFFE8EAF6)
val PillarBgMonth = Color(0xFFFFF8E1)
val PillarBgDay = Color(0xFFE0F2F1)
val PillarBgHour = Color(0xFFFCE4EC)

// ── 旧别名 (BaziResultView 等保留) ────────────────────────
val DeepIndigo = Primary
val Indigo700 = Color(0xFF283593)
val Indigo500 = Color(0xFF3F51B5)
val Indigo300 = Color(0xFF7986CB)
val Indigo100 = Color(0xFFC5CAE9)
val Indigo50 = Color(0xFFE8EAF6)
val WarmGold = Accent
val GoldLight = Color(0xFFFFD54F)
val GoldDark = Color(0xFFF9A825)
val SurfaceWhite = Surface
val SurfaceVariant = Color(0xFFF5F5F5)
val BackgroundWhite = Surface
val ErrorRed = Color(0xFFB3261E)
val SuccessGreen = Color(0xFF2E7D32)
val JieQiGreen = JieQiColor

// ── 日月食专属色 (与 iOS AppTheme.swift / 鸿蒙 Constants.ets 对齐) ──
// 天空夜色背景, 让圆盘动画有质感
val SkyDeepNight  = Color(0xFF070B22)
val SkyMidNight   = Color(0xFF121A36)
val SkyDawn       = Color(0xFF2B2347)
// 太阳: 暖金色调
val SunCore       = Color(0xFFFFE082)
val SunGlow       = Color(0xFFFFB300)
val SunHalo       = Color(0x664C2B00)
// 月亮 (日食中遮挡太阳的月亮)
val MoonDark      = Color(0xFF12131A)
val MoonRim       = Color(0xFF2A2F45)
// 月食中的月亮 (本影时呈血月)
val MoonBright    = Color(0xFFE6E6FA)
val MoonBlood     = Color(0xFFB23A2E)
// 地球本影/半影
val EarthUmbra    = Color(0xFF000814)
val EarthPenumbra = Color(0x661E2840)
// 地图色
val MapLand       = Color(0xFF223453)
val MapOcean      = Color(0xFF0B1124)
val MapGrid       = Color(0xFF172040)
val MapBorder     = Color(0xFF4C5A86)
// 食带/路径色
val PathCenter    = Color(0xFFFFC107)   // 中心线 (金黄)
val PathUmbra     = Color(0xFFEF5350)   // 本影界 (亮红)
val PathPenumbra  = Color(0xFF90CAF9)   // 半影界 (淡蓝)
val PathMarker    = Color(0xFFFFFFFF)   // 移动影点 (白)
