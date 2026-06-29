package com.sxwnl.calendar.ui.theme

import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.sp

/**
 * 与鸿蒙 common/Constants.ets AppDimens 对齐
 */
object Dimens {
    // ── 圆角 ─────────────────────────────────────────────
    val radiusSm = 8.dp
    val radiusMd = 12.dp
    val radiusLg = 16.dp
    val radiusXl = 20.dp

    // ── 内边距 ───────────────────────────────────────────
    val paddingXs = 4.dp
    val paddingSm = 8.dp
    val paddingMd = 12.dp
    val paddingLg = 16.dp
    val paddingXl = 20.dp
    val paddingXxl = 24.dp

    // ── 文字尺寸 ────────────────────────────────────────
    val fontTitle = 20.sp
    val fontSubtitle = 16.sp
    val fontBody = 14.sp
    val fontCaption = 12.sp
    val fontSmall = 10.sp

    // ── 控件尺寸 ────────────────────────────────────────
    val cardElevation = 2.dp
    val calendarCellSize = 64.dp
    val tabBarHeight = 56.dp
}

/** 鸿蒙 WEEK_NAMES */
val WeekNames: List<String> = listOf("日", "一", "二", "三", "四", "五", "六")

/** 鸿蒙 TIAN_GAN */
val TianGan: List<String> = listOf(
    "甲", "乙", "丙", "丁", "戊", "己", "庚", "辛", "壬", "癸"
)

/** 鸿蒙 DI_ZHI */
val DiZhi: List<String> = listOf(
    "子", "丑", "寅", "卯", "辰", "巳",
    "午", "未", "申", "酉", "戌", "亥"
)

/** 鸿蒙 SHENG_XIAO */
val ShengXiaoList: List<String> = listOf(
    "鼠", "牛", "虎", "兔", "龙", "蛇",
    "马", "羊", "猴", "鸡", "狗", "猪"
)
