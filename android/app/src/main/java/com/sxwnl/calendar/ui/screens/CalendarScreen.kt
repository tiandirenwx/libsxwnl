package com.sxwnl.calendar.ui.screens

import androidx.compose.foundation.background
import androidx.compose.foundation.border
import androidx.compose.foundation.clickable
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.lazy.grid.GridCells
import androidx.compose.foundation.lazy.grid.LazyVerticalGrid
import androidx.compose.foundation.lazy.grid.items
import androidx.compose.foundation.rememberScrollState
import androidx.compose.foundation.shape.CircleShape
import androidx.compose.foundation.shape.RoundedCornerShape
import androidx.compose.foundation.text.BasicTextField
import androidx.compose.foundation.text.KeyboardOptions
import androidx.compose.foundation.verticalScroll
import androidx.compose.material3.*
import androidx.compose.runtime.*
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.draw.clip
import androidx.compose.ui.graphics.Brush
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.text.TextStyle
import androidx.compose.ui.text.font.FontWeight
import androidx.compose.ui.text.input.ImeAction
import androidx.compose.ui.text.input.KeyboardType
import androidx.compose.ui.text.style.TextAlign
import androidx.compose.ui.unit.Dp
import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.sp
import androidx.compose.ui.window.Dialog
import androidx.compose.ui.window.DialogProperties
import com.sxwnl.calendar.data.CalendarRepository
import com.sxwnl.calendar.data.DayInfo
import com.sxwnl.calendar.data.DayRTS
import com.sxwnl.calendar.ui.theme.*
import com.sxwnl.calendar.util.YearUtil
import kotlinx.coroutines.launch
import java.util.Calendar

// ════════════════════════════════════════════════════════════════
//  CalendarScreen — 与鸿蒙端 CalendarPage.ets 对齐的月历页面
// ════════════════════════════════════════════════════════════════

/** 月相 / 节气事件 (用于底栏摘要) */
private data class MonthEvent(val day: Int, val time: String, val name: String)

private const val SITE_LON = 116.3833
private const val SITE_LAT = 39.9
private const val SITE_TZ = 8.0

@OptIn(ExperimentalMaterial3Api::class)
@Composable
fun CalendarScreen() {
    val scope = rememberCoroutineScope()
    val now = remember { Calendar.getInstance() }
    val todayY = now.get(Calendar.YEAR)
    val todayM = now.get(Calendar.MONTH) + 1
    val todayD = now.get(Calendar.DAY_OF_MONTH)

    var year by remember { mutableIntStateOf(todayY) }
    var month by remember { mutableIntStateOf(todayM) }
    var yearInput by remember { mutableStateOf(YearUtil.astroYearToStr(todayY)) }
    var monthInput by remember { mutableStateOf("$todayM") }
    var days by remember { mutableStateOf(emptyList<DayInfo>()) }
    var selectedDay by remember { mutableStateOf<DayInfo?>(null) }
    var showSheet by remember { mutableStateOf(false) }

    var rtsDay by remember { mutableIntStateOf(todayD) }
    var rts by remember { mutableStateOf<DayRTS?>(null) }
    var moonEvents by remember { mutableStateOf(emptyList<MonthEvent>()) }
    var jieQiEvents by remember { mutableStateOf(emptyList<MonthEvent>()) }

    fun loadRTS() {
        scope.launch {
            rts = CalendarRepository.calcDayRTS(year, month, rtsDay,
                SITE_LON, SITE_LAT, SITE_TZ)
        }
    }

    fun collectEvents(monthDays: List<DayInfo>) {
        val moon = mutableListOf<MonthEvent>()
        val jq = mutableListOf<MonthEvent>()
        for (d in monthDays) {
            if (d.yueXiangName.isNotEmpty()) {
                val label = when (d.yueXiangName) {
                    "朔" -> "朔月"
                    "望" -> "望月"
                    else -> d.yueXiangName
                }
                moon += MonthEvent(d.solarDay, extractTime(d.yueXiangTime), label)
            }
            if (d.jieQiName.isNotEmpty()) {
                jq += MonthEvent(d.solarDay, extractTime(d.jieQiTime), d.jieQiName)
            }
        }
        moonEvents = moon
        jieQiEvents = jq
    }

    fun loadMonth() {
        scope.launch {
            val md = CalendarRepository.getMonthData(year, month)
            days = md
            yearInput = YearUtil.astroYearToStr(year)
            monthInput = "$month"
            collectEvents(md)
            if (rtsDay > md.size) rtsDay = 1
            rts = CalendarRepository.calcDayRTS(year, month, rtsDay,
                SITE_LON, SITE_LAT, SITE_TZ)
        }
    }

    fun navigate(dy: Int, dm: Int) {
        var y = year + dy
        var m = month + dm
        while (m <= 0)  { m += 12; y-- }
        while (m > 12) { m -= 12; y++ }
        year = y; month = m; rtsDay = 1
        loadMonth()
    }

    fun goToday() {
        year = todayY; month = todayM; rtsDay = todayD
        loadMonth()
    }

    fun applyInput() {
        val y = YearUtil.yearStrToAstro(yearInput)
        val m = monthInput.toIntOrNull() ?: 0
        if (YearUtil.isAstroYearValid(y) && m in 1..12) {
            year = y; month = m; rtsDay = 1
            loadMonth()
        } else {
            yearInput = YearUtil.astroYearToStr(year)
            monthInput = "$month"
        }
    }

    LaunchedEffect(year, month) { loadMonth() }

    Column(Modifier.fillMaxSize().background(Background)) {
        HeaderSection(year, month, days.firstOrNull())
        NavSection(
            yearInput = yearInput, onYearChange = { yearInput = it },
            monthInput = monthInput, onMonthChange = { monthInput = it },
            onPrevYear = { navigate(-1, 0) },
            onPrevMonth = { navigate(0, -1) },
            onNextMonth = { navigate(0, 1) },
            onNextYear = { navigate(1, 0) },
            onApplyInput = ::applyInput,
            onToday = ::goToday
        )
        YearInfoBar(days.firstOrNull())
        WeekHeader()
        CalendarGrid(
            days = days,
            todayY = todayY, todayM = todayM, todayD = todayD,
            selectedDay = selectedDay,
            onTapDay = { d ->
                selectedDay = d
                rtsDay = d.solarDay
                loadRTS()
                showSheet = true
            },
            modifier = Modifier.weight(1f)
        )
        BottomInfoBar(year, month, rtsDay, rts, moonEvents, jieQiEvents)
    }

    if (showSheet && selectedDay != null) {
        DayDetailDialog(
            day = selectedDay!!,
            todayY = todayY, todayM = todayM, todayD = todayD,
            onDismiss = {
                showSheet = false
                selectedDay = null
                // 浮层关闭: 当前若为今月 -> 回到今天的日月升降, 否则用 1 号
                val resetDay = if (year == todayY && month == todayM) todayD else 1
                if (resetDay != rtsDay) {
                    rtsDay = resetDay
                    loadRTS()
                }
            }
        )
    }
}

// ─── Header ─────────────────────────────────────────────────────────────

@Composable
private fun HeaderSection(year: Int, month: Int, first: DayInfo?) {
    Box(
        Modifier
            .fillMaxWidth()
            .background(Brush.horizontalGradient(listOf(GradientStart, GradientEnd)))
            .padding(horizontal = Dimens.paddingLg, vertical = Dimens.paddingSm)
            .padding(top = Dimens.paddingMd)
    ) {
        Row(verticalAlignment = Alignment.CenterVertically) {
            Text(
                "${YearUtil.astroYearToStr(year)}年${month}月",
                fontSize = Dimens.fontTitle, fontWeight = FontWeight.Bold,
                color = OnPrimary
            )
            Spacer(Modifier.weight(1f))
            if (first != null) {
                Text("农历", fontSize = Dimens.fontCaption,
                    color = OnPrimary.copy(alpha = 0.7f))
                Spacer(Modifier.width(4.dp))
                Text("${first.yearGz}年",
                    fontSize = Dimens.fontCaption, fontWeight = FontWeight.Medium,
                    color = Accent)
            }
        }
    }
}

// ─── Nav ────────────────────────────────────────────────────────────────

@Composable
private fun NavSection(
    yearInput: String, onYearChange: (String) -> Unit,
    monthInput: String, onMonthChange: (String) -> Unit,
    onPrevYear: () -> Unit, onPrevMonth: () -> Unit,
    onNextMonth: () -> Unit, onNextYear: () -> Unit,
    onApplyInput: () -> Unit, onToday: () -> Unit
) {
    Row(
        Modifier
            .fillMaxWidth()
            .background(Brush.horizontalGradient(listOf(GradientStart, GradientEnd)))
            .padding(horizontal = Dimens.paddingLg, vertical = Dimens.paddingSm),
        verticalAlignment = Alignment.CenterVertically
    ) {
        NavIconButton("«", onPrevYear)
        NavIconButton("‹", onPrevMonth)
        Spacer(Modifier.width(6.dp))

        TextInputCompact(
            value = yearInput, onValueChange = onYearChange,
            onSubmit = onApplyInput,
            width = 78.dp,
            placeholder = "YYYY/B212"
        )
        Text("年", fontSize = Dimens.fontCaption,
            color = OnPrimary, modifier = Modifier.padding(horizontal = 2.dp))
        TextInputCompact(
            value = monthInput, onValueChange = onMonthChange,
            onSubmit = onApplyInput,
            width = 40.dp,
            keyboardType = KeyboardType.Number,
            placeholder = "M"
        )
        Text("月", fontSize = Dimens.fontCaption,
            color = OnPrimary, modifier = Modifier.padding(horizontal = 2.dp))

        Spacer(Modifier.width(6.dp))
        NavIconButton("›", onNextMonth)
        NavIconButton("»", onNextYear)

        Spacer(Modifier.weight(1f))

        Button(
            onClick = onToday,
            modifier = Modifier.height(30.dp),
            shape = RoundedCornerShape(Dimens.radiusLg),
            colors = ButtonDefaults.buttonColors(
                containerColor = Accent, contentColor = Primary
            ),
            contentPadding = PaddingValues(horizontal = 12.dp, vertical = 0.dp)
        ) {
            Text("今天", fontSize = Dimens.fontCaption, fontWeight = FontWeight.Medium)
        }
    }
}

@Composable
private fun NavIconButton(icon: String, onClick: () -> Unit) {
    Box(
        Modifier
            .padding(end = 4.dp)
            .size(30.dp)
            .clip(CircleShape)
            .background(PrimaryLight.copy(alpha = 0.5f))
            .clickable(onClick = onClick),
        contentAlignment = Alignment.Center
    ) {
        Text(icon, fontSize = 18.sp, fontWeight = FontWeight.Bold, color = OnPrimary)
    }
}

@Composable
private fun TextInputCompact(
    value: String,
    onValueChange: (String) -> Unit,
    onSubmit: () -> Unit,
    width: Dp,
    keyboardType: KeyboardType = KeyboardType.Text,
    placeholder: String = ""
) {
    Box(
        Modifier
            .width(width).height(32.dp)
            .clip(RoundedCornerShape(Dimens.radiusSm))
            .background(PrimaryLight.copy(alpha = 0.4f)),
        contentAlignment = Alignment.Center
    ) {
        if (value.isEmpty() && placeholder.isNotEmpty()) {
            Text(placeholder, color = OnPrimary.copy(alpha = 0.4f),
                fontSize = Dimens.fontBody, textAlign = TextAlign.Center)
        }
        BasicTextField(
            value = value, onValueChange = onValueChange,
            textStyle = TextStyle(
                color = OnPrimary, fontSize = Dimens.fontBody,
                textAlign = TextAlign.Center
            ),
            singleLine = true,
            modifier = Modifier.fillMaxSize().padding(horizontal = 4.dp),
            keyboardOptions = KeyboardOptions(
                keyboardType = keyboardType, imeAction = ImeAction.Done
            ),
            keyboardActions = androidx.compose.foundation.text.KeyboardActions(
                onDone = { onSubmit() }
            ),
            cursorBrush = androidx.compose.ui.graphics.SolidColor(OnPrimary)
        )
    }
}

// ─── Year info bar ──────────────────────────────────────────────────────

@Composable
private fun YearInfoBar(first: DayInfo?) {
    Row(
        Modifier
            .fillMaxWidth()
            .background(Surface)
            .padding(horizontal = Dimens.paddingLg, vertical = Dimens.paddingSm)
            .border(0.5.dp, DividerColor),
        verticalAlignment = Alignment.CenterVertically
    ) {
        Text("☰", fontSize = 14.sp, color = Accent,
            modifier = Modifier.padding(end = 6.dp))
        if (first != null) {
            Text("${first.yearGz}年", fontSize = Dimens.fontBody,
                fontWeight = FontWeight.Medium, color = Primary)
            Text("  |  ", fontSize = Dimens.fontBody, color = DividerColor)
            Text("生肖${first.shengXiao}", fontSize = Dimens.fontBody, color = Primary)
            Spacer(Modifier.weight(1f))
            if (first.huangdiYear > 0) {
                Text("黄帝${first.huangdiYear}年",
                    fontSize = Dimens.fontCaption, color = TextSecondary)
            }
        } else {
            Spacer(Modifier.weight(1f))
        }
    }
}

// ─── Week header ────────────────────────────────────────────────────────

@Composable
private fun WeekHeader() {
    Row(
        Modifier
            .fillMaxWidth()
            .background(Surface)
            .padding(horizontal = Dimens.paddingSm, vertical = Dimens.paddingSm)
    ) {
        WeekNames.forEachIndexed { i, name ->
            Text(
                name,
                modifier = Modifier.weight(1f),
                textAlign = TextAlign.Center,
                fontSize = Dimens.fontCaption,
                fontWeight = FontWeight.Medium,
                color = if (i == 0 || i == 6) WeekendColor else TextSecondary
            )
        }
    }
}

// ─── Calendar grid ──────────────────────────────────────────────────────

@Composable
private fun CalendarGrid(
    days: List<DayInfo>,
    todayY: Int, todayM: Int, todayD: Int,
    selectedDay: DayInfo?,
    onTapDay: (DayInfo) -> Unit,
    modifier: Modifier = Modifier
) {
    val offset = if (days.isNotEmpty()) days[0].weekDay else 0
    val cells: List<DayInfo?> = List(offset) { null } + days

    LazyVerticalGrid(
        columns = GridCells.Fixed(7),
        modifier = modifier.background(Surface).padding(horizontal = Dimens.paddingSm),
    ) {
        items(cells.size) { idx ->
            val d = cells[idx]
            if (d == null) {
                Box(Modifier.height(Dimens.calendarCellSize))
            } else {
                val sel = selectedDay
                val isOpen = sel != null &&
                        sel.solarYear == d.solarYear &&
                        sel.solarMonth == d.solarMonth &&
                        sel.solarDay == d.solarDay
                DayCell(d, todayY, todayM, todayD,
                    isOpen = isOpen,
                    onTap = { onTapDay(d) })
            }
        }
    }
}

@Composable
private fun DayCell(
    day: DayInfo, todayY: Int, todayM: Int, todayD: Int,
    isOpen: Boolean,
    onTap: () -> Unit
) {
    val isToday = day.solarYear == todayY && day.solarMonth == todayM && day.solarDay == todayD
    val bg = when {
        isToday -> Primary
        isOpen -> Accent.copy(alpha = 0.15f)
        else -> Color.Transparent
    }
    Column(
        Modifier
            .height(Dimens.calendarCellSize)
            .padding(2.dp)
            .clip(RoundedCornerShape(Dimens.radiusSm))
            .background(bg)
            .then(
                if (isOpen && !isToday)
                    Modifier.border(1.5.dp, Accent, RoundedCornerShape(Dimens.radiusSm))
                else Modifier
            )
            .clickable(onClick = onTap),
        horizontalAlignment = Alignment.CenterHorizontally,
        verticalArrangement = Arrangement.Center
    ) {
        Row(verticalAlignment = Alignment.CenterVertically) {
            Text(
                "${day.solarDay}",
                fontSize = Dimens.fontSubtitle,
                fontWeight = if (isToday) FontWeight.Bold else FontWeight.Normal,
                color = solarColor(day, isToday)
            )
            MoonIcon(day.yueXiangName)
        }
        Text(
            subText(day),
            fontSize = Dimens.fontSmall,
            fontWeight = if (day.jieQiName.isNotEmpty() || day.holiday.isNotEmpty())
                FontWeight.Medium else FontWeight.Normal,
            color = lunarColor(day, isToday),
            maxLines = 1
        )
    }
}

@Composable
private fun MoonIcon(yueXiangName: String) {
    if (yueXiangName.isEmpty()) return
    val (txt, color) = when (yueXiangName) {
        "朔" -> "●" to Color(0xFF505050)
        "望" -> "●" to Color(0xFFF0B000)
        "上弦" -> "◑" to Color(0xFFF0B000)
        "下弦" -> "◐" to Color(0xFFF0B000)
        else -> return
    }
    Text(txt, fontSize = 8.sp, color = color,
        modifier = Modifier.padding(start = 2.dp))
}

private fun solarColor(day: DayInfo, isToday: Boolean): Color {
    if (isToday) return OnPrimary
    if (day.isOffDay || day.holiday.isNotEmpty()) return WeekendColor
    if (day.weekDay == 0 || day.weekDay == 6) return WeekendColor
    return OnSurface
}

private fun lunarColor(day: DayInfo, isToday: Boolean): Color {
    if (day.holiday.isNotEmpty()) return WeekendColor
    if (day.jieQiName.isNotEmpty()) return JieQiColor
    if (isToday) return OnPrimary.copy(alpha = 0.8f)
    return LunarText
}

private fun subText(day: DayInfo): String {
    if (day.holiday.isNotEmpty()) return shortName(day.holiday)
    if (day.major.isNotEmpty()) return shortName(day.major)
    if (day.jieQiName.isNotEmpty()) return day.jieQiName
    if (day.lunarDayName == "初一") return day.lunarMonthName
    return day.lunarDayName
}

private fun shortName(s: String): String {
    val first = s.split(' ', ',', '，').firstOrNull().orEmpty()
    return if (first.length > 4) first.substring(0, 4) else first
}

// ─── Bottom info bar ────────────────────────────────────────────────────

@Composable
private fun BottomInfoBar(
    year: Int, month: Int, rtsDay: Int, rts: DayRTS?,
    moonEvents: List<MonthEvent>, jieQiEvents: List<MonthEvent>
) {
    Column(
        Modifier
            .fillMaxWidth()
            .background(Surface)
            .padding(Dimens.paddingMd)
    ) {
        Row(verticalAlignment = Alignment.CenterVertically) {
            Text("${YearUtil.astroYearToStr(year)}年${month}月${rtsDay}日 · 日月升降",
                fontSize = Dimens.fontCaption, fontWeight = FontWeight.Medium,
                color = Primary)
            Spacer(Modifier.weight(1f))
            Text(String.format("北京 (%.2f°E %.2f°N)", SITE_LON, SITE_LAT),
                fontSize = Dimens.fontSmall, color = TextSecondary)
        }
        Spacer(Modifier.height(6.dp))
        RTSRow3(
            "日出", rts?.sunRise ?: "--:--:--",
            "日落", rts?.sunSet ?: "--:--:--",
            "中天", rts?.sunMeridian ?: "--:--:--"
        )
        RTSRow3(
            "月出", rts?.moonRise ?: "--:--:--",
            "月落", rts?.moonSet ?: "--:--:--",
            "月中", rts?.moonMeridian ?: "--:--:--"
        )
        RTSRow2(
            "晨起天亮", rts?.civilDawn ?: "--:--:--",
            "晚上天黑", rts?.civilDusk ?: "--:--:--"
        )
        RTSRow2(
            "日照时间", rts?.dayLength ?: "--:--:--",
            "白天时间", rts?.lightLength ?: "--:--:--"
        )

        if (moonEvents.isNotEmpty() || jieQiEvents.isNotEmpty()) {
            Spacer(Modifier.height(8.dp))
            HorizontalDivider(color = DividerColor)
            Spacer(Modifier.height(6.dp))
            Text("${month}月月相与节气",
                fontSize = Dimens.fontSmall, color = TextSecondary)
            Spacer(Modifier.height(4.dp))
            EventGrid(moonEvents, PrimaryLight)
            if (jieQiEvents.isNotEmpty()) {
                EventGrid(jieQiEvents, JieQiColor)
            }
        }
    }
}

@Composable
private fun RTSRow3(l1: String, v1: String, l2: String, v2: String,
                    l3: String, v3: String) {
    Row(Modifier.fillMaxWidth().padding(vertical = 2.dp)) {
        RTSLabeled(l1, v1, Modifier.weight(1f))
        RTSLabeled(l2, v2, Modifier.weight(1f))
        RTSLabeled(l3, v3, Modifier.weight(1f))
    }
}

@Composable
private fun RTSRow2(l1: String, v1: String, l2: String, v2: String) {
    Row(Modifier.fillMaxWidth().padding(vertical = 2.dp)) {
        RTSLabeled(l1, v1, Modifier.weight(1f))
        RTSLabeled(l2, v2, Modifier.weight(1f))
    }
}

@Composable
private fun RTSLabeled(label: String, value: String, modifier: Modifier = Modifier) {
    Row(modifier, verticalAlignment = Alignment.CenterVertically) {
        Text(label, fontSize = Dimens.fontSmall, color = TextSecondary,
            modifier = Modifier.padding(end = 4.dp))
        Text(value, fontSize = Dimens.fontCaption, fontWeight = FontWeight.Medium,
            color = if (value == "--:--:--") TextSecondary else OnSurface)
    }
}

@Composable
private fun EventGrid(events: List<MonthEvent>, color: Color) {
    // 两列布局
    Column(Modifier.fillMaxWidth()) {
        events.chunked(2).forEach { pair ->
            Row(Modifier.fillMaxWidth()) {
                EventCell(pair[0], color, Modifier.weight(1f))
                if (pair.size > 1) {
                    EventCell(pair[1], color, Modifier.weight(1f))
                } else {
                    Spacer(Modifier.weight(1f))
                }
            }
        }
    }
}

@Composable
private fun EventCell(e: MonthEvent, color: Color, modifier: Modifier) {
    Row(modifier.padding(vertical = 2.dp),
        verticalAlignment = Alignment.CenterVertically) {
        Text("${padNum(e.day, 2)}日",
            fontSize = Dimens.fontSmall, color = TextSecondary,
            modifier = Modifier.width(28.dp))
        Text(e.time,
            fontSize = Dimens.fontSmall, color = OnSurface,
            modifier = Modifier.padding(start = 2.dp, end = 4.dp))
        Text(e.name,
            fontSize = Dimens.fontSmall, fontWeight = FontWeight.Medium, color = color)
    }
}

private fun padNum(n: Int, width: Int): String {
    val s = n.toString()
    return if (s.length >= width) s else " ".repeat(width - s.length) + s
}

private fun extractTime(s: String): String {
    val m = Regex("""\d{2}:\d{2}:\d{2}""").find(s)
    return m?.value ?: s.trim()
}

// ─── Detail Dialog (深色紧凑卡片, 参考 web 版样式) ──────────────────────
//   - 卡片宽 ≤ 320dp, 内边距 12dp, 信息行密集排布
//   - 纯深色背景 (Surface) + 强阴影; 不再用渐变 (避免渲染异常)
//   - 仅展示: 标题 / 历法转换三行 / 三柱纳音 / 元数据 / 事件 / 纪年
// ────────────────────────────────────────────────────────────────────────

private val PopupBg      = Color(0xFF14163A)   // 深午夜蓝, 与 statusBar 协调
private val PopupBorder  = Color(0xFF3B3F73)
private val PopupText    = Color(0xFFFFFFFF)
private val PopupSub     = Color(0xFFDCE0F2)   // 更亮、纯不透明: 避免在 #14163A 上看着发灰
private val PopupGold    = Color(0xFFFFD54F)
private val PopupGreen   = Color(0xFF80E0A7)
private val PopupRed     = Color(0xFFFF9D9D)
private val PopupDivider = Color(0xFF5A5F8C)   // 纯色分隔线 (不再 20% 透明白)

@Composable
private fun DayDetailDialog(
    day: DayInfo,
    todayY: Int, todayM: Int, todayD: Int,
    onDismiss: () -> Unit
) {
    val isToday = day.solarYear == todayY &&
                  day.solarMonth == todayM &&
                  day.solarDay == todayD
    Dialog(
        onDismissRequest = onDismiss,
        properties = DialogProperties(usePlatformDefaultWidth = false)
    ) {
        Box(
            Modifier.fillMaxWidth(),
            contentAlignment = Alignment.Center
        ) {
            Surface(
                modifier = Modifier
                    .widthIn(max = 260.dp)
                    .padding(horizontal = Dimens.paddingLg),
                shape = RoundedCornerShape(Dimens.radiusMd),
                color = PopupBg,
                tonalElevation = 8.dp,
                shadowElevation = 12.dp,
                border = androidx.compose.foundation.BorderStroke(0.5.dp, PopupBorder)
            ) {
                Column(
                    Modifier
                        .verticalScroll(rememberScrollState())
                        .padding(horizontal = 14.dp, vertical = 12.dp)
                ) {
                    PopupTitleRow(day, isToday)
                    Spacer(Modifier.height(6.dp))
                    PopupInfoLine(
                        "黄帝${day.huangdiYear}年",
                        day.weekName,
                        day.constellationName
                    )
                    PopupInfoLine(
                        "${day.yearGz}年",
                        (if (day.isLeapMonth) "闰" else "") + day.lunarMonthName,
                        "${day.lunarDayName}日"
                    )
                    PopupInfoLine(
                        "${day.yearGz}年",
                        "${day.monthGz}月",
                        "${day.dayGz}日"
                    )
                    PopupInfoLine(
                        day.yearNaYin, day.monthNaYin, day.dayNaYin,
                        color = PopupGold, bold = true
                    )

                    Spacer(Modifier.height(6.dp))
                    PopupDividerLine()
                    Spacer(Modifier.height(4.dp))

                    PopupCenterRow(
                        listOf(
                            "生肖" to day.shengXiao,
                            "建除" to day.jian12Name
                        )
                    )
                    PopupCenteredText(
                        "回历[${day.moslemYear}年${day.moslemMonth}月${day.moslemDay}日]",
                        color = PopupSub
                    )
                    PopupCenteredText(
                        "JD ${formatJD(day.julianDay)}",
                        color = PopupSub
                    )

                    val hasEvent = day.jieQiName.isNotEmpty() ||
                            day.yueXiangName.isNotEmpty() ||
                            day.holiday.isNotEmpty() ||
                            day.major.isNotEmpty() ||
                            day.minor.isNotEmpty() ||
                            day.misc.isNotEmpty()
                    if (hasEvent) {
                        Spacer(Modifier.height(6.dp))
                        PopupDividerLine()
                        Spacer(Modifier.height(4.dp))
                        if (day.jieQiName.isNotEmpty()) {
                            PopupEventLine("🌿", "节气 ${day.jieQiName}",
                                extractTime(day.jieQiTime), PopupGreen)
                        }
                        if (day.yueXiangName.isNotEmpty()) {
                            PopupEventLine("🌙", "月相 ${day.yueXiangName}",
                                extractTime(day.yueXiangTime), PopupGold)
                        }
                        if (day.holiday.isNotEmpty()) {
                            PopupEventLine("🎉", day.holiday, "", PopupRed)
                        }
                        if (day.major.isNotEmpty()) {
                            PopupEventLine("🎊", day.major, "", PopupGold)
                        }
                        if (day.minor.isNotEmpty()) {
                            PopupEventLine("·", day.minor, "", PopupSub)
                        }
                        if (day.misc.isNotEmpty()) {
                            PopupEventLine("🍂", day.misc, "", PopupGreen)
                        }
                    }
                }
            }
        }
    }
}

@Composable
private fun PopupTitleRow(d: DayInfo, isToday: Boolean) {
    Row(
        Modifier.fillMaxWidth(),
        verticalAlignment = Alignment.CenterVertically
    ) {
        Text(
            "${YearUtil.astroYearToStr(d.solarYear)}年${d.solarMonth}月${d.solarDay}日",
            fontSize = Dimens.fontSubtitle,
            fontWeight = FontWeight.Bold,
            color = PopupText,
            modifier = Modifier.weight(1f),
            textAlign = TextAlign.Center
        )
        if (isToday) {
            Box(
                Modifier
                    .clip(RoundedCornerShape(4.dp))
                    .background(Accent)
                    .padding(horizontal = 5.dp, vertical = 1.dp)
            ) {
                Text("今天", fontSize = Dimens.fontSmall,
                    fontWeight = FontWeight.Medium, color = Primary)
            }
        }
    }
}

@Composable
private fun PopupInfoLine(
    a: String, b: String, c: String,
    color: Color = PopupText,
    bold: Boolean = false
) {
    // 居中紧凑三列: 不再各占 1/3 整宽, 改为水平居中 + 小间隔 (6+6=12dp)
    Row(
        Modifier.fillMaxWidth().padding(vertical = 1.dp),
        horizontalArrangement = Arrangement.Center
    ) {
        PopupCell(a, color, bold)
        PopupCell(b, color, bold)
        PopupCell(c, color, bold)
    }
}

@Composable
private fun PopupCell(text: String, color: Color, bold: Boolean = false) {
    Text(
        text.ifEmpty { "—" },
        modifier = Modifier.padding(horizontal = 6.dp),
        fontSize = Dimens.fontCaption,
        fontWeight = if (bold) FontWeight.Medium else FontWeight.Normal,
        color = if (text.isEmpty()) PopupSub.copy(alpha = 0.4f) else color,
        maxLines = 1,
        textAlign = TextAlign.Center
    )
}

@Composable
private fun PopupCenterRow(items: List<Pair<String, String>>) {
    Row(
        Modifier.fillMaxWidth().padding(vertical = 1.dp),
        horizontalArrangement = Arrangement.Center
    ) {
        items.forEachIndexed { i, (lbl, v) ->
            if (i > 0) {
                Text("·", color = PopupSub,
                    modifier = Modifier.padding(horizontal = 6.dp),
                    fontSize = Dimens.fontCaption)
            }
            Text(lbl, color = PopupSub, fontSize = Dimens.fontSmall,
                modifier = Modifier.padding(end = 4.dp))
            Text(v, color = PopupText, fontSize = Dimens.fontCaption,
                fontWeight = FontWeight.Medium)
        }
    }
}

@Composable
private fun PopupCenteredText(
    text: String,
    color: Color = PopupText,
    size: androidx.compose.ui.unit.TextUnit = Dimens.fontSmall
) {
    Text(
        text,
        modifier = Modifier.fillMaxWidth().padding(vertical = 1.dp),
        textAlign = TextAlign.Center,
        fontSize = size,
        color = color,
        maxLines = 1
    )
}

@Composable
private fun PopupEventLine(icon: String, text: String, time: String, color: Color) {
    Row(
        Modifier.fillMaxWidth().padding(vertical = 1.dp),
        horizontalArrangement = Arrangement.Center,
        verticalAlignment = Alignment.CenterVertically
    ) {
        Text(icon, fontSize = 11.sp, color = color,
            modifier = Modifier.padding(end = 4.dp))
        Text(text, fontSize = Dimens.fontCaption,
            fontWeight = FontWeight.Medium, color = color)
        if (time.isNotEmpty()) {
            Text("  $time", fontSize = Dimens.fontSmall, color = PopupSub)
        }
    }
}

@Composable
private fun PopupDividerLine() {
    Box(
        Modifier
            .fillMaxWidth()
            .height(0.5.dp)
            .background(PopupDivider)
    )
}

private fun formatJD(jd: Double): String {
    if (jd.isNaN()) return "-"
    val full = jd.toLong() + (if (jd - jd.toLong() >= 0.5) 1 else 0)
    val d0 = full - 2451545
    return "$full($d0)"
}
