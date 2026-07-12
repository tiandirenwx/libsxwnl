package com.sxwnl.calendar.ui.screens

import androidx.compose.foundation.background
import androidx.compose.foundation.border
import androidx.compose.foundation.clickable
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.lazy.grid.GridCells
import androidx.compose.foundation.lazy.grid.LazyVerticalGrid
import androidx.compose.foundation.lazy.grid.items
import androidx.compose.foundation.horizontalScroll
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
import com.sxwnl.calendar.data.AlmanacQuote
import com.sxwnl.calendar.data.AlmanacTopic
import com.sxwnl.calendar.data.CalendarRepository
import com.sxwnl.calendar.data.DayAlmanac
import com.sxwnl.calendar.data.DayInfo
import com.sxwnl.calendar.data.DayRTS
import com.sxwnl.calendar.data.EventAdvice
import com.sxwnl.calendar.data.GeoCity
import com.sxwnl.calendar.data.GeoProvince
import com.sxwnl.calendar.data.LuckyHour
import com.sxwnl.calendar.data.ShenSha
import com.sxwnl.calendar.data.TimezoneGroup
import com.sxwnl.calendar.ui.theme.*
import com.sxwnl.calendar.util.YearUtil
import kotlinx.coroutines.delay
import kotlinx.coroutines.launch
import java.util.Calendar
import java.util.TimeZone
import kotlin.math.abs

// ════════════════════════════════════════════════════════════════
//  CalendarScreen — 与鸿蒙端 CalendarPage.ets 对齐的月历页面
// ════════════════════════════════════════════════════════════════

/** 月相 / 节气事件 (用于底栏摘要) */
private data class MonthEvent(val day: Int, val time: String, val name: String)

/** 默认地点 — NAPI 失败时兜底用 (天安门, UTC+8) */
private val DEFAULT_LOCATION = GeoCity("北京市", "天安门", 116.3833, 39.9, 8.0)

@OptIn(ExperimentalMaterial3Api::class)
@Composable
fun CalendarScreen() {
    val scope = rememberCoroutineScope()

    // ── "今天" 用可变状态: app 回到前台时刷新 (跨天/退出再进来回到今天) ──
    var todayY by remember { mutableIntStateOf(Calendar.getInstance().get(Calendar.YEAR)) }
    var todayM by remember { mutableIntStateOf(Calendar.getInstance().get(Calendar.MONTH) + 1) }
    var todayD by remember { mutableIntStateOf(Calendar.getInstance().get(Calendar.DAY_OF_MONTH)) }

    var year by remember { mutableIntStateOf(todayY) }
    var month by remember { mutableIntStateOf(todayM) }
    var yearInput by remember { mutableStateOf(YearUtil.astroYearToStr(todayY)) }
    var monthInput by remember { mutableStateOf("$todayM") }
    var days by remember { mutableStateOf(emptyList<DayInfo>()) }
    var selectedDay by remember { mutableStateOf<DayInfo?>(null) }
    var showSheet by remember { mutableStateOf(false) }

    // 老黄历底部抽屉 (ModalBottomSheet) — 对齐鸿蒙 bindSheet
    var showAlmanacSheet by remember { mutableStateOf(false) }
    var almanacSheetDay by remember { mutableStateOf<DayInfo?>(null) }

    var rtsDay by remember { mutableIntStateOf(todayD) }
    var rts by remember { mutableStateOf<DayRTS?>(null) }
    var moonEvents by remember { mutableStateOf(emptyList<MonthEvent>()) }
    var jieQiEvents by remember { mutableStateOf(emptyList<MonthEvent>()) }

    // ── 地点选择 (对齐鸿蒙端 CalendarPage.ets) ─────────────────
    //  国内: 改 location → 触发日月升降重算
    //  国际: 仅供顶栏"外地实时时间"展示, 不修改 location
    var location by remember { mutableStateOf(DEFAULT_LOCATION) }
    var provinces by remember { mutableStateOf(emptyList<GeoProvince>()) }
    var cities by remember { mutableStateOf(emptyList<GeoCity>()) }
    var provIdx by remember { mutableIntStateOf(0) }
    var cityIdx by remember { mutableIntStateOf(0) }

    var continents by remember { mutableStateOf(emptyList<String>()) }
    var countries by remember { mutableStateOf(emptyList<TimezoneGroup>()) }
    var contIdx by remember { mutableIntStateOf(0) }
    var countryIdx by remember { mutableIntStateOf(0) }
    var allTzGroups by remember { mutableStateOf(emptyList<TimezoneGroup>()) }

    // 双时钟 (本地系统时区 + 选中国际时区), 每秒刷新
    var localClock by remember { mutableStateOf("") }
    var intlClock by remember { mutableStateOf("") }

    // 老黄历静态知识 (董公总论/口诀/方位 — 全局只取一次, 懒加载常驻)
    var topics by remember { mutableStateOf(emptyList<AlmanacTopic>()) }
    var showTopics by remember { mutableStateOf(false) }

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

    fun resetToToday() {
        // 用全新 Calendar 读取, 避免 app 常驻跨天后"今天"陈旧
        val n = Calendar.getInstance()
        todayY = n.get(Calendar.YEAR)
        todayM = n.get(Calendar.MONTH) + 1
        todayD = n.get(Calendar.DAY_OF_MONTH)
        year = todayY; month = todayM; rtsDay = todayD
    }

    fun navigate(dy: Int, dm: Int) {
        var y = year + dy
        var m = month + dm
        while (m <= 0)  { m += 12; y-- }
        while (m > 12) { m -= 12; y++ }
        year = y; month = m; rtsDay = 1
        // 不显式 loadMonth — 依赖 LaunchedEffect(year, month, location) 统一驱动,
        // 避免与 LaunchedEffect 并行触发双重 IO 调用.
    }

    fun goToday() = resetToToday()

    fun applyInput() {
        val y = YearUtil.yearStrToAstro(yearInput)
        val m = monthInput.toIntOrNull() ?: 0
        if (YearUtil.isAstroYearValid(y) && m in 1..12) {
            year = y; month = m; rtsDay = 1
        } else {
            yearInput = YearUtil.astroYearToStr(year)
            monthInput = "$month"
        }
    }

    // ── 一次性初始化: 加载省/市/时区目录 + 默认地点 (单次, Unit key) ──
    //   不在这里加载月数据 — 让下面的 LaunchedEffect(year, month, location) 统一处理.
    //   这样消除"两个 effect 并行触发 loadMonth 产生竞态 (谁后写 rts 谁赢)"的问题.
    LaunchedEffect(Unit) {
        val ps = CalendarRepository.geoListProvinces()
        if (ps.isNotEmpty()) {
            provinces = ps
            val def = CalendarRepository.geoDefault() ?: DEFAULT_LOCATION
            val pi = ps.indexOfFirst { it.province == def.province }.coerceAtLeast(0)
            provIdx = pi
            val cs = CalendarRepository.geoListCities(ps[pi].province)
            cities = cs
            val ci = cs.indexOfFirst { it.district == def.district }.coerceAtLeast(0)
            cityIdx = ci
            location = if (cs.isNotEmpty()) cs[ci] else def
        }
        val tz = CalendarRepository.geoListTimezones()
        if (tz.isNotEmpty()) {
            allTzGroups = tz
            val seen = HashSet<String>()
            continents = tz.mapNotNull { if (seen.add(it.continent)) it.continent else null }
            if (continents.isNotEmpty()) {
                val asiaIdx = continents.indexOf("亚洲").coerceAtLeast(0)
                contIdx = asiaIdx
                val countriesIn = tz.filter { it.continent == continents[asiaIdx] }
                countries = countriesIn
                val chinaIdx = countriesIn.indexOfFirst { it.country == "中国" }
                    .coerceAtLeast(0)
                countryIdx = chinaIdx
            }
        }
    }

    // ── 时钟: 每秒刷新本地 + 选中国际时区 ──
    //
    //  关键点: 显示某时区的"当地小时"必须用对应 TimeZone 的 Calendar 来读字段,
    //  不能直接 timeInMillis += offset + 用本机时区 Calendar 读 — 那样读出的是
    //  设备本地时区下对应那一刻的时分秒, 而不是目标时区的本地时间.
    //
    //  做法: 共用同一个 timeInMillis (绝对 UTC 时刻), 切换 Calendar 的 timeZone,
    //  即可正确显示各自时区的当地时间.
    //
    //  LaunchedEffect 在 Composable 离开组合时自动取消, 协程随之停止, 不会
    //  泄漏定时器. 切换 countries/countryIdx 会重启循环以重算外地时间.
    LaunchedEffect(countries, countryIdx) {
        while (true) {
            val nowMs = System.currentTimeMillis()
            val localCal = Calendar.getInstance().apply { timeInMillis = nowMs }
            localClock = fmtClock(localCal)
            if (countries.isNotEmpty() && countryIdx in countries.indices) {
                val tz = countries[countryIdx].timezone
                // 用 fixed-offset TimeZone 读取, 避免设备本地时区干扰
                val intlTz = TimeZone.getTimeZone(buildGmtId(tz))
                val intlCal = Calendar.getInstance(intlTz).apply { timeInMillis = nowMs }
                intlClock = fmtClock(intlCal)
            } else {
                intlClock = ""
            }
            delay(1000)
        }
    }

    // 统一月数据加载: 任何 year/month 改变都重新拉, 用最新 location 计算 rts.
    LaunchedEffect(year, month) {
        val md = CalendarRepository.getMonthData(year, month)
        days = md
        yearInput = YearUtil.astroYearToStr(year)
        monthInput = "$month"
        collectEvents(md)
        if (rtsDay > md.size) rtsDay = 1
    }

    // 统一 RTS 加载: 任何 (year, month, rtsDay, location) 改变都重新算.
    //   - 切换日期(navigate / goToday / applyInput) → 触发
    //   - 切换城市 (onCityChange / onProvinceChange) → 触发
    //   - 点击某日 (onTapDay 设置 rtsDay) → 触发
    //   消除分散的手动 loadRTS 调用, 避免与 effect 重复.
    LaunchedEffect(year, month, rtsDay, location) {
        rts = CalendarRepository.calcDayRTS(year, month, rtsDay,
            location.longitude, location.latitude, location.timezone)
    }

    // ── app 回到前台 → 回到今天 (与鸿蒙"退出再进来回到今天"对齐) ──
    //   ON_RESUME 在首次启动也会触发 (此时已是今天, resetToToday 为空操作);
    //   切 Tab 不会触发 (Activity 生命周期不变), 不影响 Tab 间浏览状态.
    val lifecycleOwner = androidx.compose.ui.platform.LocalLifecycleOwner.current
    DisposableEffect(lifecycleOwner) {
        val observer = androidx.lifecycle.LifecycleEventObserver { _, event ->
            if (event == androidx.lifecycle.Lifecycle.Event.ON_RESUME) {
                resetToToday()
            }
        }
        lifecycleOwner.lifecycle.addObserver(observer)
        onDispose { lifecycleOwner.lifecycle.removeObserver(observer) }
    }

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
        TopLocationStrip(
            // ── 国际 ──
            continents = continents, contIdx = contIdx,
            countries = countries, countryIdx = countryIdx,
            intlClock = intlClock,
            localClock = localClock,
            onContinentChange = { idx ->
                if (idx in continents.indices) {
                    contIdx = idx
                    val filtered = allTzGroups.filter { it.continent == continents[idx] }
                    countries = filtered
                    countryIdx = 0
                }
            },
            onCountryChange = { idx ->
                if (idx in countries.indices) countryIdx = idx
            },
            // ── 国内 ──
            provinces = provinces, provIdx = provIdx,
            cities = cities, cityIdx = cityIdx, location = location,
            onProvinceChange = { idx ->
                if (idx in provinces.indices) {
                    provIdx = idx
                    scope.launch {
                        val cs = CalendarRepository.geoListCities(provinces[idx].province)
                        cities = cs
                        cityIdx = 0
                        // location 改变会自动触发 LaunchedEffect(...) 重算 rts
                        if (cs.isNotEmpty()) location = cs[0]
                    }
                }
            },
            onCityChange = { idx ->
                if (idx in cities.indices) {
                    cityIdx = idx
                    location = cities[idx]
                }
            }
        )
        YearInfoBar(days.firstOrNull())
        WeekHeader()
        CalendarGrid(
            days = days,
            todayY = todayY, todayM = todayM, todayD = todayD,
            selectedDay = selectedDay,
            onTapDay = { d ->
                selectedDay = d
                rtsDay = d.solarDay   // 触发 LaunchedEffect 重算 rts
                showSheet = true
            },
            modifier = Modifier.weight(1f)
        )
        BottomInfoBar(month, rtsDay, rts, moonEvents, jieQiEvents, localClock)
    }

    if (showSheet && selectedDay != null) {
        DayDetailDialog(
            day = selectedDay!!,
            todayY = todayY, todayM = todayM, todayD = todayD,
            onDismiss = {
                showSheet = false
                selectedDay = null
                // 浮层关闭: 今月回到今天的日月升降, 否则用 1 号; rtsDay 变化会
                // 自动触发 LaunchedEffect 重算
                val resetDay = if (year == todayY && month == todayM) todayD else 1
                if (resetDay != rtsDay) rtsDay = resetDay
            },
            onOpenAlmanacSheet = {
                // 切换到"老黄历底部抽屉": 关闭日详情弹窗, 打开 ModalBottomSheet.
                //   与鸿蒙 openAlmanacSheet -> bindSheet 语义一致 (点空白/下拉即退出).
                showSheet = false
                almanacSheetDay = selectedDay
                showAlmanacSheet = true
            }
        )
    }

    // 老黄历底部抽屉 (ModalBottomSheet) — 点空白或下拉即关闭, 不再需要后退键
    if (showAlmanacSheet && almanacSheetDay != null) {
        AlmanacSheetBottomSheet(
            day = almanacSheetDay!!,
            onDismiss = {
                showAlmanacSheet = false
                almanacSheetDay = null
            },
            onOpenTopics = {
                showAlmanacSheet = false
                showTopics = true
            }
        )
    }

    // 静态知识 Dialog (按需懒加载 + 缓存)
    if (showTopics) {
        LaunchedEffect(Unit) {
            if (topics.isEmpty()) topics = CalendarRepository.getAlmanacTopics()
        }
        AlmanacTopicsDialog(
            topics = topics,
            onDismiss = { showTopics = false }
        )
    }
}

// 双时钟显示格式 "M/d HH:mm:ss"
//  注意: 字段 (HOUR_OF_DAY 等) 是按 Calendar 实例的 timeZone 解读 timeInMillis,
//  调用方必须保证 Calendar 设置了正确的目标时区, 否则会用设备本地时区, 显示错误.
private fun fmtClock(c: Calendar): String {
    val m = c.get(Calendar.MONTH) + 1
    val d = c.get(Calendar.DAY_OF_MONTH)
    val h = c.get(Calendar.HOUR_OF_DAY)
    val mi = c.get(Calendar.MINUTE)
    val s = c.get(Calendar.SECOND)
    return "%d/%d %02d:%02d:%02d".format(m, d, h, mi, s)
}

/**
 * 把 sxwnl 的"小时时区"转成 Java TimeZone ID. 例:
 *   8.0  → "GMT+08:00" (北京)
 *   -5.0 → "GMT-05:00" (纽约 EST)
 *   5.5  → "GMT+05:30" (印度)
 *
 *  使用 "GMT±HH:MM" 而非 "Etc/GMT±H" 的原因:
 *  - "Etc/GMT" 系列的符号是反的 ("Etc/GMT-8" 才是 UTC+8), 极易写错
 *  - "GMT±HH:MM" 是 ISO 8601 风格, Java TimeZone 支持得很好
 */
private fun buildGmtId(tz: Double): String {
    val sign = if (tz >= 0) "+" else "-"
    val abs = abs(tz)
    val h = abs.toInt()
    val mm = ((abs - h) * 60).toInt()
    return "GMT%s%02d:%02d".format(sign, h, mm)
}

private fun fmtTz(tz: Double): String {
    val sign = if (tz >= 0) "+" else "-"
    val abs = abs(tz)
    val h = abs.toInt()
    val mm = ((abs - h) * 60).toInt()
    return if (mm == 0) "UTC$sign$h" else "UTC$sign$h:%02d".format(mm)
}

private fun fmtLonLat(lon: Double, lat: Double): String {
    val lonDir = if (lon >= 0) "E" else "W"
    val latDir = if (lat >= 0) "N" else "S"
    return "%.1f°%s %.1f°%s".format(abs(lon), lonDir, abs(lat), latDir)
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
                Text("${first.yearGZ}年",
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

// ─── Top location strip ────────────────────────────────────────────────
//
//  两行内联下拉, 对齐 JS 原版 indexmp.htm 顶部 (鸿蒙 CalendarPage.ets):
//
//    第 1 行: [洲 ▼] [国家 ▼]  · UTC±x · 当地实时时钟   (国际, 仅展示)
//    第 2 行: [省 ▼] [市/县 ▼] · 经X°E 纬Y°N             (国内, 改 location)
//
//  · 国内行决定日月升降算用的经纬度
//  · 国际行只影响顶栏外地时钟; 不改 location, 与 JS 一致
// ────────────────────────────────────────────────────────────────────────

@OptIn(ExperimentalMaterial3Api::class)
@Composable
private fun TopLocationStrip(
    continents: List<String>, contIdx: Int,
    countries: List<TimezoneGroup>, countryIdx: Int,
    intlClock: String,
    localClock: String,
    onContinentChange: (Int) -> Unit,
    onCountryChange: (Int) -> Unit,
    provinces: List<GeoProvince>, provIdx: Int,
    cities: List<GeoCity>, cityIdx: Int,
    location: GeoCity,
    onProvinceChange: (Int) -> Unit,
    onCityChange: (Int) -> Unit
) {
    Column(
        Modifier.fillMaxWidth().background(Surface)
            .border(0.5.dp, DividerColor)
    ) {
        // ── 国际行 ──
        Row(
            Modifier.fillMaxWidth().padding(horizontal = 8.dp, vertical = 4.dp),
            verticalAlignment = Alignment.CenterVertically
        ) {
            if (continents.isNotEmpty()) {
                DropdownText(
                    label = continents.getOrNull(contIdx) ?: "",
                    options = continents,
                    onSelect = onContinentChange
                )
                Spacer(Modifier.width(6.dp))
                DropdownText(
                    label = countries.getOrNull(countryIdx)?.country ?: "",
                    options = countries.map { it.country },
                    onSelect = onCountryChange
                )
            }
            Spacer(Modifier.weight(1f))
            val tz = countries.getOrNull(countryIdx)?.timezone
            if (tz != null) {
                Text(fmtTz(tz), fontSize = Dimens.fontSmall, color = Primary,
                    modifier = Modifier.padding(horizontal = 6.dp))
            }
            Text(intlClock, fontSize = Dimens.fontSmall,
                fontWeight = FontWeight.Medium, color = Accent)
        }
        // ── 国内行 ──
        Row(
            Modifier.fillMaxWidth().padding(horizontal = 8.dp, vertical = 4.dp),
            verticalAlignment = Alignment.CenterVertically
        ) {
            if (provinces.isNotEmpty()) {
                DropdownText(
                    label = provinces.getOrNull(provIdx)?.province ?: "",
                    options = provinces.map { it.province },
                    onSelect = onProvinceChange
                )
                Spacer(Modifier.width(6.dp))
                DropdownText(
                    label = cities.getOrNull(cityIdx)?.district ?: "",
                    options = cities.map { it.district },
                    onSelect = onCityChange,
                    primary = true
                )
            }
            Spacer(Modifier.weight(1f))
            // 国内行右侧: 经纬度 + 本地系统时钟 (与国际行 intlClock 形成对照)
            Text(fmtLonLat(location.longitude, location.latitude),
                fontSize = Dimens.fontSmall, color = TextSecondary,
                modifier = Modifier.padding(end = 6.dp))
            Text(localClock, fontSize = Dimens.fontSmall,
                fontWeight = FontWeight.Medium, color = Primary)
        }
    }
}

/** 紧凑下拉文本 — 单击展开列表; options 大于 200 项时建议改用搜索 sheet */
@Composable
private fun DropdownText(
    label: String,
    options: List<String>,
    onSelect: (Int) -> Unit,
    primary: Boolean = false
) {
    var expanded by remember { mutableStateOf(false) }
    Box {
        Row(
            Modifier
                .clip(RoundedCornerShape(Dimens.radiusSm))
                .background(PrimaryLight.copy(alpha = 0.1f))
                .clickable { expanded = true }
                .padding(horizontal = 8.dp, vertical = 4.dp),
            verticalAlignment = Alignment.CenterVertically
        ) {
            Text(
                label.ifEmpty { "—" },
                fontSize = Dimens.fontSmall,
                fontWeight = FontWeight.Medium,
                color = if (primary) Primary else OnSurface,
                maxLines = 1
            )
            Text(" ▼", fontSize = 8.sp, color = TextSecondary)
        }
        DropdownMenu(
            expanded = expanded,
            onDismissRequest = { expanded = false },
            modifier = Modifier
                .heightIn(max = 320.dp)
                .background(Surface)
        ) {
            options.forEachIndexed { idx, opt ->
                DropdownMenuItem(
                    text = {
                        Text(opt, fontSize = Dimens.fontCaption, color = OnSurface)
                    },
                    onClick = {
                        expanded = false
                        onSelect(idx)
                    }
                )
            }
        }
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
            Text("${first.yearGZ}年", fontSize = Dimens.fontBody,
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
//
//  与鸿蒙 CalendarPage.bottomInfoBar / sunMoonCompact 对齐: 紧凑布局, 不占用
//  过多高度, 保证月历主体完整可见无需拖动.
//
//  · 标题行: "{月}月{日}日 · 日月升降"  +  本地系统时钟
//  · 一行 6 列: ☀↑ / ☀ / ☀↓ / ☾↑ / ☾ / ☾↓
//  · 一行 4 列: 晨 / 日 / 昏 / 白
//  · 月相 / 节气 — 流式排列 (一行排不下自动折行)
// ────────────────────────────────────────────────────────────────────────

@OptIn(androidx.compose.foundation.layout.ExperimentalLayoutApi::class)
@Composable
private fun BottomInfoBar(
    month: Int, rtsDay: Int, rts: DayRTS?,
    moonEvents: List<MonthEvent>, jieQiEvents: List<MonthEvent>,
    localClock: String
) {
    Column(
        Modifier
            .fillMaxWidth()
            .background(Surface)
            .border(0.5.dp, DividerColor)
            .padding(horizontal = Dimens.paddingSm, vertical = 6.dp)
    ) {
        // 标题行
        Row(Modifier.fillMaxWidth(), verticalAlignment = Alignment.CenterVertically) {
            Text("${month}月${rtsDay}日 · 日月升降",
                fontSize = Dimens.fontSmall, fontWeight = FontWeight.Medium,
                color = Primary)
            Spacer(Modifier.weight(1f))
            Text(localClock, fontSize = Dimens.fontSmall, color = TextSecondary)
        }
        Spacer(Modifier.height(2.dp))

        // 一行 6 列: 太阳出/中/落 + 月亮出/中/落
        Row(Modifier.fillMaxWidth()) {
            TinyCell("☀↑", rts?.sunRise     ?: "--:--:--", Modifier.weight(1f))
            TinyCell("☀",  rts?.sunMeridian ?: "--:--:--", Modifier.weight(1f))
            TinyCell("☀↓", rts?.sunSet      ?: "--:--:--", Modifier.weight(1f))
            TinyCell("☾↑", rts?.moonRise    ?: "--:--:--", Modifier.weight(1f))
            TinyCell("☾",  rts?.moonMeridian?: "--:--:--", Modifier.weight(1f))
            TinyCell("☾↓", rts?.moonSet     ?: "--:--:--", Modifier.weight(1f))
        }
        Spacer(Modifier.height(2.dp))

        // 一行 4 列: 晨/日/昏/白
        Row(Modifier.fillMaxWidth()) {
            TinyCell("晨", rts?.civilDawn   ?: "--:--:--", Modifier.weight(1f))
            TinyCell("日", rts?.dayLength   ?: "--:--:--", Modifier.weight(1f))
            TinyCell("昏", rts?.civilDusk   ?: "--:--:--", Modifier.weight(1f))
            TinyCell("白", rts?.lightLength ?: "--:--:--", Modifier.weight(1f))
        }

        // 月相 / 节气 — 流式
        if (moonEvents.isNotEmpty()) {
            EventFlowLine("月相", moonEvents, month, PrimaryLight)
        }
        if (jieQiEvents.isNotEmpty()) {
            EventFlowLine("节气", jieQiEvents, month, JieQiColor)
        }
    }
}

/** 6/4 列等宽紧凑单元 (label + value), 居中排布 */
@Composable
private fun TinyCell(label: String, value: String, modifier: Modifier = Modifier) {
    Row(
        modifier,
        horizontalArrangement = Arrangement.Center,
        verticalAlignment = Alignment.CenterVertically
    ) {
        Text(label, fontSize = Dimens.fontSmall, color = TextSecondary,
            modifier = Modifier.padding(end = 3.dp))
        Text(value, fontSize = Dimens.fontSmall, fontWeight = FontWeight.Medium,
            color = if (value == "--:--:--") TextSecondary else OnSurface)
    }
}

/** 月相/节气一行流式: "标题: ●{名称}{月}/{日} ..." 排不下自动折行 */
@OptIn(androidx.compose.foundation.layout.ExperimentalLayoutApi::class)
@Composable
private fun EventFlowLine(
    title: String, events: List<MonthEvent>, month: Int, dotColor: Color
) {
    androidx.compose.foundation.layout.FlowRow(
        Modifier.fillMaxWidth().padding(top = 4.dp),
        horizontalArrangement = Arrangement.spacedBy(8.dp),
        verticalArrangement = Arrangement.spacedBy(2.dp)
    ) {
        Text("$title: ", fontSize = Dimens.fontSmall, fontWeight = FontWeight.Medium,
            color = Primary)
        events.forEach { e ->
            Row(verticalAlignment = Alignment.CenterVertically) {
                Text("●", fontSize = 8.sp, color = dotColor,
                    modifier = Modifier.padding(end = 2.dp))
                Text("${e.name}${month}/${e.day}",
                    fontSize = Dimens.fontSmall, color = OnSurface)
            }
        }
    }
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
    onDismiss: () -> Unit,
    onOpenAlmanacSheet: () -> Unit
) {
    val isToday = day.solarYear == todayY &&
                  day.solarMonth == todayM &&
                  day.solarDay == todayD

    // 老黄历摘要按需加载 (弹窗内只展示一行速览, 全量走底部抽屉)
    var almanac by remember(day.solarYear, day.solarMonth, day.solarDay) {
        mutableStateOf<DayAlmanac?>(null)
    }
    LaunchedEffect(day.solarYear, day.solarMonth, day.solarDay) {
        almanac = CalendarRepository.getAlmanac(
            day.solarYear, day.solarMonth, day.solarDay)
    }

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
                    .widthIn(max = 320.dp)
                    .padding(horizontal = Dimens.paddingLg),
                shape = RoundedCornerShape(Dimens.radiusMd),
                color = PopupBg,
                tonalElevation = 8.dp,
                shadowElevation = 12.dp,
                border = androidx.compose.foundation.BorderStroke(0.5.dp, PopupBorder)
            ) {
                // 紧凑内容 (不再内嵌全量老黄历, 无需滚动) — 卡片小, 点空白即退出
                Column(
                    Modifier
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
                        "${day.yearGZ}年",
                        day.lunarMonthName,
                        "${day.lunarDayName}日"
                    )
                    PopupInfoLine(
                        "${day.yearGZ}年",
                        "${day.monthGZ}月",
                        "${day.dayGZ}日"
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

                    // ── 老黄历速览 (一行摘要, 对齐鸿蒙 AlmanacInlineRow) ──
                    val a = almanac
                    if (a != null) {
                        Spacer(Modifier.height(6.dp))
                        PopupDividerLine()
                        Spacer(Modifier.height(4.dp))
                        AlmanacInlineRow(a)
                    }

                    // ── 老黄历详情入口 → 打开底部抽屉 (点空白/下拉即退出) ──
                    Spacer(Modifier.height(6.dp))
                    Row(
                        Modifier
                            .fillMaxWidth()
                            .clip(RoundedCornerShape(Dimens.radiusSm))
                            .clickable(onClick = onOpenAlmanacSheet)
                            .padding(vertical = 6.dp),
                        horizontalArrangement = Arrangement.Center,
                        verticalAlignment = Alignment.CenterVertically
                    ) {
                        Text("📜 老黄历详情 ›",
                            fontSize = Dimens.fontCaption,
                            fontWeight = FontWeight.Medium,
                            color = PopupGold)
                    }
                }
            }
        }
    }
}

// ─── 老黄历速览 (popup 内嵌入, 对齐鸿蒙 AlmanacInlineRow) ───────────────
//   第 1 行: 二十八宿 + 十二神(黄道/黑道) + 冲煞
//   第 2 行: 宜 (前 3 项) / 忌 (前 3 项) 色块徽章
@Composable
private fun AlmanacInlineRow(a: DayAlmanac) {
    Column(Modifier.fillMaxWidth()) {
        Row(
            Modifier.fillMaxWidth(),
            horizontalArrangement = Arrangement.Center,
            verticalAlignment = Alignment.CenterVertically
        ) {
            Text("📜", fontSize = 11.sp, color = PopupGold,
                modifier = Modifier.padding(end = 4.dp))
            Text("${a.xiu}${a.xiuZheng}${a.xiuAnimal}",
                fontSize = Dimens.fontCaption, fontWeight = FontWeight.Medium,
                color = if (a.xiuLuck == "吉") PopupGreen else PopupRed)
            Text("·", fontSize = Dimens.fontSmall, color = PopupSub,
                modifier = Modifier.padding(horizontal = 4.dp))
            Text("${a.twelveGod}${a.huangHei}",
                fontSize = Dimens.fontCaption, fontWeight = FontWeight.Medium,
                color = if (a.isHuangDao) PopupGold else PopupSub)
            Text("·", fontSize = Dimens.fontSmall, color = PopupSub,
                modifier = Modifier.padding(horizontal = 4.dp))
            Text("冲${a.chongShengXiao}煞${a.sha}",
                fontSize = Dimens.fontCaption, color = PopupSub)
        }
        if (a.yi.isNotEmpty() || a.ji.isNotEmpty()) {
            Spacer(Modifier.height(3.dp))
            Row(
                Modifier.fillMaxWidth(),
                horizontalArrangement = Arrangement.Center,
                verticalAlignment = Alignment.CenterVertically
            ) {
                if (a.yi.isNotEmpty()) {
                    Text("宜", fontSize = Dimens.fontCaption, fontWeight = FontWeight.Medium,
                        color = PopupText,
                        modifier = Modifier
                            .background(PopupGreen, RoundedCornerShape(4.dp))
                            .padding(horizontal = 5.dp, vertical = 1.dp))
                    Text(" " + a.yi.take(3).joinToString(" "),
                        fontSize = Dimens.fontCaption, color = PopupSub,
                        modifier = Modifier.padding(start = 4.dp, end = 6.dp),
                        maxLines = 1)
                }
                if (a.ji.isNotEmpty()) {
                    Text("忌", fontSize = Dimens.fontCaption, fontWeight = FontWeight.Medium,
                        color = PopupText,
                        modifier = Modifier
                            .background(PopupRed, RoundedCornerShape(4.dp))
                            .padding(horizontal = 5.dp, vertical = 1.dp))
                    Text(" " + a.ji.take(3).joinToString(" "),
                        fontSize = Dimens.fontCaption, color = PopupSub,
                        modifier = Modifier.padding(start = 4.dp),
                        maxLines = 1)
                }
            }
        }
    }
}

// ─── 老黄历底部抽屉 (ModalBottomSheet, 对齐鸿蒙 bindSheet) ───────────────
//   · 点遮罩/下拉即关闭, 不再需要后退键
//   · 内含全量 AlmanacPanel + "董公择日经典知识"入口
@OptIn(ExperimentalMaterial3Api::class)
@Composable
private fun AlmanacSheetBottomSheet(
    day: DayInfo,
    onDismiss: () -> Unit,
    onOpenTopics: () -> Unit
) {
    val sheetState = rememberModalBottomSheetState(skipPartiallyExpanded = true)
    var almanac by remember(day.solarYear, day.solarMonth, day.solarDay) {
        mutableStateOf<DayAlmanac?>(null)
    }
    LaunchedEffect(day.solarYear, day.solarMonth, day.solarDay) {
        almanac = CalendarRepository.getAlmanac(
            day.solarYear, day.solarMonth, day.solarDay)
    }

    ModalBottomSheet(
        onDismissRequest = onDismiss,
        sheetState = sheetState,
        containerColor = PopupBg,
        shape = RoundedCornerShape(topStart = Dimens.radiusLg, topEnd = Dimens.radiusLg),
        dragHandle = { BottomSheetDefaults.DragHandle() }
    ) {
        Column(
            Modifier
                .fillMaxWidth()
                .heightIn(max = 620.dp)
                .verticalScroll(rememberScrollState())
                .padding(horizontal = Dimens.paddingLg, vertical = Dimens.paddingSm)
        ) {
            // 日期标题
            Text(
                "${YearUtil.astroYearToStr(day.solarYear)}年${day.solarMonth}月${day.solarDay}日" +
                " · ${day.lunarMonthName}${day.lunarDayName}",
                fontSize = Dimens.fontSubtitle, fontWeight = FontWeight.Bold,
                color = PopupText,
                modifier = Modifier.fillMaxWidth(),
                textAlign = TextAlign.Center
            )
            Spacer(Modifier.height(8.dp))

            val a = almanac
            if (a != null) {
                AlmanacPanel(a)
            } else {
                Text("数据加载中...",
                    fontSize = Dimens.fontBody, color = PopupSub,
                    modifier = Modifier.fillMaxWidth().padding(vertical = 24.dp),
                    textAlign = TextAlign.Center)
            }

            // 董公择日经典知识入口 (与鸿蒙 topicsEntrySection 对齐)
            Spacer(Modifier.height(8.dp))
            Row(
                Modifier
                    .fillMaxWidth()
                    .clip(RoundedCornerShape(Dimens.radiusSm))
                    .background(PopupBorder.copy(alpha = 0.3f))
                    .clickable(onClick = onOpenTopics)
                    .padding(Dimens.paddingMd),
                verticalAlignment = Alignment.CenterVertically
            ) {
                Text("📖 董公择日 · 总论 · 口诀 · 方位",
                    fontSize = Dimens.fontBody, fontWeight = FontWeight.Medium,
                    color = PopupGold, modifier = Modifier.weight(1f))
                Text("›", fontSize = 20.sp, color = PopupGold)
            }
            Spacer(Modifier.height(Dimens.paddingLg))
        }
    }
}

// ─── 老黄历面板 ────────────────────────────────────────────────────────
//
//  紧凑展示: 二十八宿/黄道/冲煞 + 五吉神 + 彭祖 + 神煞 + 宜忌 + 用事 + 吉时.
//  在 AlmanacSheetBottomSheet 内调用, 滚动容器由外层提供.

private val ZhiNames = listOf("子","丑","寅","卯","辰","巳",
                              "午","未","申","酉","戌","亥")

@Composable
private fun AlmanacPanel(a: DayAlmanac) {
    Text("老黄历", fontSize = Dimens.fontCaption,
        fontWeight = FontWeight.Medium, color = PopupGold,
        modifier = Modifier.fillMaxWidth().padding(bottom = 4.dp),
        textAlign = TextAlign.Center)

    PopupInfoLine(
        "${a.xiu}宿(${a.xiuZheng}${a.xiuAnimal})",
        a.twelveGod, a.huangHei,
        color = if (a.isHuangDao) PopupGreen else PopupRed,
        bold = true
    )
    PopupInfoLine(
        "冲${a.chongShengXiao}",
        a.chongGanZhi, "煞${a.sha}"
    )
    // 五方位 (喜/财/福 + 阳贵/阴贵), 与鸿蒙 AlmanacComponents 对齐
    PopupDirRow5(listOf(
        "喜"   to a.xiShenFang,
        "财"   to a.caiShenFang,
        "福"   to a.fuShenFang,
        "阳贵" to a.yangGuiFang,
        "阴贵" to a.yinGuiFang
    ))

    if (a.pengZuGan.isNotEmpty() || a.pengZuZhi.isNotEmpty()) {
        Spacer(Modifier.height(2.dp))
        PopupCenteredText(a.pengZuGan, color = PopupSub)
        PopupCenteredText(a.pengZuZhi, color = PopupSub)
    }

    // 神煞: 完整展示, 按权重降序 (与鸿蒙一致). 同类型可能 20+ 项, 让弹窗
    //   滚动容器承接溢出, 不再做 take(6) 截断 — 老黄历参考价值在"全部命中".
    if (a.shenSha.isNotEmpty()) {
        Spacer(Modifier.height(4.dp))
        val lucky = a.shenSha.filter { it.isLucky }.sortedByDescending { it.weight }
        val bad   = a.shenSha.filter { !it.isLucky }.sortedByDescending { it.weight }
        if (lucky.isNotEmpty()) ShenShaRow("吉神", lucky, PopupGreen)
        if (bad.isNotEmpty())   ShenShaRow("凶煞", bad,   PopupRed)
    }

    // 宜 / 忌
    if (a.yi.isNotEmpty() || a.ji.isNotEmpty()) {
        Spacer(Modifier.height(4.dp))
        if (a.yi.isNotEmpty()) TextLine("宜", a.yi.toList(), PopupGreen)
        if (a.ji.isNotEmpty()) TextLine("忌", a.ji.toList(), PopupRed)
    }

    // 吉时
    if (a.luckyHours.isNotEmpty()) {
        Spacer(Modifier.height(4.dp))
        Row(Modifier.fillMaxWidth()) {
            Text("吉时", fontSize = Dimens.fontSmall, color = PopupSub,
                modifier = Modifier.padding(end = 6.dp))
            Text(a.luckyHours.joinToString("  ") { lh ->
                val z = if (lh.zhi in 0..11) ZhiNames[lh.zhi] else ""
                "${lh.name}($z)"
            }, fontSize = Dimens.fontSmall, color = PopupGold)
        }
    }

    // 用事择吉 — 完整展示 (典型 4-8 条, 不再截断)
    if (a.events.isNotEmpty()) {
        Spacer(Modifier.height(4.dp))
        Text("用事择吉", fontSize = Dimens.fontSmall, color = PopupSub,
            modifier = Modifier.fillMaxWidth(), textAlign = TextAlign.Center)
        for (e in a.events) EventRow(e)
    }

    // 备注 (节气/月度提示)
    if (a.notes.isNotEmpty()) {
        Spacer(Modifier.height(4.dp))
        for (n in a.notes) {
            PopupCenteredText("· $n", color = PopupSub,
                size = Dimens.fontSmall)
        }
    }

    // 董公择日要诀语录 — 完整展示 (与鸿蒙一致, 多条逐一显示)
    if (a.quotes.isNotEmpty()) {
        Spacer(Modifier.height(4.dp))
        for (q in a.quotes) {
            QuoteCard(q)
            Spacer(Modifier.height(2.dp))
        }
    }
}

/**
 * 老黄历静态知识 Dialog — 按 category 分组列出 (董公总论/择日基础理论/口诀/方位 等).
 *
 * 与鸿蒙 AlmanacTopicsSheet 对齐: 顶部分类导航(默认显示首类), 滚动展示该类
 * 全部条目, 标题 + 正文双层结构. topics 通常 ≤ 64 条, 单次渲染无性能问题.
 */
@Composable
private fun AlmanacTopicsDialog(
    topics: List<AlmanacTopic>,
    onDismiss: () -> Unit
) {
    // 按 category 保序去重 — 显示顺序与 C++ 注册顺序一致, 不做字典序重排.
    val categories: List<String> = remember(topics) {
        val seen = LinkedHashSet<String>()
        topics.forEach { seen.add(it.category.ifEmpty { "未分类" }) }
        seen.toList()
    }
    var selectedCat by remember(categories) { mutableStateOf(categories.firstOrNull() ?: "") }

    Dialog(
        onDismissRequest = onDismiss,
        properties = DialogProperties(usePlatformDefaultWidth = false)
    ) {
        Box(Modifier.fillMaxWidth(), contentAlignment = Alignment.Center) {
            Surface(
                modifier = Modifier
                    .widthIn(max = 340.dp)
                    .padding(horizontal = Dimens.paddingLg),
                shape = RoundedCornerShape(Dimens.radiusMd),
                color = PopupBg,
                tonalElevation = 8.dp,
                shadowElevation = 12.dp,
                border = androidx.compose.foundation.BorderStroke(0.5.dp, PopupBorder)
            ) {
                Column(Modifier.padding(horizontal = 12.dp, vertical = 12.dp)) {
                    Text(
                        "📜 老黄历经典 · 董公择日",
                        fontSize = Dimens.fontBody,
                        fontWeight = FontWeight.Bold,
                        color = PopupGold,
                        modifier = Modifier.fillMaxWidth(),
                        textAlign = TextAlign.Center
                    )
                    Spacer(Modifier.height(8.dp))

                    if (topics.isEmpty()) {
                        Text(
                            "暂无数据",
                            fontSize = Dimens.fontCaption,
                            color = PopupSub,
                            modifier = Modifier.fillMaxWidth().padding(vertical = 16.dp),
                            textAlign = TextAlign.Center
                        )
                    } else {
                        // 横向滚动的分类导航条
                        Row(
                            Modifier
                                .fillMaxWidth()
                                .horizontalScroll(rememberScrollState()),
                            horizontalArrangement = Arrangement.spacedBy(6.dp)
                        ) {
                            for (cat in categories) {
                                val sel = cat == selectedCat
                                Box(
                                    Modifier
                                        .clip(RoundedCornerShape(Dimens.radiusSm))
                                        .background(
                                            if (sel) PopupGold.copy(alpha = 0.18f)
                                            else Color.Transparent
                                        )
                                        .clickable { selectedCat = cat }
                                        .padding(horizontal = 10.dp, vertical = 4.dp)
                                ) {
                                    Text(cat,
                                        fontSize = Dimens.fontSmall,
                                        fontWeight = if (sel) FontWeight.Bold else FontWeight.Normal,
                                        color = if (sel) PopupGold else PopupSub)
                                }
                            }
                        }

                        Spacer(Modifier.height(6.dp))
                        PopupDividerLine()
                        Spacer(Modifier.height(6.dp))

                        // 内容滚动区 (限高 ~480dp, 留出对话框周边留白)
                        Column(
                            Modifier
                                .heightIn(max = 480.dp)
                                .verticalScroll(rememberScrollState())
                        ) {
                            val filtered = topics.filter {
                                (it.category.ifEmpty { "未分类" }) == selectedCat
                            }
                            for ((idx, t) in filtered.withIndex()) {
                                if (idx > 0) Spacer(Modifier.height(8.dp))
                                Text(t.title,
                                    fontSize = Dimens.fontCaption,
                                    fontWeight = FontWeight.Bold,
                                    color = PopupGold)
                                Spacer(Modifier.height(2.dp))
                                Text(t.body,
                                    fontSize = Dimens.fontSmall,
                                    color = PopupText,
                                    lineHeight = 18.sp)
                            }
                        }
                    }

                    Spacer(Modifier.height(10.dp))
                    Box(
                        Modifier
                            .fillMaxWidth()
                            .clip(RoundedCornerShape(Dimens.radiusSm))
                            .clickable(onClick = onDismiss)
                            .padding(vertical = 8.dp),
                        contentAlignment = Alignment.Center
                    ) {
                        Text("关闭",
                            fontSize = Dimens.fontCaption,
                            fontWeight = FontWeight.Medium,
                            color = PopupSub)
                    }
                }
            }
        }
    }
}

/**
 * 紧凑 5 列方位 (喜/财/福/阳贵/阴贵), 等间距均分.
 *
 * label 是中文方位简称, dir 是 "东南"/"西北" 等具体方位串; dir 为空时显示 "—".
 */
@Composable
private fun PopupDirRow5(items: List<Pair<String, String>>) {
    Row(
        Modifier.fillMaxWidth().padding(vertical = 1.dp),
        horizontalArrangement = Arrangement.SpaceEvenly
    ) {
        for ((label, dir) in items) {
            Text(
                "$label${dir.ifEmpty { "—" }}",
                fontSize = Dimens.fontSmall,
                fontWeight = FontWeight.Medium,
                color = if (dir.isEmpty()) PopupSub.copy(alpha = 0.4f) else PopupGold,
                maxLines = 1
            )
        }
    }
}

@Composable
private fun ShenShaRow(label: String, items: List<ShenSha>, color: Color) {
    Row(Modifier.fillMaxWidth().padding(vertical = 1.dp)) {
        Text(label, fontSize = Dimens.fontSmall, color = PopupSub,
            modifier = Modifier.padding(end = 6.dp))
        Text(items.joinToString("  ") { it.name },
            fontSize = Dimens.fontSmall, color = color,
            modifier = Modifier.weight(1f))
    }
}

@Composable
private fun TextLine(label: String, items: List<String>, color: Color) {
    Row(Modifier.fillMaxWidth().padding(vertical = 1.dp)) {
        Text(label, fontSize = Dimens.fontCaption,
            fontWeight = FontWeight.Medium, color = color,
            modifier = Modifier.padding(end = 6.dp))
        Text(items.joinToString("、"),
            fontSize = Dimens.fontCaption, color = PopupText,
            modifier = Modifier.weight(1f))
    }
}

@Composable
private fun EventRow(e: EventAdvice) {
    val marker = if (e.suitable) "✓" else "✗"
    val color = if (e.suitable) PopupGreen else PopupRed
    Row(
        Modifier.fillMaxWidth().padding(vertical = 1.dp),
        verticalAlignment = Alignment.CenterVertically
    ) {
        Text(marker, fontSize = Dimens.fontCaption,
            fontWeight = FontWeight.Bold, color = color,
            modifier = Modifier.padding(end = 4.dp))
        Text(e.event, fontSize = Dimens.fontCaption,
            color = PopupText, modifier = Modifier.padding(end = 4.dp))
        if (e.reason.isNotEmpty()) {
            Text("(${e.reason})", fontSize = Dimens.fontSmall, color = PopupSub)
        }
    }
}

@Composable
private fun QuoteCard(q: AlmanacQuote) {
    val accentColor = when (q.luck) {
        "吉" -> PopupGreen
        "凶" -> PopupRed
        "混" -> PopupGold
        else -> PopupSub
    }
    Column(
        Modifier
            .fillMaxWidth()
            .clip(RoundedCornerShape(Dimens.radiusSm))
            .background(PopupBorder.copy(alpha = 0.3f))
            .padding(8.dp)
    ) {
        Row(verticalAlignment = Alignment.CenterVertically) {
            Text(q.title, fontSize = Dimens.fontSmall,
                fontWeight = FontWeight.Medium, color = accentColor,
                modifier = Modifier.weight(1f))
            if (q.luck.isNotEmpty()) {
                Text(q.luck, fontSize = Dimens.fontSmall,
                    fontWeight = FontWeight.Bold, color = accentColor)
            }
        }
        if (q.source.isNotEmpty()) {
            Text("— ${q.source}", fontSize = 9.sp, color = PopupSub,
                modifier = Modifier.padding(top = 1.dp))
        }
        if (q.body.isNotEmpty()) {
            Text(q.body, fontSize = Dimens.fontSmall, color = PopupText,
                modifier = Modifier.padding(top = 4.dp), maxLines = 6)
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
