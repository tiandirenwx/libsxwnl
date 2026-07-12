package com.sxwnl.calendar.ui.components

import androidx.compose.foundation.background
import androidx.compose.foundation.clickable
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.lazy.LazyColumn
import androidx.compose.foundation.shape.RoundedCornerShape
import androidx.compose.foundation.text.BasicTextField
import androidx.compose.foundation.text.KeyboardActions
import androidx.compose.foundation.text.KeyboardOptions
import androidx.compose.material3.*
import androidx.compose.runtime.*
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.draw.clip
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.graphics.SolidColor
import androidx.compose.ui.text.TextStyle
import androidx.compose.ui.text.font.FontWeight
import androidx.compose.ui.text.input.ImeAction
import androidx.compose.ui.text.input.KeyboardType
import androidx.compose.ui.text.style.TextAlign
import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.sp
import com.sxwnl.calendar.data.CalendarRepository
import com.sxwnl.calendar.data.LunarMonth
import com.sxwnl.calendar.data.ReverseItem
import com.sxwnl.calendar.ui.theme.*
import com.sxwnl.calendar.util.BaziCalc
import com.sxwnl.calendar.util.BaziCalc.BirthInputMode
import com.sxwnl.calendar.util.BaziCalc.DateSelection
import kotlinx.coroutines.launch

// ════════════════════════════════════════════════════════════════
//  BaziDateTimePicker — 与鸿蒙 DateTimePicker.ets 对齐
//  三 tab: 公历 / 农历 / 四柱反推
// ════════════════════════════════════════════════════════════════

@OptIn(ExperimentalMaterial3Api::class)
@Composable
fun BaziDateTimePicker(
    initial: DateSelection,
    onDismiss: () -> Unit,
    onConfirm: (DateSelection) -> Unit
) {
    val sheetState = rememberModalBottomSheetState(skipPartiallyExpanded = true)
    ModalBottomSheet(
        onDismissRequest = onDismiss,
        sheetState = sheetState,
        containerColor = Surface,
        shape = RoundedCornerShape(
            topStart = Dimens.radiusLg, topEnd = Dimens.radiusLg
        )
    ) {
        PickerContent(initial, onDismiss, onConfirm)
    }
}

@Composable
private fun PickerContent(
    initial: DateSelection,
    onDismiss: () -> Unit,
    onConfirm: (DateSelection) -> Unit
) {
    val scope = rememberCoroutineScope()

    var inputMode by remember {
        mutableStateOf(
            if (initial.inputMode == BirthInputMode.LUNAR) BirthInputMode.LUNAR
            else BirthInputMode.SOLAR
        )
    }
    var year by remember { mutableIntStateOf(initial.year) }
    var yearText by remember { mutableStateOf("${initial.year}") }
    var month by remember { mutableIntStateOf(initial.month) }
    var day by remember { mutableIntStateOf(initial.day) }
    var hour by remember { mutableIntStateOf(initial.hour) }
    var minute by remember { mutableIntStateOf(initial.minute) }
    var isLeap by remember { mutableStateOf(initial.isLeap) }
    var isSpec by remember { mutableStateOf(initial.isSpec) }
    var lunarMonths by remember { mutableStateOf(emptyList<LunarMonth>()) }
    var validDays by remember { mutableStateOf(emptyList<Int>()) }

    // ── 四柱反推 ──
    var revIdx by remember { mutableStateOf(intArrayOf(0, 2, 0, 0).copyOf()) }
    var revHourUnknown by remember { mutableStateOf(false) }
    var revStartText by remember { mutableStateOf("1900") }
    var revEndText by remember { mutableStateOf("2100") }
    var revResults by remember { mutableStateOf(emptyList<ReverseItem>()) }
    var revSearched by remember { mutableStateOf(false) }
    var pillarPickerIdx by remember { mutableStateOf<Int?>(null) }

    val tabs = listOf("公历", "农历", "四柱反推")
    val pillarNames = listOf("年柱", "月柱", "日柱", "时柱")

    suspend fun refreshLunarMonths() {
        lunarMonths = CalendarRepository.getLunarMonths(year)
    }

    suspend fun refreshValidDays() {
        val ds: List<Int>
        if (inputMode == BirthInputMode.LUNAR) {
            val n = CalendarRepository.getLunarMonthDays(year, month, isLeap, isSpec)
                .let { if (it == 29 || it == 30) it else 30 }
            ds = (1..n).toList()
        } else {
            val d = CalendarRepository.getSolarMonthValidDays(year, month)
            ds = d.ifEmpty {
                (1..BaziCalc.solarDaysInMonth(year, month)).toList()
            }
        }
        validDays = ds
        if (!ds.contains(day)) {
            var pick = ds.firstOrNull() ?: 1
            for (d in ds) if (d <= day) pick = d
            day = pick
        }
    }

    LaunchedEffect(Unit) {
        if (inputMode == BirthInputMode.LUNAR) refreshLunarMonths()
        refreshValidDays()
    }

    fun applyYearText(v: String) {
        yearText = v
        val y = BaziCalc.parseYear(v)
        if (y != null) {
            year = y
            scope.launch {
                if (inputMode == BirthInputMode.LUNAR) refreshLunarMonths()
                refreshValidDays()
            }
        }
    }

    fun confirm() {
        onConfirm(
            DateSelection(
                inputMode = inputMode,
                year = year, month = month, day = day,
                hour = hour, minute = minute,
                isLeap = isLeap, isSpec = isSpec
            )
        )
        onDismiss()
    }

    fun switchTab(target: BirthInputMode) {
        inputMode = target
        scope.launch {
            if (target == BirthInputMode.LUNAR) refreshLunarMonths()
            if (target != BirthInputMode.REVERSE) refreshValidDays()
        }
    }

    Column(Modifier.fillMaxWidth()) {
        // 顶部条
        Row(
            Modifier.fillMaxWidth()
                .padding(horizontal = Dimens.paddingLg, vertical = Dimens.paddingMd),
            verticalAlignment = Alignment.CenterVertically
        ) {
            Text(
                "取消", fontSize = Dimens.fontBody, color = TextSecondary,
                modifier = Modifier.clickable(onClick = onDismiss)
            )
            Spacer(Modifier.weight(1f))
            Row {
                tabs.forEachIndexed { i, t ->
                    val target = when (i) {
                        0 -> BirthInputMode.SOLAR
                        1 -> BirthInputMode.LUNAR
                        else -> BirthInputMode.REVERSE
                    }
                    val active = inputMode == target
                    Text(
                        t,
                        fontSize = Dimens.fontCaption,
                        fontWeight = if (active) FontWeight.Bold else FontWeight.Normal,
                        color = if (active) OnPrimary else TextSecondary,
                        modifier = Modifier
                            .padding(horizontal = 4.dp)
                            .clip(RoundedCornerShape(Dimens.radiusSm))
                            .background(if (active) Accent else Color.Transparent)
                            .clickable { switchTab(target) }
                            .padding(horizontal = 12.dp, vertical = 6.dp)
                    )
                }
            }
            Spacer(Modifier.weight(1f))
            if (inputMode != BirthInputMode.REVERSE) {
                Text(
                    "确认",
                    fontSize = Dimens.fontBody,
                    fontWeight = FontWeight.Bold,
                    color = Accent,
                    modifier = Modifier.clickable { confirm() }
                )
            } else {
                Box(Modifier.width(36.dp))
            }
        }
        HorizontalDivider(color = DividerColor)

        if (inputMode == BirthInputMode.REVERSE) {
            ReverseBody(
                revIdx = revIdx, revHourUnknown = revHourUnknown,
                revStartText = revStartText, revEndText = revEndText,
                revResults = revResults, revSearched = revSearched,
                pillarNames = pillarNames,
                onPickPillar = { pillarPickerIdx = it },
                onHourUnknownChange = { revHourUnknown = it },
                onStartChange = { revStartText = it },
                onEndChange = { revEndText = it },
                onSearch = {
                    val sy = BaziCalc.parseYear(revStartText)
                    val ey = BaziCalc.parseYear(revEndText)
                    if (sy != null && ey != null && sy <= ey) {
                        val g = IntArray(4); val z = IntArray(4)
                        for (i in 0..3) {
                            val idx = ((revIdx[i] % 60) + 60) % 60
                            g[i] = idx % 10
                            z[i] = idx % 12
                        }
                        val hg = if (revHourUnknown) -1 else g[3]
                        val hz = if (revHourUnknown) -1 else z[3]
                        scope.launch {
                            revResults = CalendarRepository.baziReverse(
                                g[0], z[0], g[1], z[1], g[2], z[2],
                                hg, hz, sy, ey
                            )
                            revSearched = true
                        }
                    } else {
                        revResults = emptyList()
                        revSearched = true
                    }
                },
                onApply = { item ->
                    onConfirm(
                        DateSelection(
                            inputMode = BirthInputMode.SOLAR,
                            year = item.year, month = item.month, day = item.day,
                            hour = if (item.hour >= 0) item.hour else 12,
                            minute = 0,
                            isLeap = false, isSpec = false
                        )
                    )
                    onDismiss()
                }
            )
        } else {
            DateWheels(
                inputMode = inputMode,
                yearText = yearText, onYearChange = ::applyYearText,
                year = year,
                lunarMonths = lunarMonths,
                month = month, isLeap = isLeap,
                onMonthChange = { newMonth, newLeap, newSpec ->
                    month = newMonth; isLeap = newLeap; isSpec = newSpec
                    scope.launch { refreshValidDays() }
                },
                validDays = validDays, day = day, onDayChange = { day = it },
                hour = hour, onHourChange = { hour = it },
                minute = minute, onMinuteChange = { minute = it }
            )
        }
    }

    // 反推柱选择 sheet
    val pickerIdx = pillarPickerIdx
    if (pickerIdx != null) {
        val list: List<Int> = when (pickerIdx) {
            1 -> BaziCalc.validMonths(revIdx[0]).toList()
            3 -> BaziCalc.validHours(revIdx[2]).toList()
            else -> (0 until 60).toList()
        }
        PillarPickerSheet(
            title = pillarNames[pickerIdx],
            list = list,
            currentIdx = ((revIdx[pickerIdx] % 60) + 60) % 60,
            onDismiss = { pillarPickerIdx = null },
            onSelected = { picked ->
                val arr = revIdx.copyOf()
                arr[pickerIdx] = picked
                if (pickerIdx == 0) {
                    arr[1] = BaziCalc.remapByZhi(BaziCalc.validMonths(arr[0]), arr[1])
                }
                if (pickerIdx == 2) {
                    arr[3] = BaziCalc.remapByZhi(BaziCalc.validHours(arr[2]), arr[3])
                }
                revIdx = arr
                pillarPickerIdx = null
            }
        )
    }
}

// ─── 公历/农历: 年输入 + 月/日/时/分 ─────────────────────────

@Composable
private fun DateWheels(
    inputMode: BirthInputMode,
    yearText: String, onYearChange: (String) -> Unit, year: Int,
    lunarMonths: List<LunarMonth>,
    month: Int, isLeap: Boolean,
    onMonthChange: (Int, Boolean, Boolean) -> Unit,
    validDays: List<Int>, day: Int, onDayChange: (Int) -> Unit,
    hour: Int, onHourChange: (Int) -> Unit,
    minute: Int, onMinuteChange: (Int) -> Unit
) {
    Column(Modifier.fillMaxWidth().padding(bottom = Dimens.paddingLg)) {
        // 年份输入
        Row(
            Modifier.fillMaxWidth()
                .padding(horizontal = Dimens.paddingLg)
                .padding(top = Dimens.paddingMd),
            verticalAlignment = Alignment.CenterVertically
        ) {
            Text("年份", fontSize = Dimens.fontCaption, color = TextSecondary)
            Box(
                Modifier
                    .padding(start = 8.dp)
                    .weight(1f)
                    .height(38.dp)
                    .clip(RoundedCornerShape(Dimens.radiusSm))
                    .background(Background),
                contentAlignment = Alignment.CenterStart
            ) {
                if (yearText.isEmpty()) {
                    Text(
                        "如 1990 / B212 / -211",
                        color = TextSecondary,
                        fontSize = Dimens.fontBody,
                        modifier = Modifier.padding(horizontal = 10.dp)
                    )
                }
                BasicTextField(
                    value = yearText,
                    onValueChange = onYearChange,
                    textStyle = TextStyle(
                        color = OnSurface, fontSize = Dimens.fontBody
                    ),
                    singleLine = true,
                    modifier = Modifier.fillMaxSize().padding(horizontal = 10.dp),
                    keyboardOptions = KeyboardOptions(imeAction = ImeAction.Done),
                    keyboardActions = KeyboardActions(),
                    cursorBrush = SolidColor(OnSurface)
                )
            }
            Text(
                BaziCalc.formatYear(year),
                fontSize = Dimens.fontCaption, color = Primary,
                modifier = Modifier.padding(start = 8.dp)
            )
        }

        // 滚轮区
        Row(
            Modifier
                .fillMaxWidth()
                .height(180.dp)
                .padding(horizontal = Dimens.paddingSm, vertical = Dimens.paddingSm)
        ) {
            // 月
            if (inputMode == BirthInputMode.LUNAR) {
                val names = lunarMonths.map { it.name }
                val curIdx = lunarMonths.indexOfFirst {
                    it.month == month && it.isLeap == isLeap
                }.coerceAtLeast(0)
                WheelPicker(
                    items = names,
                    selectedIndex = curIdx,
                    onSelectedChange = { idx ->
                        val m = lunarMonths.getOrNull(idx) ?: return@WheelPicker
                        onMonthChange(m.month, m.isLeap, m.isSpec)
                    },
                    modifier = Modifier.weight(1f)
                )
            } else {
                WheelPicker(
                    items = (1..12).map { "${it}月" },
                    selectedIndex = (month - 1).coerceIn(0, 11),
                    onSelectedChange = { idx -> onMonthChange(idx + 1, false, false) },
                    modifier = Modifier.weight(1f)
                )
            }

            // 日 — 农历模式显示"初一/廿九"式名称, 公历显示"1日"
            WheelPicker(
                items = validDays.map {
                    if (inputMode == BirthInputMode.LUNAR) BaziCalc.lunarDayName(it)
                    else "${it}日"
                },
                selectedIndex = validDays.indexOf(day).coerceAtLeast(0),
                onSelectedChange = { idx ->
                    validDays.getOrNull(idx)?.let(onDayChange)
                },
                modifier = Modifier.weight(1f)
            )

            // 时
            WheelPicker(
                items = (0..23).map { "${it}时" },
                selectedIndex = hour.coerceIn(0, 23),
                onSelectedChange = onHourChange,
                modifier = Modifier.weight(1f)
            )

            // 分
            WheelPicker(
                items = (0..59).map { "${it}分" },
                selectedIndex = minute.coerceIn(0, 59),
                onSelectedChange = onMinuteChange,
                modifier = Modifier.weight(1f)
            )
        }
    }
}

// ─── 反推 ──────────────────────────────────────────────────

@Composable
private fun ReverseBody(
    revIdx: IntArray,
    revHourUnknown: Boolean,
    revStartText: String,
    revEndText: String,
    revResults: List<ReverseItem>,
    revSearched: Boolean,
    pillarNames: List<String>,
    onPickPillar: (Int) -> Unit,
    onHourUnknownChange: (Boolean) -> Unit,
    onStartChange: (String) -> Unit,
    onEndChange: (String) -> Unit,
    onSearch: () -> Unit,
    onApply: (ReverseItem) -> Unit
) {
    Column(
        Modifier
            .fillMaxWidth()
            .padding(horizontal = Dimens.paddingLg, vertical = Dimens.paddingMd)
    ) {
        Row(Modifier.fillMaxWidth()) {
            pillarNames.forEachIndexed { i, name ->
                Column(
                    Modifier
                        .weight(1f)
                        .padding(end = if (i < 3) 6.dp else 0.dp)
                        .clip(RoundedCornerShape(Dimens.radiusSm))
                        .background(Background)
                        .clickable(enabled = !(i == 3 && revHourUnknown)) {
                            onPickPillar(i)
                        }
                        .height(70.dp)
                        .padding(vertical = 4.dp),
                    horizontalAlignment = Alignment.CenterHorizontally,
                    verticalArrangement = Arrangement.Center
                ) {
                    Text(name, fontSize = Dimens.fontSmall, color = TextSecondary)
                    Text(
                        if (revHourUnknown && i == 3) "未知"
                        else BaziCalc.jiaZiName(revIdx[i]),
                        fontSize = 20.sp,
                        fontWeight = FontWeight.Bold,
                        color = if (revHourUnknown && i == 3) TextSecondary else Primary,
                        modifier = Modifier.padding(top = 4.dp)
                    )
                    Text("▾", fontSize = Dimens.fontSmall, color = TextSecondary)
                }
            }
        }

        Row(
            Modifier.fillMaxWidth().padding(top = Dimens.paddingMd),
            verticalAlignment = Alignment.CenterVertically
        ) {
            Switch(
                checked = revHourUnknown,
                onCheckedChange = onHourUnknownChange,
                colors = SwitchDefaults.colors(checkedTrackColor = Accent)
            )
            Text(
                "时辰未知",
                fontSize = Dimens.fontCaption, color = TextSecondary,
                modifier = Modifier.padding(start = 4.dp)
            )
        }

        Row(
            Modifier.fillMaxWidth().padding(top = Dimens.paddingMd),
            verticalAlignment = Alignment.CenterVertically
        ) {
            Text("搜索范围", fontSize = Dimens.fontCaption, color = TextSecondary)
            CompactTextField(
                value = revStartText, onValueChange = onStartChange,
                modifier = Modifier.padding(start = 6.dp).width(80.dp).height(36.dp)
            )
            Text(
                "—", fontSize = Dimens.fontCaption, color = TextSecondary,
                modifier = Modifier.padding(horizontal = 6.dp)
            )
            CompactTextField(
                value = revEndText, onValueChange = onEndChange,
                modifier = Modifier.width(80.dp).height(36.dp)
            )
            Spacer(Modifier.weight(1f))
            Button(
                onClick = onSearch,
                colors = ButtonDefaults.buttonColors(
                    containerColor = Accent, contentColor = PrimaryDark
                ),
                shape = RoundedCornerShape(Dimens.radiusSm),
                modifier = Modifier.height(36.dp),
                contentPadding = PaddingValues(horizontal = 14.dp)
            ) { Text("反推", fontSize = Dimens.fontCaption) }
        }

        if (revSearched) {
            Text(
                "匹配结果（${revResults.size} 条" +
                    (if (revResults.size >= 60) "，已截断" else "") +
                    "），点击选用",
                fontSize = Dimens.fontSmall, color = TextSecondary,
                modifier = Modifier.padding(top = Dimens.paddingSm)
            )
            LazyColumn(
                Modifier.fillMaxWidth().height(160.dp)
            ) {
                items(revResults.size) { i ->
                    val it = revResults[i]
                    Row(
                        Modifier
                            .fillMaxWidth()
                            .clickable { onApply(it) }
                            .padding(vertical = 8.dp),
                        verticalAlignment = Alignment.CenterVertically
                    ) {
                        Text(
                            buildString {
                                append(BaziCalc.formatYear(it.year))
                                append("${it.month}月${it.day}日")
                                if (it.hour >= 0) append(" ${it.hour}时")
                            },
                            fontSize = Dimens.fontCaption, color = OnSurface,
                            modifier = Modifier.weight(1f)
                        )
                        Text(
                            "选用 ›",
                            fontSize = Dimens.fontCaption, color = Accent
                        )
                    }
                    HorizontalDivider(color = DividerColor, thickness = 0.5.dp)
                }
            }
        }
    }
}

@Composable
private fun CompactTextField(
    value: String, onValueChange: (String) -> Unit,
    modifier: Modifier = Modifier
) {
    Box(
        modifier
            .clip(RoundedCornerShape(Dimens.radiusSm))
            .background(Background)
            .padding(horizontal = 6.dp),
        contentAlignment = Alignment.CenterStart
    ) {
        BasicTextField(
            value = value, onValueChange = onValueChange,
            textStyle = TextStyle(color = OnSurface, fontSize = Dimens.fontCaption),
            singleLine = true,
            modifier = Modifier.fillMaxWidth(),
            keyboardOptions = KeyboardOptions(keyboardType = KeyboardType.Text),
            cursorBrush = SolidColor(OnSurface)
        )
    }
}

@OptIn(ExperimentalMaterial3Api::class)
@Composable
private fun PillarPickerSheet(
    title: String,
    list: List<Int>,
    currentIdx: Int,
    onDismiss: () -> Unit,
    onSelected: (Int) -> Unit
) {
    val sheetState = rememberModalBottomSheetState()
    var pick by remember(currentIdx) {
        mutableIntStateOf(list.indexOf(currentIdx).coerceAtLeast(0))
    }
    ModalBottomSheet(
        onDismissRequest = onDismiss,
        sheetState = sheetState,
        containerColor = Surface
    ) {
        Column {
            Row(
                Modifier
                    .fillMaxWidth()
                    .padding(Dimens.paddingLg),
                verticalAlignment = Alignment.CenterVertically
            ) {
                Text(
                    "取消", color = TextSecondary,
                    modifier = Modifier.clickable(onClick = onDismiss)
                )
                Spacer(Modifier.weight(1f))
                Text(title, fontSize = Dimens.fontBody,
                    fontWeight = FontWeight.Bold, color = OnSurface)
                Spacer(Modifier.weight(1f))
                Text(
                    "确定", color = Accent,
                    modifier = Modifier.clickable {
                        list.getOrNull(pick)?.let(onSelected)
                    }
                )
            }
            HorizontalDivider(color = DividerColor)
            WheelPicker(
                items = list.map { BaziCalc.jiaZiName(it) },
                selectedIndex = pick,
                onSelectedChange = { pick = it },
                modifier = Modifier.fillMaxWidth()
            )
        }
    }
}

