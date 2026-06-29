package com.sxwnl.calendar.ui.screens

import androidx.compose.foundation.background
import androidx.compose.foundation.border
import androidx.compose.foundation.clickable
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.rememberScrollState
import androidx.compose.foundation.shape.CircleShape
import androidx.compose.foundation.shape.RoundedCornerShape
import androidx.compose.foundation.text.BasicTextField
import androidx.compose.foundation.text.KeyboardActions
import androidx.compose.foundation.text.KeyboardOptions
import androidx.compose.foundation.verticalScroll
import androidx.compose.material3.*
import androidx.compose.runtime.*
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.draw.clip
import androidx.compose.ui.graphics.Brush
import androidx.compose.ui.graphics.SolidColor
import androidx.compose.ui.text.TextStyle
import androidx.compose.ui.text.font.FontWeight
import androidx.compose.ui.text.input.ImeAction
import androidx.compose.ui.text.style.TextAlign
import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.sp
import com.sxwnl.calendar.data.CalendarRepository
import com.sxwnl.calendar.data.YearCalJieQi
import com.sxwnl.calendar.data.YearCalMonth
import com.sxwnl.calendar.ui.theme.*
import com.sxwnl.calendar.util.YearUtil
import kotlinx.coroutines.launch
import java.util.Calendar

// ════════════════════════════════════════════════════════════════
//  YearCalendarScreen — 与鸿蒙端 YearCalendarPage.ets 对齐的年历页面
// ════════════════════════════════════════════════════════════════

@Composable
fun YearCalendarScreen() {
    val scope = rememberCoroutineScope()
    val todayYear = remember { Calendar.getInstance().get(Calendar.YEAR) }
    var currentYear by remember { mutableIntStateOf(todayYear) }
    var yearInput by remember { mutableStateOf(YearUtil.astroYearToStr(todayYear)) }
    var months by remember { mutableStateOf(emptyList<YearCalMonth>()) }
    var loading by remember { mutableStateOf(false) }

    fun reload() {
        loading = true
        scope.launch {
            val data = CalendarRepository.getYearCalendar(currentYear)
            months = data
            yearInput = YearUtil.astroYearToStr(currentYear)
            loading = false
        }
    }

    fun shiftYear(delta: Int) {
        currentYear += delta
        reload()
    }

    fun goCurYear() {
        currentYear = todayYear
        reload()
    }

    fun applyYearInput() {
        val y = YearUtil.yearStrToAstro(yearInput)
        if (YearUtil.isAstroYearValid(y)) {
            currentYear = y
            reload()
        } else {
            yearInput = YearUtil.astroYearToStr(currentYear)
        }
    }

    LaunchedEffect(Unit) { reload() }

    Column(Modifier.fillMaxSize().background(Background)) {
        YearHeaderSection(currentYear)
        YearNavSection(
            yearInput = yearInput, onYearChange = { yearInput = it },
            onMinus10 = { shiftYear(-10) },
            onMinus1 = { shiftYear(-1) },
            onPlus1 = { shiftYear(1) },
            onPlus10 = { shiftYear(10) },
            onApply = ::applyYearInput,
            onToday = ::goCurYear
        )
        SummaryBar(months)

        if (loading) {
            Column(
                Modifier.fillMaxSize(),
                verticalArrangement = Arrangement.Center,
                horizontalAlignment = Alignment.CenterHorizontally
            ) {
                CircularProgressIndicator(color = Primary, modifier = Modifier.size(40.dp))
                Text(
                    "正在计算 ${YearUtil.astroYearToStr(currentYear)} 年历...",
                    fontSize = Dimens.fontCaption, color = TextSecondary,
                    modifier = Modifier.padding(top = 8.dp)
                )
            }
        } else {
            Column(
                Modifier
                    .fillMaxSize()
                    .verticalScroll(rememberScrollState())
                    .padding(horizontal = Dimens.paddingMd)
                    .padding(top = Dimens.paddingMd),
                verticalArrangement = Arrangement.spacedBy(Dimens.paddingSm)
            ) {
                months.forEach { m ->
                    MonthRow(m)
                }
                Text(
                    "数据来源: 寿星万年历内核(C++) · 节气精确到秒",
                    fontSize = Dimens.fontSmall, color = TextSecondary,
                    textAlign = TextAlign.Center,
                    modifier = Modifier
                        .fillMaxWidth()
                        .padding(vertical = Dimens.paddingLg)
                )
            }
        }
    }
}

// ─── Header ─────────────────────────────────────────────────────────────

@Composable
private fun YearHeaderSection(year: Int) {
    Box(
        Modifier
            .fillMaxWidth()
            .background(Brush.horizontalGradient(listOf(GradientStart, GradientEnd)))
            .padding(horizontal = Dimens.paddingLg, vertical = Dimens.paddingSm)
            .padding(top = Dimens.paddingMd)
    ) {
        Column(verticalArrangement = Arrangement.spacedBy(2.dp)) {
            Text(
                "${YearUtil.astroYearToStr(year)} 年",
                fontSize = Dimens.fontTitle, fontWeight = FontWeight.Bold,
                color = OnPrimary
            )
            Text(
                "农历年视图 (按朔月排列)",
                fontSize = Dimens.fontSmall,
                color = OnPrimary.copy(alpha = 0.7f)
            )
        }
    }
}

// ─── Nav ────────────────────────────────────────────────────────────────

@Composable
private fun YearNavSection(
    yearInput: String, onYearChange: (String) -> Unit,
    onMinus10: () -> Unit, onMinus1: () -> Unit,
    onPlus1: () -> Unit, onPlus10: () -> Unit,
    onApply: () -> Unit, onToday: () -> Unit
) {
    Row(
        Modifier
            .fillMaxWidth()
            .background(Brush.horizontalGradient(listOf(GradientStart, GradientEnd)))
            .padding(horizontal = Dimens.paddingLg, vertical = Dimens.paddingSm),
        verticalAlignment = Alignment.CenterVertically
    ) {
        YearNavIconButton("«", onMinus10)
        YearNavIconButton("‹", onMinus1)
        Spacer(Modifier.width(6.dp))

        Box(
            Modifier
                .width(80.dp).height(32.dp)
                .clip(RoundedCornerShape(Dimens.radiusSm))
                .background(PrimaryLight.copy(alpha = 0.4f)),
            contentAlignment = Alignment.Center
        ) {
            if (yearInput.isEmpty()) {
                Text(
                    "YYYY/B212", color = OnPrimary.copy(alpha = 0.4f),
                    fontSize = Dimens.fontBody
                )
            }
            BasicTextField(
                value = yearInput, onValueChange = onYearChange,
                textStyle = TextStyle(
                    color = OnPrimary, fontSize = Dimens.fontBody,
                    textAlign = TextAlign.Center
                ),
                singleLine = true,
                modifier = Modifier.fillMaxSize().padding(horizontal = 4.dp),
                keyboardOptions = KeyboardOptions(imeAction = ImeAction.Done),
                keyboardActions = KeyboardActions(onDone = { onApply() }),
                cursorBrush = SolidColor(OnPrimary)
            )
        }

        Text(
            "年", fontSize = Dimens.fontCaption, color = OnPrimary,
            modifier = Modifier.padding(start = 2.dp, end = 6.dp)
        )

        YearNavIconButton("›", onPlus1)
        YearNavIconButton("»", onPlus10)

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
            Text("本年", fontSize = Dimens.fontCaption, fontWeight = FontWeight.Medium)
        }
    }
}

@Composable
private fun YearNavIconButton(icon: String, onClick: () -> Unit) {
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

// ─── Summary ────────────────────────────────────────────────────────────

@Composable
private fun SummaryBar(months: List<YearCalMonth>) {
    val hasLeap = months.any { it.isLeap }
    val leapName = months.firstOrNull { it.isLeap }?.monthName ?: ""
    val totalJq = months.sumOf { it.jieQi.size }

    Row(
        Modifier
            .fillMaxWidth()
            .background(Surface)
            .padding(horizontal = Dimens.paddingLg, vertical = Dimens.paddingSm)
            .border(0.5.dp, DividerColor),
        verticalAlignment = Alignment.CenterVertically
    ) {
        Text(
            "本年共 ${months.size} 个农历月",
            fontSize = Dimens.fontCaption, color = TextSecondary
        )
        if (hasLeap) {
            Text(
                "  · 含 $leapName",
                fontSize = Dimens.fontCaption, fontWeight = FontWeight.Medium,
                color = JieQiColor
            )
        }
        Spacer(Modifier.weight(1f))
        Text(
            "24 节气 / $totalJq 项",
            fontSize = Dimens.fontSmall, color = TextSecondary
        )
    }
}

// ─── Month Row ──────────────────────────────────────────────────────────

@Composable
private fun MonthRow(m: YearCalMonth) {
    Column(
        Modifier
            .fillMaxWidth()
            .clip(RoundedCornerShape(Dimens.radiusMd))
            .background(Surface)
            .padding(horizontal = Dimens.paddingLg, vertical = Dimens.paddingSm)
    ) {
        Row(
            verticalAlignment = Alignment.CenterVertically,
            modifier = Modifier
                .fillMaxWidth()
                .padding(top = Dimens.paddingMd, bottom = Dimens.paddingSm)
        ) {
            Text(
                m.monthName,
                fontSize = Dimens.fontSubtitle, fontWeight = FontWeight.Bold,
                color = if (m.isLeap) JieQiColor else Primary,
                modifier = Modifier.width(70.dp)
            )

            Column(
                Modifier.padding(start = Dimens.paddingMd),
                verticalArrangement = Arrangement.spacedBy(0.dp)
            ) {
                Text(
                    "朔 ${m.shuoGZ}",
                    fontSize = Dimens.fontBody, fontWeight = FontWeight.Medium,
                    color = OnSurface
                )
                Text(
                    formatShuoDate(m),
                    fontSize = Dimens.fontSmall, color = TextSecondary
                )
            }

            Spacer(Modifier.weight(1f))

            Box(
                Modifier
                    .clip(RoundedCornerShape(4.dp))
                    .background(if (m.dayCount > 29) Primary else TextSecondary)
                    .padding(horizontal = 6.dp, vertical = 2.dp)
                    .padding(end = 0.dp)
            ) {
                Text(
                    if (m.dayCount > 29) "大" else "小",
                    fontSize = Dimens.fontCaption, color = OnPrimary
                )
            }

            Text(
                "${m.dayCount} 天",
                fontSize = Dimens.fontSmall, color = TextSecondary,
                modifier = Modifier.padding(start = 4.dp)
            )
        }

        if (m.jieQi.isNotEmpty()) {
            Column(
                Modifier
                    .fillMaxWidth()
                    .clip(RoundedCornerShape(Dimens.radiusSm))
                    .background(Background)
                    .padding(
                        horizontal = Dimens.paddingMd,
                        vertical = Dimens.paddingSm
                    )
                    .padding(top = Dimens.paddingXs)
            ) {
                m.jieQi.forEach { jq ->
                    JieQiRow(jq, m)
                }
            }
        }
    }
}

@Composable
private fun JieQiRow(jq: YearCalJieQi, m: YearCalMonth) {
    Row(
        Modifier
            .fillMaxWidth()
            .padding(vertical = 2.dp),
        verticalAlignment = Alignment.CenterVertically
    ) {
        Text(
            jq.dayName,
            fontSize = Dimens.fontCaption, color = TextSecondary,
            modifier = Modifier.width(36.dp)
        )
        Text(
            jq.gz,
            fontSize = Dimens.fontCaption, fontWeight = FontWeight.Medium,
            color = OnSurface,
            modifier = Modifier.width(40.dp)
        )
        Text(
            jq.name,
            fontSize = Dimens.fontBody, fontWeight = FontWeight.Medium,
            color = JieQiColor,
            modifier = Modifier.width(48.dp)
        )
        Text(
            formatJieQiDate(jq, m),
            fontSize = Dimens.fontCaption, color = TextSecondary,
            modifier = Modifier.weight(1f)
        )
    }
}

// ─── Helpers ────────────────────────────────────────────────────────────

private fun formatShuoDate(m: YearCalMonth): String =
    "${YearUtil.astroYearToStr(m.solarYear)}-${pad2(m.solarMonth)}-${pad2(m.solarDay)}"

private fun formatJieQiDate(jq: YearCalJieQi, m: YearCalMonth): String {
    val date = "${pad2(jq.solarMonth)}-${pad2(jq.solarDay)}"
    return if (jq.time.isNotEmpty()) "$date ${jq.time}"
    else "${YearUtil.astroYearToStr(m.solarYear)}-$date"
}

private fun pad2(n: Int): String = if (n < 10) "0$n" else "$n"
