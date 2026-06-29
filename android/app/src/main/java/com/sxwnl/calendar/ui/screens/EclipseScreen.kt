package com.sxwnl.calendar.ui.screens

import androidx.compose.foundation.BorderStroke
import androidx.compose.foundation.background
import androidx.compose.foundation.border
import androidx.compose.foundation.clickable
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.lazy.LazyColumn
import androidx.compose.foundation.lazy.items
import androidx.compose.foundation.rememberScrollState
import androidx.compose.foundation.shape.RoundedCornerShape
import androidx.compose.foundation.text.BasicTextField
import androidx.compose.foundation.text.KeyboardOptions
import androidx.compose.foundation.verticalScroll
import androidx.compose.material3.*
import androidx.compose.runtime.*
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.draw.clip
import androidx.compose.ui.draw.shadow
import androidx.compose.ui.graphics.Brush
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.graphics.SolidColor
import androidx.compose.ui.text.TextStyle
import androidx.compose.ui.text.font.FontFamily
import androidx.compose.ui.text.font.FontWeight
import androidx.compose.ui.text.input.KeyboardType
import androidx.compose.ui.text.style.TextAlign
import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.sp
import com.sxwnl.calendar.data.*
import com.sxwnl.calendar.ui.components.*
import com.sxwnl.calendar.ui.theme.*
import com.sxwnl.calendar.util.EclipseUtil
import com.sxwnl.calendar.util.EclipseShareUtil
import androidx.compose.material.icons.Icons
import androidx.compose.material.icons.filled.Share
import androidx.compose.material.icons.filled.Event
import androidx.compose.ui.platform.LocalContext
import androidx.compose.ui.platform.LocalView
import kotlinx.coroutines.launch
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.withContext
import java.util.Calendar

// ════════════════════════════════════════════════════════════════
//  EclipseScreen — 优雅版日月食工具
//
//  ① 顶部 banner + 搜索面板 (年/月起/数量, 日食/月食 segmented)
//  ② 卡片列表 — 加 emoji + 渐变徽章
//  ③ 点击 → ModalBottomSheet, 三 tab:
//      日食: 食带地图 | 本地观测 | 数据
//      月食: 过程动画 | 数据
// ════════════════════════════════════════════════════════════════

@OptIn(ExperimentalMaterial3Api::class)
@Composable
fun EclipseScreen() {
    val scope = rememberCoroutineScope()
    val now = remember { Calendar.getInstance() }

    var isSolar by remember { mutableStateOf(true) }
    var yearStr by remember { mutableStateOf("${now.get(Calendar.YEAR)}") }
    var monthStr by remember { mutableStateOf("1") }
    var countStr by remember { mutableStateOf("8") }

    var solarItems by remember { mutableStateOf<List<SolarEclipseItem>>(emptyList()) }
    var lunarItems by remember { mutableStateOf<List<LunarEclipseItem>>(emptyList()) }
    var loading by remember { mutableStateOf(false) }
    var searched by remember { mutableStateOf(false) }

    var selectedSolar by remember { mutableStateOf<SolarEclipseItem?>(null) }
    var selectedLunar by remember { mutableStateOf<LunarEclipseItem?>(null) }

    fun search() {
        loading = true
        searched = true
        scope.launch {
            val y = yearStr.toIntOrNull() ?: now.get(Calendar.YEAR)
            val m = monthStr.toIntOrNull()?.coerceIn(1, 12) ?: 1
            val c = countStr.toIntOrNull()?.coerceIn(1, 30) ?: 8
            if (isSolar) {
                solarItems = CalendarRepository.searchSolarEclipses(y, m, c)
                lunarItems = emptyList()
            } else {
                lunarItems = CalendarRepository.searchLunarEclipses(y, m, c)
                solarItems = emptyList()
            }
            loading = false
        }
    }

    Column(Modifier.fillMaxSize().background(Background)) {
        HeaderBanner()
        SearchPanel(
            isSolar = isSolar, onSolarChange = { isSolar = it },
            yearStr = yearStr, onYearChange = { yearStr = it },
            monthStr = monthStr, onMonthChange = { monthStr = it },
            countStr = countStr, onCountChange = { countStr = it },
            loading = loading, onSearch = ::search
        )

        if (searched) {
            if (isSolar) SolarList(solarItems) { selectedSolar = it }
            else LunarList(lunarItems) { selectedLunar = it }
        } else {
            EmptyHint("点击「开始搜索」查询日月食")
        }
    }

    selectedSolar?.let { item ->
        ModalBottomSheet(
            onDismissRequest = { selectedSolar = null },
            containerColor = Surface,
            dragHandle = { BottomSheetDefaults.DragHandle() }
        ) {
            SolarDetailSheet(item)
        }
    }
    selectedLunar?.let { item ->
        ModalBottomSheet(
            onDismissRequest = { selectedLunar = null },
            containerColor = Surface,
            dragHandle = { BottomSheetDefaults.DragHandle() }
        ) {
            LunarDetailSheet(item)
        }
    }
}

// ─── 顶部 banner ────────────────────────────────────────────

@Composable
private fun HeaderBanner() {
    Column(
        Modifier.fillMaxWidth()
            .background(Brush.linearGradient(listOf(GradientStart, GradientEnd)))
            .padding(top = Dimens.paddingXl, bottom = Dimens.paddingLg),
        horizontalAlignment = Alignment.CenterHorizontally
    ) {
        Text("\uD83C\uDF11 日月食工具",
            fontSize = Dimens.fontTitle, fontWeight = FontWeight.Bold, color = OnPrimary)
        Spacer(Modifier.height(4.dp))
        Text("基于现代天文算法 · 精确到秒级",
            fontSize = Dimens.fontCaption, color = OnPrimary.copy(alpha = 0.85f))
    }
}

// ─── 搜索面板 ──────────────────────────────────────────────

@Composable
private fun SearchPanel(
    isSolar: Boolean, onSolarChange: (Boolean) -> Unit,
    yearStr: String, onYearChange: (String) -> Unit,
    monthStr: String, onMonthChange: (String) -> Unit,
    countStr: String, onCountChange: (String) -> Unit,
    loading: Boolean, onSearch: () -> Unit
) {
    Column(
        Modifier.padding(Dimens.paddingMd).fillMaxWidth()
            .shadow(3.dp, RoundedCornerShape(Dimens.radiusLg))
            .clip(RoundedCornerShape(Dimens.radiusLg))
            .background(Surface)
            .padding(Dimens.paddingMd),
        verticalArrangement = Arrangement.spacedBy(10.dp)
    ) {
        Row(Modifier.fillMaxWidth().clip(RoundedCornerShape(20.dp))
            .background(Background)) {
            listOf(true to "日食 ☀️", false to "月食 🌕").forEach { (v, label) ->
                val sel = isSolar == v
                Box(
                    Modifier.weight(1f).height(36.dp)
                        .clip(RoundedCornerShape(20.dp))
                        .background(if (sel) Primary else Color.Transparent)
                        .clickable { onSolarChange(v) },
                    contentAlignment = Alignment.Center
                ) {
                    Text(label, color = if (sel) OnPrimary else TextSecondary,
                        fontSize = Dimens.fontBody,
                        fontWeight = if (sel) FontWeight.Bold else FontWeight.Normal)
                }
            }
        }

        Row(verticalAlignment = Alignment.CenterVertically,
            horizontalArrangement = Arrangement.spacedBy(6.dp)) {
            Text("起始", fontSize = Dimens.fontCaption, color = TextSecondary)
            NumBox(yearStr, onYearChange, Modifier.weight(1.4f))
            Text("年", fontSize = Dimens.fontBody, color = OnSurface)
            NumBox(monthStr, onMonthChange, Modifier.weight(0.7f))
            Text("月起", fontSize = Dimens.fontBody, color = OnSurface)
            NumBox(countStr, onCountChange, Modifier.weight(0.7f))
            Text("条", fontSize = Dimens.fontBody, color = OnSurface)
        }

        Button(
            onClick = onSearch, enabled = !loading,
            modifier = Modifier.fillMaxWidth().height(42.dp),
            shape = RoundedCornerShape(21.dp),
            colors = ButtonDefaults.buttonColors(
                containerColor = Primary, contentColor = OnPrimary)
        ) {
            if (loading) {
                CircularProgressIndicator(
                    Modifier.size(18.dp), color = OnPrimary, strokeWidth = 2.dp)
                Spacer(Modifier.width(8.dp))
            }
            Text(if (loading) "搜索中…" else "开始搜索", fontSize = Dimens.fontBody)
        }
    }
}

// ─── 列表 ─────────────────────────────────────────────────

@Composable
private fun SolarList(items: List<SolarEclipseItem>, onClick: (SolarEclipseItem) -> Unit) {
    if (items.isEmpty()) {
        EmptyHint("起始时间范围内没有日食")
        return
    }
    LazyColumn(
        Modifier.fillMaxWidth().padding(horizontal = Dimens.paddingMd),
        verticalArrangement = Arrangement.spacedBy(8.dp),
        contentPadding = PaddingValues(bottom = 24.dp)
    ) {
        items(items) { item -> SolarCard(item, onClick) }
    }
}

@Composable
private fun LunarList(items: List<LunarEclipseItem>, onClick: (LunarEclipseItem) -> Unit) {
    if (items.isEmpty()) {
        EmptyHint("起始时间范围内没有月食")
        return
    }
    LazyColumn(
        Modifier.fillMaxWidth().padding(horizontal = Dimens.paddingMd),
        verticalArrangement = Arrangement.spacedBy(8.dp),
        contentPadding = PaddingValues(bottom = 24.dp)
    ) {
        items(items) { item -> LunarCard(item, onClick) }
    }
}

@Composable
private fun EmptyHint(text: String) {
    Box(
        Modifier.fillMaxWidth().padding(Dimens.paddingXl),
        contentAlignment = Alignment.Center
    ) {
        Text(text, color = TextSecondary, fontSize = Dimens.fontBody)
    }
}

// ─── 卡片 ─────────────────────────────────────────────────

@Composable
private fun SolarCard(item: SolarEclipseItem, onClick: (SolarEclipseItem) -> Unit) {
    val color = EclipseUtil.solarBadgeColor(item.type)
    Column(
        Modifier.fillMaxWidth()
            .shadow(2.dp, RoundedCornerShape(Dimens.radiusMd))
            .clip(RoundedCornerShape(Dimens.radiusMd))
            .background(Surface)
            .clickable { onClick(item) }
            .padding(Dimens.paddingMd),
        verticalArrangement = Arrangement.spacedBy(8.dp)
    ) {
        Row(verticalAlignment = Alignment.CenterVertically) {
            Text(EclipseUtil.solarEmoji(item.type), fontSize = 22.sp)
            Spacer(Modifier.width(8.dp))
            Column(Modifier.weight(1f)) {
                Text("%04d-%02d-%02d %02d:%02d TD"
                    .format(item.year, item.month, item.day, item.hour, item.minute),
                    fontSize = Dimens.fontBody, fontWeight = FontWeight.Bold, color = OnSurface)
                Text(EclipseUtil.jdTdToLocal(item.jd, 8.0, withDate = true) + "  UTC+08",
                    fontSize = 11.sp, color = TextTertiary, fontFamily = FontFamily.Monospace)
            }
            GradientBadge(item.typeName, color)
        }

        Row(verticalAlignment = Alignment.CenterVertically,
            horizontalArrangement = Arrangement.spacedBy(8.dp)) {
            Text("食分", fontSize = Dimens.fontCaption, color = TextSecondary)
            LinearProgressIndicator(
                progress = { EclipseUtil.magnitudeProgress(item.magnitude) },
                modifier = Modifier.weight(1f).height(6.dp).clip(RoundedCornerShape(3.dp)),
                color = color, trackColor = DividerColor
            )
            Text("%.3f".format(item.magnitude),
                fontSize = Dimens.fontCaption, fontWeight = FontWeight.Bold, color = OnSurface)
        }

        Row(horizontalArrangement = Arrangement.spacedBy(16.dp)) {
            ParamSm("γ", "%+.3f".format(item.gamma))
            ParamSm("最长", EclipseUtil.formatDuration(item.duration))
            ParamSm("食带宽",
                if (item.pathWidth > 0) "%.0f km".format(item.pathWidth) else "—")
            if (item.hasCenter) {
                ParamSm("中心点", "%.1f°,%.1f°".format(item.centerLon, item.centerLat))
            }
        }
        Text("点击 → 食带地图 · 本地观测 · 数据",
            fontSize = 11.sp, color = Primary, fontWeight = FontWeight.Medium)
    }
}

@Composable
private fun LunarCard(item: LunarEclipseItem, onClick: (LunarEclipseItem) -> Unit) {
    val color = EclipseUtil.lunarBadgeColor(item.type)
    Column(
        Modifier.fillMaxWidth()
            .shadow(2.dp, RoundedCornerShape(Dimens.radiusMd))
            .clip(RoundedCornerShape(Dimens.radiusMd))
            .background(Surface)
            .clickable { onClick(item) }
            .padding(Dimens.paddingMd),
        verticalArrangement = Arrangement.spacedBy(8.dp)
    ) {
        Row(verticalAlignment = Alignment.CenterVertically) {
            Text(EclipseUtil.lunarEmoji(item.type), fontSize = 22.sp)
            Spacer(Modifier.width(8.dp))
            Column(Modifier.weight(1f)) {
                Text("%04d-%02d-%02d %02d:%02d TD"
                    .format(item.year, item.month, item.day, item.hour, item.minute),
                    fontSize = Dimens.fontBody, fontWeight = FontWeight.Bold, color = OnSurface)
                Text(EclipseUtil.jdTdToLocal(item.jd, 8.0, withDate = true) + "  UTC+08",
                    fontSize = 11.sp, color = TextTertiary, fontFamily = FontFamily.Monospace)
            }
            GradientBadge(item.typeName, color)
        }
        if (item.type != "B") {
            Row(verticalAlignment = Alignment.CenterVertically,
                horizontalArrangement = Arrangement.spacedBy(8.dp)) {
                Text("食分", fontSize = Dimens.fontCaption, color = TextSecondary)
                LinearProgressIndicator(
                    progress = { EclipseUtil.magnitudeProgress(item.magnitude) },
                    modifier = Modifier.weight(1f).height(6.dp).clip(RoundedCornerShape(3.dp)),
                    color = color, trackColor = DividerColor
                )
                Text("%.3f".format(item.magnitude),
                    fontSize = Dimens.fontCaption, fontWeight = FontWeight.Bold, color = OnSurface)
            }
        } else {
            Text("半影食 · 月亮仅穿过地球半影, 不进入本影",
                fontSize = 11.sp, color = TextSecondary)
        }
        Text("点击 → 月食过程动画 · 数据",
            fontSize = 11.sp, color = Primary, fontWeight = FontWeight.Medium)
    }
}

// ─── Solar 详情 Sheet (3 tabs) ─────────────────────────────

@OptIn(ExperimentalMaterial3Api::class)
@Composable
private fun SolarDetailSheet(item: SolarEclipseItem) {
    var tab by remember { mutableStateOf(0) }
    val context = LocalContext.current
    val rootView = LocalView.current
    val scope = rememberCoroutineScope()
    Column(Modifier.fillMaxWidth().padding(horizontal = Dimens.paddingMd, vertical = 8.dp)) {
        // 头部
        Row(verticalAlignment = Alignment.CenterVertically) {
            Text(EclipseUtil.solarEmoji(item.type), fontSize = 26.sp)
            Spacer(Modifier.width(8.dp))
            Column(Modifier.weight(1f)) {
                Text("${item.typeName}日食", fontSize = 16.sp,
                    fontWeight = FontWeight.Bold, color = Primary)
                Text(EclipseUtil.jdTdToLocal(item.jd, 8.0, withDate = true)
                    + " 北京时", fontSize = 11.sp, color = TextSecondary,
                    fontFamily = FontFamily.Monospace)
                Text("Saros ${item.saros} · #${item.sarosMember} · ${item.season}",
                    fontSize = 10.sp, color = TextTertiary)
            }
            GradientBadge(item.typeName, EclipseUtil.solarBadgeColor(item.type))
        }
        Spacer(Modifier.height(8.dp))

        // 分享工具栏
        Row(horizontalArrangement = Arrangement.spacedBy(8.dp),
            verticalAlignment = Alignment.CenterVertically) {
            ShareIconButton("截图分享", Icons.Default.Share) {
                scope.launch {
                    val bmp = EclipseShareUtil.captureView(rootView)
                    val file = withContext(Dispatchers.IO) {
                        EclipseShareUtil.saveBitmap(context, bmp,
                            "solar_${item.year}_${item.month}_${item.day}.png")
                    }
                    EclipseShareUtil.shareImage(context, file)
                }
            }
            ShareIconButton("分享日程 ICS", Icons.Default.Event) {
                scope.launch {
                    val title = "${item.typeName} (Saros ${item.saros})"
                    val desc = """
                        类型: ${item.typeName} (${item.type})
                        食分: %.3f
                        食带宽: ${if (item.pathWidth > 0) "%.0f km".format(item.pathWidth) else "—"}
                        Saros: ${item.saros} (#${item.sarosMember})
                        ${item.season}
                    """.trimIndent().format(item.magnitude)
                    val begin = item.jd - 1.5 / 24
                    val end = item.jd + 1.5 / 24
                    val ics = EclipseShareUtil.buildICS(
                        uid = "solar-${item.year}${item.month}${item.day}@sxwnl",
                        summary = title, description = desc,
                        startJd = begin, endJd = end)
                    val file = withContext(Dispatchers.IO) {
                        EclipseShareUtil.writeICS(context, ics,
                            "solar_${item.year}_${item.month}_${item.day}.ics")
                    }
                    EclipseShareUtil.shareFile(context, file)
                }
            }
        }
        Spacer(Modifier.height(8.dp))

        // tab bar
        TabBar(listOf("🗺  食带地图", "🌞  本地观测", "📊  数据"), tab) { tab = it }
        Spacer(Modifier.height(10.dp))

        when (tab) {
            0 -> SolarMapTab(item)
            1 -> SolarLocalTab(item)
            2 -> SolarDataTab(item)
        }
        Spacer(Modifier.height(24.dp))
    }
}

@Composable
private fun SolarMapTab(item: SolarEclipseItem) {
    val scope = rememberCoroutineScope()
    var path by remember(item.jd) { mutableStateOf<SolarEclipsePath?>(null) }
    var ditu0 by remember { mutableStateOf<DoubleArray>(DoubleArray(0)) }
    var ditu1 by remember { mutableStateOf<DoubleArray>(DoubleArray(0)) }
    var ditu2 by remember { mutableStateOf<DoubleArray>(DoubleArray(0)) }
    var useBig by remember { mutableStateOf(true) }
    var currentJd by remember(item.jd) { mutableStateOf(0.0) }

    LaunchedEffect(item.jd) {
        scope.launch {
            if (ditu0.isEmpty()) ditu0 = CalendarRepository.getWorldMapDitu0()
            if (ditu1.isEmpty()) ditu1 = CalendarRepository.getWorldMapDitu1()
            if (ditu2.isEmpty()) ditu2 = CalendarRepository.getWorldMapDitu2()
            val p = CalendarRepository.getSolarEclipsePath(item.year, item.month, item.day)
            path = p
            if (p != null) currentJd = p.maxEclipseJd
        }
    }

    val p = path
    Column(verticalArrangement = Arrangement.spacedBy(10.dp)) {
        if (p == null) {
            LoadingBox("加载食带数据…")
        } else {
            Row(verticalAlignment = Alignment.CenterVertically,
                horizontalArrangement = Arrangement.spacedBy(8.dp)) {
                Text("地图精度", fontSize = 11.sp, color = TextSecondary)
                Switch(checked = useBig, onCheckedChange = { useBig = it })
                Text(if (useBig) "大图 (ditu1+2)" else "小图 (ditu0)",
                    fontSize = 11.sp, color = TextSecondary)
            }
            val mapToShow = if (useBig && ditu1.isNotEmpty()) ditu1 else ditu0
            val borderToShow = if (useBig && ditu2.isNotEmpty()) ditu2 else null
            SolarPathMapCanvas(mapToShow, borderToShow, p, currentJd)

            // 沿中心线播放
            val pts = p.centerLine
            if (pts.isNotEmpty()) {
                PlaybackController(
                    jdBegin = pts.first().jd,
                    jdEnd = pts.last().jd,
                    currentJd = currentJd,
                    onChange = { currentJd = it }
                )
            }

            EclipseTimelineBar(
                marks = buildSolarMarks(p),
                currentJd = currentJd,
                onSeek = { currentJd = it }
            )

            // 图例
            Column(verticalArrangement = Arrangement.spacedBy(4.dp)) {
                LegendRow(PathCenter, "中心线 (全/环食带核心)")
                LegendRow(PathUmbra, "本影南北界")
                LegendRow(PathPenumbra, "半影南北界 (偏食可见区)")
                LegendRow(MapLand, "海岸线 (内置 ditu0)")
            }
        }
    }
}

@Composable
private fun SolarLocalTab(item: SolarEclipseItem) {
    val scope = rememberCoroutineScope()
    var city by remember { mutableStateOf(Cities.DEFAULT) }
    var local by remember(item.jd, city.name) { mutableStateOf<LocalSolarEclipse?>(null) }
    var loading by remember(item.jd, city.name) { mutableStateOf(true) }
    var currentJd by remember(item.jd, city.name) { mutableStateOf(0.0) }

    LaunchedEffect(item.jd, city.name) {
        loading = true
        scope.launch {
            val r = CalendarRepository.getLocalSolarEclipse(
                item.year, item.month, item.day, city.lon, city.lat, 60)
            local = r
            if (r != null && r.frames.isNotEmpty()) {
                currentJd = r.tMax.takeIf { it > 0 } ?: r.frames.first().jd
            }
            loading = false
        }
    }

    Column(verticalArrangement = Arrangement.spacedBy(10.dp)) {
        LocationStrip(selected = city, onSelect = { city = it })

        val l = local
        when {
            loading -> LoadingBox("计算本地观测帧…")
            l == null || l.frames.isEmpty() -> NoticeBox(
                "在 ${city.name} 看不到此次日食",
                "本地不在月球阴影覆盖范围内, 或日食发生于地平线以下."
            )
            else -> {
                SolarLocalDiscCanvas(
                    frames = l.frames, currentJd = currentJd, tzHours = city.tz)

                PlaybackController(
                    jdBegin = l.frames.first().jd,
                    jdEnd = l.frames.last().jd,
                    currentJd = currentJd,
                    onChange = { currentJd = it }
                )

                EclipseTimelineBar(
                    marks = buildLocalSolarMarks(l),
                    currentJd = currentJd,
                    onSeek = { currentJd = it },
                    tzHours = city.tz
                )

                ParamGrid(buildList {
                    add("当地最大食分" to "%.3f".format(l.maxMagnitude))
                    add("月日视径比" to "%.3f".format(l.moonSunRatio))
                    add("日出" to EclipseUtil.jdTdToLocal(l.tSunrise, city.tz))
                    add("日没" to EclipseUtil.jdTdToLocal(l.tSunset, city.tz))
                    add("类型" to l.type.ifEmpty { "—" })
                    add("帧数" to "${l.frames.size}")
                })
            }
        }
    }
}

@Composable
private fun SolarDataTab(item: SolarEclipseItem) {
    val scope = rememberCoroutineScope()
    var path by remember(item.jd) { mutableStateOf<SolarEclipsePath?>(null) }
    LaunchedEffect(item.jd) {
        scope.launch {
            path = CalendarRepository.getSolarEclipsePath(item.year, item.month, item.day)
        }
    }
    Column(
        Modifier.verticalScroll(rememberScrollState()),
        verticalArrangement = Arrangement.spacedBy(10.dp)
    ) {
        ParamGrid(listOf(
            "日期" to "%04d-%02d-%02d".format(item.year, item.month, item.day),
            "TD 时刻" to "%02d:%02d".format(item.hour, item.minute),
            "北京时" to EclipseUtil.jdTdToLocal(item.jd, 8.0, withDate = false),
            "食分" to "%.4f".format(item.magnitude),
            "γ" to "%+.4f".format(item.gamma),
            "最长时长" to EclipseUtil.formatDuration(item.duration),
            "食带宽" to if (item.pathWidth > 0) "%.1f km".format(item.pathWidth) else "—",
            "类型代码" to item.type,
            "中心点" to if (item.hasCenter)
                "${EclipseUtil.lonLabel(item.centerLon)}, ${EclipseUtil.latLabel(item.centerLat)}"
            else "无中心",
            "Saros 周期" to "#${item.saros}  ·  序号 ${item.sarosMember}",
            "食季" to item.season.ifEmpty { "—" }
        ))

        path?.let { p ->
            HorizontalDivider(color = DividerColor)
            Text("全球路径关键点 (TD)",
                fontSize = Dimens.fontSubtitle, fontWeight = FontWeight.Bold, color = Primary)
            EventPoint("偏食始 P₁", p.partialStart)
            EventPoint("中心始 C₁", p.centralStart)
            EventPoint("食甚 Greatest", p.greatestEclipse)
            EventPoint("中心终 C₂", p.centralEnd)
            EventPoint("偏食终 P₄", p.partialEnd)

            Row(horizontalArrangement = Arrangement.spacedBy(16.dp)) {
                ParamSm("中心线", "${p.centerLine.size} 点")
                ParamSm("本影界", "${p.umbraNorth.size + p.umbraSouth.size} 点")
                ParamSm("半影界", "${p.penumbraNorth.size + p.penumbraSouth.size} 点")
            }
        }
    }
}

private fun buildSolarMarks(p: SolarEclipsePath): List<TimelineMark> = listOfNotNull(
    p.partialStart.takeIf { it.jd > 0 }?.let {
        TimelineMark("P1", it.jd, PathPenumbra) },
    p.centralStart.takeIf { it.jd > 0 && it.longitude < 99 }?.let {
        TimelineMark("C1", it.jd, PathUmbra) },
    p.greatestEclipse.takeIf { it.jd > 0 }?.let {
        TimelineMark("MAX", it.jd, PathCenter) },
    p.centralEnd.takeIf { it.jd > 0 && it.longitude < 99 }?.let {
        TimelineMark("C2", it.jd, PathUmbra) },
    p.partialEnd.takeIf { it.jd > 0 }?.let {
        TimelineMark("P4", it.jd, PathPenumbra) }
)

private fun buildLocalSolarMarks(l: LocalSolarEclipse): List<TimelineMark> = listOfNotNull(
    l.tC1.takeIf { it > 0 }?.let { TimelineMark("初亏", it, PathPenumbra) },
    l.tC2.takeIf { it > 0 }?.let { TimelineMark("食既", it, PathUmbra) },
    l.tMax.takeIf { it > 0 }?.let { TimelineMark("食甚", it, PathCenter) },
    l.tC3.takeIf { it > 0 }?.let { TimelineMark("生光", it, PathUmbra) },
    l.tC4.takeIf { it > 0 }?.let { TimelineMark("复圆", it, PathPenumbra) }
)

// ─── Lunar 详情 Sheet (2 tabs) ─────────────────────────────

@Composable
private fun LunarDetailSheet(item: LunarEclipseItem) {
    var tab by remember { mutableStateOf(0) }
    val context = LocalContext.current
    val rootView = LocalView.current
    val scope = rememberCoroutineScope()
    Column(Modifier.fillMaxWidth().padding(horizontal = Dimens.paddingMd, vertical = 8.dp)) {
        Row(verticalAlignment = Alignment.CenterVertically) {
            Text(EclipseUtil.lunarEmoji(item.type), fontSize = 26.sp)
            Spacer(Modifier.width(8.dp))
            Column(Modifier.weight(1f)) {
                Text("${item.typeName}月食", fontSize = 16.sp,
                    fontWeight = FontWeight.Bold, color = Primary)
                Text(EclipseUtil.jdTdToLocal(item.jd, 8.0, withDate = true)
                    + " 北京时", fontSize = 11.sp, color = TextSecondary,
                    fontFamily = FontFamily.Monospace)
                Text("Saros ${item.saros} · #${item.sarosMember} · ${item.season}",
                    fontSize = 10.sp, color = TextTertiary)
            }
            GradientBadge(item.typeName, EclipseUtil.lunarBadgeColor(item.type))
        }
        Spacer(Modifier.height(8.dp))

        Row(horizontalArrangement = Arrangement.spacedBy(8.dp)) {
            ShareIconButton("截图分享", Icons.Default.Share) {
                scope.launch {
                    val bmp = EclipseShareUtil.captureView(rootView)
                    val file = withContext(Dispatchers.IO) {
                        EclipseShareUtil.saveBitmap(context, bmp,
                            "lunar_${item.year}_${item.month}_${item.day}.png")
                    }
                    EclipseShareUtil.shareImage(context, file)
                }
            }
            ShareIconButton("分享日程 ICS", Icons.Default.Event) {
                scope.launch {
                    val title = "${item.typeName} (Saros ${item.saros})"
                    val desc = """
                        类型: ${item.typeName} (${item.type})
                        食分: %.3f
                        Saros: ${item.saros} (#${item.sarosMember})
                        ${item.season}
                    """.trimIndent().format(item.magnitude)
                    val begin = item.jd - 2.0 / 24
                    val end = item.jd + 2.0 / 24
                    val ics = EclipseShareUtil.buildICS(
                        uid = "lunar-${item.year}${item.month}${item.day}@sxwnl",
                        summary = title, description = desc,
                        startJd = begin, endJd = end)
                    val file = withContext(Dispatchers.IO) {
                        EclipseShareUtil.writeICS(context, ics,
                            "lunar_${item.year}_${item.month}_${item.day}.ics")
                    }
                    EclipseShareUtil.shareFile(context, file)
                }
            }
        }
        Spacer(Modifier.height(8.dp))

        TabBar(listOf("🌙  过程动画", "📊  数据"), tab) { tab = it }
        Spacer(Modifier.height(10.dp))

        when (tab) {
            0 -> LunarAnimTab(item)
            1 -> LunarDataTab(item)
        }
        Spacer(Modifier.height(24.dp))
    }
}

@Composable
private fun LunarAnimTab(item: LunarEclipseItem) {
    val scope = rememberCoroutineScope()
    var detail by remember(item.jd) { mutableStateOf<LunarEclipseDetail?>(null) }
    var currentJd by remember(item.jd) { mutableStateOf(0.0) }
    var loading by remember(item.jd) { mutableStateOf(true) }

    LaunchedEffect(item.jd) {
        loading = true
        scope.launch {
            val d = CalendarRepository.getLunarEclipseDetail(
                item.year, item.month, item.day, 120)
            detail = d
            if (d != null && d.frames.isNotEmpty()) {
                currentJd = d.tMax.takeIf { it > 0 } ?: d.frames.first().jd
            }
            loading = false
        }
    }

    val d = detail
    Column(verticalArrangement = Arrangement.spacedBy(10.dp)) {
        when {
            loading -> LoadingBox("加载月食帧…")
            d == null || d.frames.isEmpty() -> NoticeBox("无月食帧数据", "")
            else -> {
                LunarDiscCanvas(frames = d.frames, currentJd = currentJd, tzHours = 8.0)

                PlaybackController(
                    jdBegin = d.frames.first().jd,
                    jdEnd = d.frames.last().jd,
                    currentJd = currentJd,
                    onChange = { currentJd = it }
                )

                EclipseTimelineBar(
                    marks = buildLunarMarks(d),
                    currentJd = currentJd,
                    onSeek = { currentJd = it }
                )
            }
        }
    }
}

@Composable
private fun LunarDataTab(item: LunarEclipseItem) {
    val scope = rememberCoroutineScope()
    var d by remember(item.jd) { mutableStateOf<LunarEclipseDetail?>(null) }
    LaunchedEffect(item.jd) {
        scope.launch {
            d = CalendarRepository.getLunarEclipseDetail(item.year, item.month, item.day, 120)
        }
    }
    Column(
        Modifier.verticalScroll(rememberScrollState()),
        verticalArrangement = Arrangement.spacedBy(10.dp)
    ) {
        ParamGrid(listOf(
            "日期" to "%04d-%02d-%02d".format(item.year, item.month, item.day),
            "TD 时刻" to "%02d:%02d".format(item.hour, item.minute),
            "北京时" to EclipseUtil.jdTdToLocal(item.jd, 8.0, withDate = false),
            "类型" to item.typeName,
            "食分" to if (item.type == "B") "—" else "%.4f".format(item.magnitude),
            "Saros 周期" to "#${item.saros}  ·  序号 ${item.sarosMember}",
            "食季" to item.season.ifEmpty { "—" }
        ))
        d?.let { detail ->
            HorizontalDivider(color = DividerColor)
            Text("过程时刻 (TD)", fontSize = Dimens.fontSubtitle,
                fontWeight = FontWeight.Bold, color = Primary)
            listOf(
                "半影食始 P₁" to detail.tP1,
                "初亏 U₁" to detail.tU1,
                "食既 U₂" to detail.tU2,
                "食甚 Max" to detail.tMax,
                "生光 U₃" to detail.tU3,
                "复圆 U₄" to detail.tU4,
                "半影食终 P₄" to detail.tP4
            ).forEach { (label, jd) ->
                if (jd > 0) TimePointRow(label, jd, 8.0)
            }
        }
    }
}

private fun buildLunarMarks(d: LunarEclipseDetail): List<TimelineMark> = listOfNotNull(
    d.tP1.takeIf { it > 0 }?.let { TimelineMark("P1", it, PathPenumbra) },
    d.tU1.takeIf { it > 0 }?.let { TimelineMark("U1", it, PathUmbra) },
    d.tU2.takeIf { it > 0 }?.let { TimelineMark("U2", it, MoonBlood) },
    d.tMax.takeIf { it > 0 }?.let { TimelineMark("MAX", it, PathCenter) },
    d.tU3.takeIf { it > 0 }?.let { TimelineMark("U3", it, MoonBlood) },
    d.tU4.takeIf { it > 0 }?.let { TimelineMark("U4", it, PathUmbra) },
    d.tP4.takeIf { it > 0 }?.let { TimelineMark("P4", it, PathPenumbra) }
)

// ─── 通用原子 ─────────────────────────────────────────────

@Composable
private fun TabBar(tabs: List<String>, selected: Int, onSelect: (Int) -> Unit) {
    Row(
        Modifier.fillMaxWidth().clip(RoundedCornerShape(50))
            .background(Background)
            .padding(4.dp)
    ) {
        tabs.forEachIndexed { i, label ->
            val sel = selected == i
            Box(
                Modifier.weight(1f).height(36.dp)
                    .clip(RoundedCornerShape(50))
                    .background(if (sel) Primary else Color.Transparent)
                    .clickable { onSelect(i) },
                contentAlignment = Alignment.Center
            ) {
                Text(label, color = if (sel) OnPrimary else TextSecondary,
                    fontSize = 12.sp,
                    fontWeight = if (sel) FontWeight.Bold else FontWeight.Medium)
            }
        }
    }
}

@Composable
private fun GradientBadge(text: String, color: Color) {
    Box(
        Modifier.clip(RoundedCornerShape(50))
            .background(Brush.horizontalGradient(listOf(color, color.copy(alpha = 0.75f))))
            .padding(horizontal = 10.dp, vertical = 3.dp)
    ) {
        Text(text, color = Color.White, fontSize = 11.sp, fontWeight = FontWeight.Bold)
    }
}

@Composable
private fun ParamSm(label: String, value: String) {
    Column {
        Text(label, fontSize = 10.sp, color = TextTertiary)
        Text(value, fontSize = 12.sp, color = OnSurface, fontWeight = FontWeight.Medium)
    }
}

@Composable
private fun ParamGrid(pairs: List<Pair<String, String>>) {
    Column(verticalArrangement = Arrangement.spacedBy(4.dp)) {
        pairs.chunked(2).forEach { row ->
            Row(horizontalArrangement = Arrangement.spacedBy(10.dp)) {
                row.forEach { (k, v) ->
                    Row(Modifier.weight(1f).clip(RoundedCornerShape(6.dp))
                        .background(Background)
                        .padding(horizontal = 8.dp, vertical = 5.dp),
                        verticalAlignment = Alignment.CenterVertically) {
                        Text(k, fontSize = 11.sp, color = TextSecondary,
                            modifier = Modifier.width(82.dp))
                        Text(v, fontSize = 12.sp, color = OnSurface,
                            fontWeight = FontWeight.Medium,
                            fontFamily = FontFamily.Monospace)
                    }
                }
                if (row.size == 1) Spacer(Modifier.weight(1f))
            }
        }
    }
}

@Composable
private fun EventPoint(label: String, p: EclipseGeoPoint) {
    Row(verticalAlignment = Alignment.CenterVertically,
        modifier = Modifier.fillMaxWidth().clip(RoundedCornerShape(6.dp))
            .background(Background).padding(horizontal = 10.dp, vertical = 6.dp)) {
        Text(label, fontSize = Dimens.fontCaption, color = TextSecondary,
            modifier = Modifier.width(110.dp))
        Column(Modifier.weight(1f)) {
            Text(if (p.longitude < 99) "%.2f°,  %.2f°".format(p.longitude, p.latitude) else "—",
                fontSize = Dimens.fontCaption, color = OnSurface,
                fontFamily = FontFamily.Monospace)
            Text(EclipseUtil.jdToDateTime(p.jd),
                fontSize = 11.sp, color = TextTertiary, fontFamily = FontFamily.Monospace)
        }
    }
}

@Composable
private fun TimePointRow(label: String, jd: Double, tzHours: Double) {
    Row(verticalAlignment = Alignment.CenterVertically,
        modifier = Modifier.fillMaxWidth().clip(RoundedCornerShape(6.dp))
            .background(Background).padding(horizontal = 10.dp, vertical = 6.dp)) {
        Text(label, fontSize = Dimens.fontCaption, color = TextSecondary,
            modifier = Modifier.width(110.dp))
        Column(Modifier.weight(1f)) {
            Text(EclipseUtil.jdToDateTime(jd),
                fontSize = Dimens.fontCaption, color = OnSurface,
                fontFamily = FontFamily.Monospace)
            Text(EclipseUtil.jdTdToLocal(jd, tzHours, withDate = true) + "  " +
                EclipseUtil.tzLabel(tzHours),
                fontSize = 11.sp, color = TextTertiary, fontFamily = FontFamily.Monospace)
        }
    }
}

@Composable
private fun ShareIconButton(
    label: String,
    icon: androidx.compose.ui.graphics.vector.ImageVector,
    onClick: () -> Unit
) {
    Row(
        Modifier.clip(RoundedCornerShape(8.dp))
            .background(Background)
            .clickable { onClick() }
            .padding(horizontal = 10.dp, vertical = 6.dp),
        verticalAlignment = Alignment.CenterVertically,
        horizontalArrangement = Arrangement.spacedBy(6.dp)
    ) {
        Icon(icon, null, Modifier.size(14.dp), tint = Primary)
        Text(label, fontSize = 11.sp, color = Primary, fontWeight = FontWeight.SemiBold)
    }
}

@Composable
private fun LegendRow(color: Color, label: String) {
    Row(verticalAlignment = Alignment.CenterVertically,
        horizontalArrangement = Arrangement.spacedBy(8.dp)) {
        Box(Modifier.size(width = 22.dp, height = 4.dp)
            .clip(RoundedCornerShape(2.dp)).background(color))
        Text(label, fontSize = 11.sp, color = TextSecondary)
    }
}

@Composable
private fun LoadingBox(text: String) {
    Box(Modifier.fillMaxWidth().padding(24.dp), contentAlignment = Alignment.Center) {
        Column(horizontalAlignment = Alignment.CenterHorizontally,
            verticalArrangement = Arrangement.spacedBy(8.dp)) {
            CircularProgressIndicator(Modifier.size(24.dp), strokeWidth = 2.dp)
            Text(text, fontSize = 12.sp, color = TextSecondary)
        }
    }
}

@Composable
private fun NoticeBox(title: String, subtitle: String) {
    Column(
        Modifier.fillMaxWidth().clip(RoundedCornerShape(8.dp))
            .background(Background).padding(16.dp),
        horizontalAlignment = Alignment.CenterHorizontally,
        verticalArrangement = Arrangement.spacedBy(4.dp)
    ) {
        Text("⊘  $title", color = TextSecondary, fontSize = 14.sp,
            fontWeight = FontWeight.Medium)
        if (subtitle.isNotEmpty()) {
            Text(subtitle, color = TextTertiary, fontSize = 11.sp,
                textAlign = TextAlign.Center)
        }
    }
}

@Composable
private fun NumBox(value: String, onChange: (String) -> Unit, modifier: Modifier = Modifier) {
    Box(
        modifier.height(38.dp)
            .clip(RoundedCornerShape(Dimens.radiusSm))
            .background(Background)
            .border(BorderStroke(1.dp, DividerColor), RoundedCornerShape(Dimens.radiusSm)),
        contentAlignment = Alignment.Center
    ) {
        BasicTextField(
            value = value, onValueChange = onChange, singleLine = true,
            textStyle = TextStyle(color = OnSurface, fontSize = Dimens.fontBody,
                textAlign = TextAlign.Center),
            modifier = Modifier.fillMaxWidth().padding(horizontal = 6.dp),
            keyboardOptions = KeyboardOptions(keyboardType = KeyboardType.Number),
            cursorBrush = SolidColor(Primary)
        )
    }
}
