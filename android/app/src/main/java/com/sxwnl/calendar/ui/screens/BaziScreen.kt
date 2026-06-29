package com.sxwnl.calendar.ui.screens

import android.widget.Toast
import androidx.compose.foundation.background
import androidx.compose.foundation.border
import androidx.compose.foundation.clickable
import androidx.compose.foundation.layout.*
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
import androidx.compose.ui.graphics.SolidColor
import androidx.compose.ui.platform.LocalContext
import androidx.compose.ui.text.TextStyle
import androidx.compose.ui.text.font.FontWeight
import androidx.compose.ui.text.input.KeyboardType
import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.sp
import com.sxwnl.calendar.data.BaziParams
import com.sxwnl.calendar.data.CalendarRepository
import com.sxwnl.calendar.data.Cities
import com.sxwnl.calendar.data.LunarMonth
import com.sxwnl.calendar.ui.components.BaziDateTimePicker
import com.sxwnl.calendar.ui.components.BaziResultArg
import com.sxwnl.calendar.ui.components.BaziResultView
import com.sxwnl.calendar.ui.components.CityPickerSheet
import com.sxwnl.calendar.ui.theme.*
import com.sxwnl.calendar.util.BaziCalc
import com.sxwnl.calendar.util.BaziCalc.BirthInputMode
import com.sxwnl.calendar.util.BaziCalc.DateSelection
import kotlinx.coroutines.launch
import androidx.compose.material3.Switch
import androidx.compose.material3.SwitchDefaults
import kotlin.math.abs

// ════════════════════════════════════════════════════════════════
//  BaziScreen — 与鸿蒙 BaziPage.ets 对齐
//
//  顶部 banner + 输入卡片 + 「排盘」按钮 → 弹出 BaziResultView
// ════════════════════════════════════════════════════════════════

@Composable
fun BaziScreen() {
    val ctx = LocalContext.current
    val scope = rememberCoroutineScope()

    var name by remember { mutableStateOf("") }
    var gender by remember { mutableIntStateOf(0) }   // 0 男 1 女
    var lifa by remember { mutableIntStateOf(BaziCalc.LIFA_DING_QI) }
    var astEnabled by remember { mutableStateOf(false) }

    // 出生地: 默认「北京市 / 天安门」, 手动开关关闭=省/区县两级下拉, 打开=手动输入经纬度
    val defaultCity = Cities.DEFAULT
    var manualGeo by remember { mutableStateOf(false) }
    var regionIdx by remember { mutableIntStateOf(0) }
    var cityIdx by remember { mutableIntStateOf(0) }
    var longitude by remember { mutableDoubleStateOf(defaultCity.lon) }
    var latitude by remember { mutableDoubleStateOf(defaultCity.lat) }
    var lonText by remember { mutableStateOf(defaultCity.lon.toString()) }
    var latText by remember { mutableStateOf(defaultCity.lat.toString()) }

    var sel by remember {
        mutableStateOf(DateSelection(
            inputMode = BirthInputMode.SOLAR,
            year = 1990, month = 1, day = 1,
            hour = 12, minute = 0))
    }
    var picked by remember { mutableStateOf(false) }
    var lunarMonths by remember { mutableStateOf(emptyList<LunarMonth>()) }

    var showPicker by remember { mutableStateOf(false) }
    var resultArg by remember { mutableStateOf<BaziResultArg?>(null) }

    LaunchedEffect(Unit) {
        lunarMonths = CalendarRepository.getLunarMonths(sel.year)
    }

    Column(
        Modifier
            .fillMaxSize()
            .background(Background)
            .verticalScroll(rememberScrollState())
            .padding(bottom = 32.dp)
    ) {
        // 顶部 banner
        Column(
            Modifier
                .fillMaxWidth()
                .background(
                    Brush.linearGradient(
                        colors = listOf(GradientStart, GradientEnd)
                    )
                )
                .padding(top = Dimens.paddingXl, bottom = Dimens.paddingLg),
            horizontalAlignment = Alignment.CenterHorizontally
        ) {
            Text("☯ 八字排盘",
                fontSize = Dimens.fontTitle, fontWeight = FontWeight.Bold,
                color = OnPrimary)
        }

        // 输入卡片
        Column(
            Modifier
                .padding(Dimens.paddingMd)
                .fillMaxWidth()
                .shadow(6.dp, RoundedCornerShape(Dimens.radiusLg))
                .clip(RoundedCornerShape(Dimens.radiusLg))
                .background(Surface)
                .padding(Dimens.paddingLg),
            verticalArrangement = Arrangement.spacedBy(Dimens.paddingMd)
        ) {
            // 姓名 + 性别
            Row(verticalAlignment = Alignment.CenterVertically) {
                CompactInput(
                    value = name, onValueChange = { name = it },
                    placeholder = "姓名(选填)",
                    modifier = Modifier
                        .weight(1f)
                        .height(44.dp)
                        .clip(RoundedCornerShape(Dimens.radiusSm))
                        .background(Background)
                )
                Row(
                    Modifier
                        .padding(start = 10.dp)
                        .clip(RoundedCornerShape(Dimens.radiusSm))
                        .border(1.dp, DividerColor, RoundedCornerShape(Dimens.radiusSm))
                ) {
                    GenderChip("男", 0, gender) { gender = 0 }
                    GenderChip("女", 1, gender) { gender = 1 }
                }
            }

            // 出生时间
            Row(
                Modifier
                    .fillMaxWidth()
                    .height(48.dp)
                    .clip(RoundedCornerShape(Dimens.radiusSm))
                    .background(Background)
                    .clickable { showPicker = true }
                    .padding(horizontal = 14.dp),
                verticalAlignment = Alignment.CenterVertically
            ) {
                Text(
                    recordText(picked, sel, lunarMonths),
                    fontSize = Dimens.fontSubtitle,
                    color = if (picked) OnSurface else TextSecondary,
                    modifier = Modifier.weight(1f)
                )
                Text("▾", fontSize = Dimens.fontBody, color = TextSecondary)
            }

            // 出生地: 省/区县联动弹窗(默认) 或 手动输入经纬度(开关切换)
            BirthPlaceSection(
                manualGeo = manualGeo,
                onManualGeoChange = { manualGeo = it },
                regionIdx = regionIdx, cityIdx = cityIdx,
                longitude = longitude, latitude = latitude,
                lonText = lonText, latText = latText,
                onCitySelected = { r, c, city ->
                    regionIdx = r; cityIdx = c
                    longitude = city.lon; latitude = city.lat
                    lonText = city.lon.toString(); latText = city.lat.toString()
                },
                onLonTextChange = { lonText = it; it.toDoubleOrNull()?.let { v -> longitude = v } },
                onLatTextChange = { latText = it; it.toDoubleOrNull()?.let { v -> latitude = v } }
            )

            // 历法
            Row(verticalAlignment = Alignment.CenterVertically) {
                Text("历法", fontSize = Dimens.fontCaption, color = TextSecondary,
                    modifier = Modifier.width(48.dp))
                Row(Modifier.weight(1f)) {
                    BaziCalc.LIFA_OPTIONS.forEachIndexed { i, opt ->
                        val sel0 = lifa == opt.value
                        Text(
                            opt.label,
                            fontSize = Dimens.fontCaption,
                            fontWeight = if (sel0) FontWeight.Bold else FontWeight.Normal,
                            color = if (sel0) OnPrimary else TextSecondary,
                            modifier = Modifier
                                .weight(1f)
                                .padding(end = if (i < BaziCalc.LIFA_OPTIONS.size - 1) 6.dp else 0.dp)
                                .height(34.dp)
                                .clip(RoundedCornerShape(Dimens.radiusSm))
                                .background(if (sel0) Primary else Background)
                                .clickable { lifa = opt.value }
                                .wrapContentSize()
                        )
                    }
                }
            }

            // 真太阳时
            Row(verticalAlignment = Alignment.CenterVertically) {
                Checkbox(
                    checked = astEnabled, onCheckedChange = { astEnabled = it },
                    colors = CheckboxDefaults.colors(checkedColor = Primary),
                    modifier = Modifier.size(20.dp)
                )
                Text("使用真太阳时",
                    fontSize = Dimens.fontBody, color = OnSurface,
                    modifier = Modifier.padding(start = 6.dp))
            }

            // 排盘按钮
            Button(
                onClick = {
                    scope.launch {
                        // 八字时间一律按北京时间处理, 经纬度仅用于真太阳时方程修正
                        val result = CalendarRepository.calcBazi(
                            BaziParams(
                                name = name.ifEmpty { "匿名" },
                                gender = gender == 1,
                                lifa = lifa,
                                astEnabled = astEnabled,
                                longitude = longitude, latitude = latitude,
                                inputMode = sel.inputMode.rawValue,
                                year = sel.year, month = sel.month, day = sel.day,
                                hour = sel.hour, minute = sel.minute,
                                isLeap = if (sel.inputMode == BirthInputMode.LUNAR)
                                    sel.isLeap else false,
                                isSpec = if (sel.inputMode == BirthInputMode.LUNAR)
                                    sel.isSpec else false
                            )
                        )
                        if (result == null) {
                            Toast.makeText(ctx, "排盘失败", Toast.LENGTH_SHORT).show()
                        } else {
                            val lifaLabel = BaziCalc.LIFA_OPTIONS
                                .firstOrNull { it.value == lifa }?.label ?: "定气"
                            resultArg = BaziResultArg(
                                result = result,
                                birthYear = sel.year,
                                astEnabled = astEnabled,
                                longitude = longitude, latitude = latitude,
                                lifaLabel = lifaLabel
                            )
                        }
                    }
                },
                modifier = Modifier
                    .fillMaxWidth()
                    .height(48.dp),
                shape = RoundedCornerShape(Dimens.radiusLg),
                colors = ButtonDefaults.buttonColors(
                    containerColor = Accent, contentColor = PrimaryDark
                )
            ) {
                Text("排  盘",
                    fontSize = Dimens.fontSubtitle, fontWeight = FontWeight.Bold)
            }
        }
    }

    if (showPicker) {
        BaziDateTimePicker(
            initial = sel,
            onDismiss = { showPicker = false },
            onConfirm = { newSel ->
                sel = newSel
                scope.launch {
                    lunarMonths = CalendarRepository.getLunarMonths(newSel.year)
                }
                picked = true
            }
        )
    }

    val ra = resultArg
    if (ra != null) {
        BaziResultView(arg = ra, onDismiss = { resultArg = null })
    }
}

@Composable
private fun GenderChip(
    label: String, value: Int, current: Int, onClick: () -> Unit
) {
    val active = current == value
    Box(
        Modifier
            .width(48.dp).height(44.dp)
            .background(if (active) Primary else Surface)
            .clickable(onClick = onClick),
        contentAlignment = Alignment.Center
    ) {
        Text(label,
            fontSize = Dimens.fontBody,
            fontWeight = if (active) FontWeight.Bold else FontWeight.Normal,
            color = if (active) OnPrimary else TextSecondary)
    }
}

@Composable
private fun CompactInput(
    value: String,
    onValueChange: (String) -> Unit,
    modifier: Modifier = Modifier,
    placeholder: String = "",
    keyboardOptions: KeyboardOptions = KeyboardOptions.Default
) {
    Box(modifier.padding(horizontal = 10.dp), contentAlignment = Alignment.CenterStart) {
        if (placeholder.isNotEmpty() && value.isEmpty()) {
            Text(placeholder, color = TextSecondary, fontSize = Dimens.fontBody)
        }
        BasicTextField(
            value = value, onValueChange = onValueChange,
            singleLine = true,
            textStyle = TextStyle(color = OnSurface, fontSize = Dimens.fontBody),
            modifier = Modifier.fillMaxWidth(),
            keyboardOptions = keyboardOptions,
            cursorBrush = SolidColor(OnSurface)
        )
    }
}

private fun recordText(
    picked: Boolean, sel: DateSelection, months: List<LunarMonth>
): String {
    if (!picked) return "请选择出生时间"
    val p = BaziParams(
        inputMode = sel.inputMode.rawValue,
        year = sel.year, month = sel.month, day = sel.day,
        hour = sel.hour, minute = sel.minute,
        isLeap = sel.isLeap, isSpec = sel.isSpec
    )
    return BaziCalc.formatRecord(p, months)
}

// ─── 出生地选择 ────────────────────────────────────────────
//   manualGeo=false: 整行点击 → 底部弹出两栏(省/区县)联动选择器 CityPickerSheet
//   manualGeo=true:  直接输入经纬度
//   时间一律按北京时间处理, 经纬度仅用于真太阳时方程修正
@Composable
private fun BirthPlaceSection(
    manualGeo: Boolean,
    onManualGeoChange: (Boolean) -> Unit,
    regionIdx: Int,
    cityIdx: Int,
    longitude: Double, latitude: Double,
    lonText: String, latText: String,
    onCitySelected: (regionIdx: Int, cityIdx: Int, city: com.sxwnl.calendar.data.City) -> Unit,
    onLonTextChange: (String) -> Unit,
    onLatTextChange: (String) -> Unit
) {
    var showCityPicker by remember { mutableStateOf(false) }

    Column(verticalArrangement = Arrangement.spacedBy(6.dp)) {
        Row(verticalAlignment = Alignment.CenterVertically) {
            Text("出生地", fontSize = Dimens.fontCaption, color = TextSecondary,
                modifier = Modifier.width(48.dp))
            Text(
                birthPlaceSummary(manualGeo, regionIdx, cityIdx, longitude, latitude),
                fontSize = Dimens.fontSmall, color = TextSecondary,
                maxLines = 1,
                modifier = Modifier.weight(1f)
            )
            Text("手动", fontSize = Dimens.fontSmall, color = TextSecondary,
                modifier = Modifier.padding(end = 4.dp))
            Switch(
                checked = manualGeo,
                onCheckedChange = onManualGeoChange,
                colors = SwitchDefaults.colors(checkedThumbColor = Primary)
            )
        }
        if (manualGeo) {
            Row(verticalAlignment = Alignment.CenterVertically) {
                Text("经度", fontSize = Dimens.fontSmall, color = TextSecondary)
                CompactInput(
                    value = lonText, onValueChange = onLonTextChange,
                    keyboardOptions = KeyboardOptions(keyboardType = KeyboardType.Decimal),
                    modifier = Modifier
                        .padding(start = 4.dp, end = 8.dp)
                        .weight(1f).height(38.dp)
                        .clip(RoundedCornerShape(Dimens.radiusSm))
                        .background(Background)
                )
                Text("纬度", fontSize = Dimens.fontSmall, color = TextSecondary)
                CompactInput(
                    value = latText, onValueChange = onLatTextChange,
                    keyboardOptions = KeyboardOptions(keyboardType = KeyboardType.Decimal),
                    modifier = Modifier
                        .padding(start = 4.dp)
                        .weight(1f).height(38.dp)
                        .clip(RoundedCornerShape(Dimens.radiusSm))
                        .background(Background)
                )
            }
        } else {
            // 触发器: 点击整行打开两栏弹窗
            Row(
                Modifier
                    .fillMaxWidth().height(38.dp)
                    .clip(RoundedCornerShape(Dimens.radiusSm))
                    .background(Background)
                    .clickable { showCityPicker = true }
                    .padding(horizontal = 14.dp),
                verticalAlignment = Alignment.CenterVertically
            ) {
                Text(
                    regionCityLabel(regionIdx, cityIdx),
                    fontSize = Dimens.fontCaption,
                    color = OnSurface,
                    maxLines = 1,
                    modifier = Modifier.weight(1f)
                )
                Text("▾", fontSize = Dimens.fontSmall, color = TextSecondary)
            }
        }
    }

    if (showCityPicker) {
        CityPickerSheet(
            initRegionIdx = regionIdx,
            initCityIdx = cityIdx,
            onConfirm = { r, c, city ->
                onCitySelected(r, c, city)
                showCityPicker = false
            },
            onDismiss = { showCityPicker = false }
        )
    }
}

private fun regionCityLabel(regionIdx: Int, cityIdx: Int): String {
    val r = Cities.REGIONS[regionIdx]
    val c = r.cities[cityIdx]
    return "${r.name}  ▸  ${c.name}"
}

private fun birthPlaceSummary(
    manualGeo: Boolean, regionIdx: Int, cityIdx: Int,
    longitude: Double, latitude: Double
): String {
    val ew = if (longitude >= 0) "E" else "W"
    val ns = if (latitude >= 0) "N" else "S"
    val lon = "%.2f".format(abs(longitude))
    val lat = "%.2f".format(abs(latitude))
    if (manualGeo) return "$lon°$ew, $lat°$ns"
    val r = Cities.REGIONS[regionIdx]
    val c = r.cities[cityIdx]
    return "${r.name} · ${c.name} · $lon°$ew, $lat°$ns"
}
