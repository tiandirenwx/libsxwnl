package com.sxwnl.calendar.ui.components

import androidx.compose.foundation.background
import androidx.compose.foundation.clickable
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.lazy.LazyColumn
import androidx.compose.foundation.lazy.itemsIndexed
import androidx.compose.foundation.lazy.rememberLazyListState
import androidx.compose.material3.ExperimentalMaterial3Api
import androidx.compose.material3.HorizontalDivider
import androidx.compose.material3.ModalBottomSheet
import androidx.compose.material3.Text
import androidx.compose.material3.TextButton
import androidx.compose.material3.rememberModalBottomSheetState
import androidx.compose.runtime.Composable
import androidx.compose.runtime.LaunchedEffect
import androidx.compose.runtime.getValue
import androidx.compose.runtime.mutableIntStateOf
import androidx.compose.runtime.remember
import androidx.compose.runtime.setValue
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.text.font.FontWeight
import androidx.compose.ui.text.style.TextAlign
import androidx.compose.ui.unit.dp
import com.sxwnl.calendar.data.Cities
import com.sxwnl.calendar.data.City
import com.sxwnl.calendar.ui.theme.*
import kotlin.math.abs

// ════════════════════════════════════════════════════════════════
//  CityPickerSheet — 出生地选择 (与鸿蒙 CityPickerDialog 对齐)
//
//  ModalBottomSheet + Row { LazyColumn(省) | LazyColumn(区县) }
//  左栏选省时联动右栏, 顶部「确定」回调最终结果, 「取消」直接关闭.
// ════════════════════════════════════════════════════════════════

@OptIn(ExperimentalMaterial3Api::class)
@Composable
fun CityPickerSheet(
    initRegionIdx: Int,
    initCityIdx: Int,
    onConfirm: (regionIdx: Int, cityIdx: Int, city: City) -> Unit,
    onDismiss: () -> Unit
) {
    var tmpRegionIdx by remember { mutableIntStateOf(initRegionIdx.coerceIn(0, Cities.REGIONS.lastIndex)) }
    var tmpCityIdx by remember { mutableIntStateOf(
        initCityIdx.coerceIn(0, Cities.REGIONS[tmpRegionIdx].cities.lastIndex)
    ) }
    val sheetState = rememberModalBottomSheetState(skipPartiallyExpanded = true)

    val regionListState = rememberLazyListState()
    val cityListState = rememberLazyListState()

    // 首次显示滚动到当前位置
    LaunchedEffect(Unit) {
        if (tmpRegionIdx > 0) regionListState.scrollToItem(tmpRegionIdx)
        if (tmpCityIdx > 0) cityListState.scrollToItem(tmpCityIdx)
    }
    // 切换省份时右栏回滚到顶
    LaunchedEffect(tmpRegionIdx) { cityListState.scrollToItem(0) }

    ModalBottomSheet(
        onDismissRequest = onDismiss,
        sheetState = sheetState,
        containerColor = Surface,
        dragHandle = null
    ) {
        Column(Modifier.fillMaxWidth()) {
            // 标题栏
            Row(
                Modifier.fillMaxWidth(),
                verticalAlignment = Alignment.CenterVertically
            ) {
                TextButton(onClick = onDismiss) {
                    Text("取消", color = TextSecondary, fontSize = Dimens.fontBody)
                }
                Text(
                    "选择出生地",
                    modifier = Modifier.weight(1f),
                    textAlign = TextAlign.Center,
                    color = OnSurface,
                    fontSize = Dimens.fontSubtitle,
                    fontWeight = FontWeight.Bold
                )
                TextButton(onClick = {
                    val r = Cities.REGIONS[tmpRegionIdx]
                    val c = r.cities[tmpCityIdx]
                    onConfirm(tmpRegionIdx, tmpCityIdx, c)
                }) {
                    Text("确定", color = Primary, fontSize = Dimens.fontBody, fontWeight = FontWeight.Bold)
                }
            }
            HorizontalDivider(color = DividerColor, thickness = 0.5.dp)

            // 当前选中预览
            Text(
                summary(tmpRegionIdx, tmpCityIdx),
                color = TextSecondary,
                fontSize = Dimens.fontSmall,
                modifier = Modifier
                    .fillMaxWidth()
                    .padding(horizontal = 16.dp, vertical = 8.dp),
                maxLines = 1
            )

            // 左省 + 右市
            Row(Modifier.fillMaxWidth().height(360.dp)) {
                // 左栏: 省/直辖市
                LazyColumn(
                    state = regionListState,
                    modifier = Modifier
                        .weight(1f)
                        .fillMaxHeight()
                        .background(Surface)
                ) {
                    itemsIndexed(
                        Cities.REGIONS,
                        key = { idx, r -> "r_${idx}_${r.name}" }
                    ) { idx, r ->
                        val active = idx == tmpRegionIdx
                        Row(
                            Modifier
                                .fillMaxWidth()
                                .height(44.dp)
                                .background(if (active) Background else Surface)
                                .clickable {
                                    if (idx != tmpRegionIdx) {
                                        tmpRegionIdx = idx
                                        tmpCityIdx = 0
                                    }
                                }
                                .padding(start = 14.dp, end = 8.dp),
                            verticalAlignment = Alignment.CenterVertically
                        ) {
                            Text(
                                r.name,
                                fontSize = Dimens.fontCaption,
                                color = if (active) Primary else OnSurface,
                                fontWeight = if (active) FontWeight.Bold else FontWeight.Normal
                            )
                        }
                    }
                }

                // 右栏: 区/县 (随左栏切换)
                LazyColumn(
                    state = cityListState,
                    modifier = Modifier
                        .weight(2f)
                        .fillMaxHeight()
                        .background(Background)
                ) {
                    itemsIndexed(
                        Cities.REGIONS[tmpRegionIdx].cities,
                        key = { idx, c -> "c_${tmpRegionIdx}_${idx}_${c.name}" }
                    ) { idx, c ->
                        val active = idx == tmpCityIdx
                        Row(
                            Modifier
                                .fillMaxWidth()
                                .height(44.dp)
                                .background(if (active) Surface else Background)
                                .clickable { tmpCityIdx = idx }
                                .padding(start = 14.dp, end = 8.dp),
                            verticalAlignment = Alignment.CenterVertically
                        ) {
                            Text(
                                c.name,
                                fontSize = Dimens.fontCaption,
                                color = if (active) Primary else OnSurface,
                                fontWeight = if (active) FontWeight.Bold else FontWeight.Normal
                            )
                        }
                    }
                }
            }

            Spacer(Modifier.height(16.dp))
        }
    }
}

private fun summary(regionIdx: Int, cityIdx: Int): String {
    val r = Cities.REGIONS[regionIdx]
    val c = r.cities[cityIdx]
    val ew = if (c.lon >= 0) "E" else "W"
    val ns = if (c.lat >= 0) "N" else "S"
    val lon = "%.2f".format(abs(c.lon))
    val lat = "%.2f".format(abs(c.lat))
    return "${r.name} · ${c.name} · $lon°$ew, $lat°$ns"
}
