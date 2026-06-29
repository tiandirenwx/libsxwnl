package com.sxwnl.calendar.ui.components

import androidx.compose.foundation.BorderStroke
import androidx.compose.foundation.background
import androidx.compose.foundation.border
import androidx.compose.foundation.clickable
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.lazy.LazyColumn
import androidx.compose.foundation.lazy.items
import androidx.compose.foundation.shape.RoundedCornerShape
import androidx.compose.foundation.text.BasicTextField
import androidx.compose.foundation.text.KeyboardOptions
import androidx.compose.material3.Text
import androidx.compose.runtime.*
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.draw.clip
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.graphics.SolidColor
import androidx.compose.ui.text.TextStyle
import androidx.compose.ui.text.font.FontFamily
import androidx.compose.ui.text.font.FontWeight
import androidx.compose.ui.text.input.KeyboardType
import androidx.compose.ui.text.style.TextAlign
import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.sp
import com.sxwnl.calendar.data.Cities
import com.sxwnl.calendar.data.City
import com.sxwnl.calendar.ui.theme.*
import com.sxwnl.calendar.util.EclipseUtil

/**
 * 紧凑的观测点选择器: 横向滚动城市按钮 + "自定义" 弹出经纬度输入.
 */
@Composable
fun LocationStrip(
    selected: City,
    onSelect: (City) -> Unit,
    modifier: Modifier = Modifier
) {
    var showCustom by remember { mutableStateOf(false) }

    Column(modifier.fillMaxWidth()) {
        Row(verticalAlignment = Alignment.CenterVertically,
            horizontalArrangement = Arrangement.spacedBy(6.dp)) {
            Text("📍 观测地", fontSize = 12.sp, color = TextSecondary)
            Spacer(Modifier.weight(1f))
            Text("${selected.name}  ${EclipseUtil.lonLabel(selected.lon)}  " +
                 EclipseUtil.latLabel(selected.lat) + "  ${EclipseUtil.tzLabel(selected.tz)}",
                fontSize = 11.sp, color = OnSurface, fontFamily = FontFamily.Monospace)
        }
        Spacer(Modifier.height(6.dp))
        LazyColumn(Modifier.heightIn(min = 38.dp, max = 38.dp)) {
            item {
                Row(horizontalArrangement = Arrangement.spacedBy(6.dp)) {
                    Cities.PRESET.forEach { c ->
                        CityChip(c, c.name == selected.name) { onSelect(c) }
                    }
                    Box(
                        Modifier.height(34.dp)
                            .clip(RoundedCornerShape(50))
                            .background(SkyDawn)
                            .clickable { showCustom = true }
                            .padding(horizontal = 12.dp),
                        contentAlignment = Alignment.Center
                    ) {
                        Text("✎ 自定义", color = Color.White, fontSize = 12.sp,
                            fontWeight = FontWeight.Medium)
                    }
                }
            }
        }
    }

    if (showCustom) {
        CustomLocationDialog(
            initial = selected,
            onDismiss = { showCustom = false },
            onConfirm = { c ->
                onSelect(c)
                showCustom = false
            }
        )
    }
}

@Composable
private fun CityChip(city: City, selected: Boolean, onClick: () -> Unit) {
    val bg = if (selected) Primary else Surface
    val fg = if (selected) Color.White else OnSurface
    Box(
        Modifier.height(34.dp)
            .clip(RoundedCornerShape(50))
            .background(bg)
            .border(BorderStroke(1.dp, if (selected) Primary else DividerColor),
                RoundedCornerShape(50))
            .clickable(onClick = onClick)
            .padding(horizontal = 12.dp),
        contentAlignment = Alignment.Center
    ) {
        Text(city.name, color = fg, fontSize = 12.sp,
            fontWeight = if (selected) FontWeight.Bold else FontWeight.Normal)
    }
}

@Composable
private fun CustomLocationDialog(
    initial: City,
    onDismiss: () -> Unit,
    onConfirm: (City) -> Unit
) {
    var lon by remember { mutableStateOf("%.4f".format(initial.lon)) }
    var lat by remember { mutableStateOf("%.4f".format(initial.lat)) }
    var tz  by remember { mutableStateOf("%.1f".format(initial.tz)) }

    androidx.compose.material3.AlertDialog(
        onDismissRequest = onDismiss,
        confirmButton = {
            androidx.compose.material3.TextButton(
                onClick = {
                    val L = lon.toDoubleOrNull() ?: initial.lon
                    val W = lat.toDoubleOrNull() ?: initial.lat
                    val Z = tz.toDoubleOrNull() ?: initial.tz
                    onConfirm(City("自定义", L, W, Z))
                }
            ) { Text("确定", color = Primary) }
        },
        dismissButton = {
            androidx.compose.material3.TextButton(onClick = onDismiss) {
                Text("取消")
            }
        },
        title = { Text("自定义观测地") },
        text = {
            Column(verticalArrangement = Arrangement.spacedBy(8.dp)) {
                NumRow("经度 (东+ 西-)", lon, { lon = it })
                NumRow("纬度 (北+ 南-)", lat, { lat = it })
                NumRow("时区 (东+, 北京=8)", tz, { tz = it })
            }
        }
    )
}

@Composable
private fun NumRow(label: String, value: String, onChange: (String) -> Unit) {
    Row(verticalAlignment = Alignment.CenterVertically,
        horizontalArrangement = Arrangement.spacedBy(8.dp)) {
        Text(label, fontSize = 12.sp, color = TextSecondary,
            modifier = Modifier.width(118.dp))
        Box(
            Modifier.weight(1f).height(38.dp)
                .clip(RoundedCornerShape(8.dp))
                .background(Background)
                .border(BorderStroke(1.dp, DividerColor), RoundedCornerShape(8.dp))
                .padding(horizontal = 8.dp),
            contentAlignment = Alignment.CenterStart
        ) {
            BasicTextField(
                value = value, onValueChange = onChange, singleLine = true,
                textStyle = TextStyle(color = OnSurface, fontSize = 14.sp,
                    textAlign = TextAlign.End),
                modifier = Modifier.fillMaxWidth(),
                keyboardOptions = KeyboardOptions(keyboardType = KeyboardType.Number),
                cursorBrush = SolidColor(Primary)
            )
        }
    }
}
