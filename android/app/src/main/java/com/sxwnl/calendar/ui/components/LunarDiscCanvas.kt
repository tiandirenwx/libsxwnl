package com.sxwnl.calendar.ui.components

import androidx.compose.foundation.Canvas
import androidx.compose.foundation.background
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.shape.RoundedCornerShape
import androidx.compose.material3.Text
import androidx.compose.runtime.*
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.draw.clip
import androidx.compose.ui.geometry.Offset
import androidx.compose.ui.graphics.Brush
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.graphics.drawscope.DrawScope
import androidx.compose.ui.graphics.drawscope.Stroke
import androidx.compose.ui.text.font.FontFamily
import androidx.compose.ui.text.font.FontWeight
import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.sp
import com.sxwnl.calendar.data.LunarEclipseFrame
import com.sxwnl.calendar.ui.theme.*
import com.sxwnl.calendar.util.EclipseUtil
import kotlin.math.sqrt

/**
 * 月食过程动画:
 *   - 中央: 地球本影圆 + 半影圆 (固定)
 *   - 月球: 按帧 moveTo(moonX, moonY); coverage > 0.6 时染血月色
 */
@Composable
fun LunarDiscCanvas(
    frames: Array<LunarEclipseFrame>,
    currentJd: Double,
    tzHours: Double = 8.0,
    modifier: Modifier = Modifier
) {
    if (frames.isEmpty()) return

    val frameJds = remember(frames) { DoubleArray(frames.size) { frames[it].jd } }
    val (idx, t) = EclipseUtil.frameLerp(frameJds, currentJd)
    val idx2 = (idx + 1).coerceAtMost(frames.size - 1)
    val f1 = frames[idx]; val f2 = frames[idx2]

    val moonX = EclipseUtil.lerp(f1.moonX, f2.moonX, t)
    val moonY = EclipseUtil.lerp(f1.moonY, f2.moonY, t)
    val moonR = EclipseUtil.lerp(f1.moonRadius, f2.moonRadius, t)
    val umbra = EclipseUtil.lerp(f1.umbraRadius, f2.umbraRadius, t)
    val penum = EclipseUtil.lerp(f1.penumbraRadius, f2.penumbraRadius, t)
    val coverage = EclipseUtil.lerp(f1.coverage, f2.coverage, t)

    Box(
        modifier.fillMaxWidth().aspectRatio(1f)
            .clip(RoundedCornerShape(16.dp))
            .background(Brush.radialGradient(listOf(SkyMidNight, SkyDeepNight)))
    ) {
        Canvas(Modifier.fillMaxSize()) {
            val cx = size.width / 2
            val cy = size.height / 2

            // 把半影圆映射到屏幕 38% 半径
            val penumPx = size.minDimension * 0.38f
            val scale = (penumPx / penum).toFloat()
            val umbraPx = (umbra * scale).toFloat()
            val moonPxR = (moonR * scale).toFloat()
            val mx = cx + (moonX * scale).toFloat()
            val my = cy + (moonY * scale).toFloat()

            drawStars(seed = 99)

            // 半影 (径向渐变)
            drawCircle(
                brush = Brush.radialGradient(
                    listOf(EarthPenumbra.copy(alpha = 0.55f),
                           EarthPenumbra.copy(alpha = 0.05f)),
                    center = Offset(cx, cy), radius = penumPx
                ),
                radius = penumPx, center = Offset(cx, cy)
            )
            // 半影边界
            drawCircle(EarthPenumbra.copy(alpha = 0.8f),
                radius = penumPx, center = Offset(cx, cy),
                style = Stroke(width = 1f))

            // 本影 (实心暗)
            drawCircle(
                brush = Brush.radialGradient(
                    listOf(EarthUmbra, EarthUmbra.copy(alpha = 0.8f)),
                    center = Offset(cx, cy), radius = umbraPx
                ),
                radius = umbraPx, center = Offset(cx, cy)
            )
            // 本影边界
            drawCircle(Color(0xFF000000),
                radius = umbraPx, center = Offset(cx, cy),
                style = Stroke(width = 1f))

            // 月球: coverage 越高越偏血月; 全食时 (coverage>1) 完全血色
            val moonColor = when {
                coverage <= 0 -> MoonBright
                coverage >= 1 -> MoonBlood
                else -> {
                    val c = coverage.toFloat()
                    Color(
                        red = MoonBright.red * (1 - c) + MoonBlood.red * c,
                        green = MoonBright.green * (1 - c) + MoonBlood.green * c,
                        blue = MoonBright.blue * (1 - c) + MoonBlood.blue * c
                    )
                }
            }
            drawCircle(moonColor, radius = moonPxR, center = Offset(mx, my))
            // 月球缘
            drawCircle(Color(0x55000000),
                radius = moonPxR, center = Offset(mx, my),
                style = Stroke(width = 1.2f))
        }

        Column(Modifier.align(Alignment.TopStart).padding(10.dp)) {
            DataChip("食分", "%.3f".format(coverage.coerceAtLeast(0.0)))
        }
        Column(Modifier.align(Alignment.TopEnd).padding(10.dp)) {
            DataChip("时刻", EclipseUtil.jdTdToLocal(currentJd, tzHours, withDate = false))
        }
        // 中央标签: 地影
        Text("地球本影/半影", color = Color.White.copy(alpha = 0.35f),
            fontSize = 10.sp,
            modifier = Modifier.align(Alignment.BottomCenter).padding(bottom = 8.dp))
    }
}

private fun DrawScope.drawStars(seed: Int) {
    var s = seed
    for (i in 0..50) {
        s = (s * 1103515245 + 12345) and 0x7fffffff
        val x = (s % 1000) / 1000f * size.width
        s = (s * 1103515245 + 12345) and 0x7fffffff
        val y = (s % 1000) / 1000f * size.height
        s = (s * 1103515245 + 12345) and 0x7fffffff
        val a = 0.20f + (s % 100) / 220f
        drawCircle(Color.White.copy(alpha = a), radius = 1.0f, center = Offset(x, y))
    }
}

@Composable
private fun DataChip(label: String, value: String) {
    Row(
        Modifier
            .clip(RoundedCornerShape(50))
            .background(Color(0xCC000000))
            .padding(horizontal = 10.dp, vertical = 4.dp),
        verticalAlignment = Alignment.CenterVertically
    ) {
        Text(label, color = Color.White.copy(alpha = 0.7f), fontSize = 10.sp)
        Spacer(Modifier.width(6.dp))
        Text(value, color = Color.White, fontSize = 12.sp,
            fontWeight = FontWeight.Bold, fontFamily = FontFamily.Monospace)
    }
}
