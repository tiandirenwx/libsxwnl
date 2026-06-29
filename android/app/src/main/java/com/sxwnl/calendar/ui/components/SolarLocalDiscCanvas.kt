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
import com.sxwnl.calendar.data.EclipseFrame
import com.sxwnl.calendar.ui.theme.*
import com.sxwnl.calendar.util.EclipseUtil

/**
 * 本地日食观测: 渲染当前 jd 对应帧的太阳/月亮圆盘.
 *
 * 半径单位为弧度. 屏幕缩放: 太阳半径占画布 32%, 月亮按比例.
 * 月亮覆盖太阳形成"咬月"效果.
 */
@Composable
fun SolarLocalDiscCanvas(
    frames: Array<EclipseFrame>,
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
    val sunR  = EclipseUtil.lerp(f1.sunRadius, f2.sunRadius, t)
    val mag   = EclipseUtil.lerp(f1.magnitude, f2.magnitude, t)

    Box(
        modifier.fillMaxWidth().aspectRatio(1f)
            .clip(RoundedCornerShape(16.dp))
            .background(Brush.radialGradient(listOf(SkyMidNight, SkyDeepNight)))
    ) {
        Canvas(Modifier.fillMaxSize()) {
            val cx = size.width / 2
            val cy = size.height / 2
            val sunPx = size.minDimension * 0.32f
            val scale = (sunPx / sunR).toFloat()
            val mPx = (moonR * scale).toFloat()
            val mxPx = (moonX * scale).toFloat()
            val myPx = (moonY * scale).toFloat()

            drawStars(seed = 42)

            // 日冕光晕 (食甚附近显著)
            val haloAlpha = (0.18f + 0.32f * mag.toFloat()).coerceIn(0.18f, 0.55f)
            drawCircle(SunHalo.copy(alpha = haloAlpha),
                radius = sunPx * 1.7f, center = Offset(cx, cy))
            drawCircle(SunHalo.copy(alpha = haloAlpha * 0.6f),
                radius = sunPx * 1.3f, center = Offset(cx, cy))

            // 太阳本体
            drawCircle(
                brush = Brush.radialGradient(
                    listOf(SunCore, SunGlow),
                    center = Offset(cx, cy),
                    radius = sunPx
                ),
                radius = sunPx, center = Offset(cx, cy)
            )

            // 月亮 (黑盘 + 边缘)
            val moonCx = cx + mxPx
            val moonCy = cy + myPx
            drawCircle(MoonDark, radius = mPx, center = Offset(moonCx, moonCy))
            drawCircle(MoonRim, radius = mPx,
                center = Offset(moonCx, moonCy),
                style = Stroke(width = 1.5f))
        }

        Column(Modifier.align(Alignment.TopStart).padding(10.dp)) {
            DataChip("食分", "%.3f".format(mag.coerceAtLeast(0.0)))
        }
        Column(Modifier.align(Alignment.TopEnd).padding(10.dp)) {
            DataChip("时刻", EclipseUtil.jdTdToLocal(currentJd, tzHours, withDate = false))
        }
    }
}

private fun DrawScope.drawStars(seed: Int) {
    var s = seed
    for (i in 0..40) {
        s = (s * 1103515245 + 12345) and 0x7fffffff
        val x = (s % 1000) / 1000f * size.width
        s = (s * 1103515245 + 12345) and 0x7fffffff
        val y = (s % 1000) / 1000f * size.height
        s = (s * 1103515245 + 12345) and 0x7fffffff
        val a = 0.18f + (s % 100) / 250f
        drawCircle(Color.White.copy(alpha = a), radius = 1.1f, center = Offset(x, y))
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
