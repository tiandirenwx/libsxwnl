package com.sxwnl.calendar.ui.components

import androidx.compose.foundation.Canvas
import androidx.compose.foundation.background
import androidx.compose.foundation.gestures.detectDragGestures
import androidx.compose.foundation.gestures.detectTapGestures
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.shape.RoundedCornerShape
import androidx.compose.material3.Text
import androidx.compose.runtime.*
import androidx.compose.ui.Modifier
import androidx.compose.ui.draw.clip
import androidx.compose.ui.geometry.Offset
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.graphics.drawscope.Stroke
import androidx.compose.ui.input.pointer.pointerInput
import androidx.compose.ui.platform.LocalDensity
import androidx.compose.ui.text.font.FontFamily
import androidx.compose.ui.text.font.FontWeight
import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.sp
import com.sxwnl.calendar.ui.theme.*
import com.sxwnl.calendar.util.EclipseUtil

/** 一个时间轴节点 */
data class TimelineMark(val label: String, val jd: Double, val color: Color)

/**
 * 横向时间轴; 用于日食 5 点(C1 / C2 / Max / C3 / C4) 或月食 7 点(P1..P4).
 * 支持点击节点跳转 + 拖动 scrubber.
 */
@Composable
fun EclipseTimelineBar(
    marks: List<TimelineMark>,
    currentJd: Double,
    onSeek: (Double) -> Unit,
    tzHours: Double = 8.0,
    modifier: Modifier = Modifier
) {
    val filtered = marks.filter { it.jd > 0 }.sortedBy { it.jd }
    if (filtered.isEmpty()) return

    val minJd = filtered.first().jd
    val maxJd = filtered.last().jd
    val span = (maxJd - minJd).coerceAtLeast(1.0 / 1440.0) // 至少 1 分钟

    val density = LocalDensity.current
    var widthPx by remember { mutableStateOf(1f) }

    fun jdToX(jd: Double): Float =
        ((jd - minJd) / span).coerceIn(0.0, 1.0).toFloat() * widthPx
    fun xToJd(x: Float): Double =
        (minJd + (x / widthPx).coerceIn(0f, 1f) * span)

    Column(modifier.fillMaxWidth()) {
        // 当前时间显示
        Text(
            EclipseUtil.jdTdToLocal(currentJd, tzHours, withDate = false) +
                " " + EclipseUtil.tzLabel(tzHours),
            fontSize = 13.sp, fontFamily = FontFamily.Monospace,
            color = OnSurface, fontWeight = FontWeight.Bold,
            modifier = Modifier.padding(start = 4.dp, bottom = 4.dp)
        )

        // 画布区
        Box(
            Modifier.fillMaxWidth().height(56.dp)
                .clip(RoundedCornerShape(8.dp))
                .background(SkyMidNight)
                .pointerInput(filtered, span) {
                    detectTapGestures { off -> onSeek(xToJd(off.x)) }
                }
                .pointerInput(filtered, span) {
                    detectDragGestures { change, _ -> onSeek(xToJd(change.position.x)) }
                }
        ) {
            Canvas(Modifier.fillMaxSize()) {
                widthPx = size.width
                val midY = size.height / 2

                drawLine(
                    color = Color(0x33FFFFFF),
                    start = Offset(8f, midY), end = Offset(size.width - 8f, midY),
                    strokeWidth = 2f
                )

                filtered.forEach { m ->
                    val x = jdToX(m.jd)
                    drawCircle(m.color, radius = 5f, center = Offset(x, midY))
                    drawCircle(Color.White, radius = 5f, center = Offset(x, midY),
                        style = Stroke(width = 1.5f))
                }

                val cx = jdToX(currentJd)
                drawLine(
                    color = PathMarker,
                    start = Offset(cx, 4f), end = Offset(cx, size.height - 4f),
                    strokeWidth = 2.5f
                )
                drawCircle(PathMarker, radius = 7f, center = Offset(cx, midY))
                drawCircle(Color.Black, radius = 7f, center = Offset(cx, midY),
                    style = Stroke(width = 1.5f))
            }
        }

        // 节点标签条
        BoxWithConstraints(Modifier.fillMaxWidth().height(28.dp).padding(top = 2.dp)) {
            val widthDp = maxWidth
            val widthPxLabel = with(density) { widthDp.toPx() }
            filtered.forEach { m ->
                val x = ((m.jd - minJd) / span).toFloat().coerceIn(0f, 1f) * widthPxLabel
                val labelOffsetDp = with(density) { x.toDp() } - 22.dp
                Column(
                    Modifier.absoluteOffset(x = labelOffsetDp, y = 0.dp).width(44.dp)
                ) {
                    Text(
                        m.label, fontSize = 10.sp, color = m.color,
                        fontWeight = FontWeight.Bold,
                        modifier = Modifier.fillMaxWidth(),
                        textAlign = androidx.compose.ui.text.style.TextAlign.Center
                    )
                    Text(
                        EclipseUtil.jdTdToLocal(m.jd, tzHours, withDate = false).substring(0, 5),
                        fontSize = 9.sp, color = TextTertiary,
                        fontFamily = FontFamily.Monospace,
                        modifier = Modifier.fillMaxWidth(),
                        textAlign = androidx.compose.ui.text.style.TextAlign.Center
                    )
                }
            }
        }
    }
}
