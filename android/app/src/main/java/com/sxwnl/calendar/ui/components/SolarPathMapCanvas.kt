package com.sxwnl.calendar.ui.components

import androidx.compose.foundation.Canvas
import androidx.compose.foundation.background
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.shape.RoundedCornerShape
import androidx.compose.runtime.*
import androidx.compose.ui.Modifier
import androidx.compose.ui.draw.clip
import androidx.compose.ui.geometry.Offset
import androidx.compose.ui.geometry.Size
import androidx.compose.ui.graphics.Brush
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.graphics.Path
import androidx.compose.ui.graphics.PathEffect
import androidx.compose.ui.graphics.drawscope.DrawScope
import androidx.compose.ui.graphics.drawscope.Stroke
import androidx.compose.ui.graphics.drawscope.translate
import androidx.compose.ui.unit.dp
import com.sxwnl.calendar.data.EclipseGeoPoint
import com.sxwnl.calendar.data.SolarEclipsePath
import com.sxwnl.calendar.ui.theme.*
import kotlin.math.abs

/**
 * 世界地图 + 日食食带 + 动画影点.
 *
 * 投影: 等距圆柱 (Plate Carrée), 经度 [-π,π] → x [0,w], 纬度 [π/2,-π/2] → y [0,h].
 * - ditu0/1/2 都是弧度坐标, 段间 1e7 标记 Move-To.
 * - centerLine / umbra* / penumbra* 都是度数 (来自 capi 已乘 radd).
 *
 * @param worldMapRad 海岸线弧度 (ditu0 或 ditu1) lon/lat 交替, 段间 1e7
 * @param bordersRad  可选: ditu2 国界线 (弧度), null 时不绘制
 * @param path 日食路径数据
 * @param currentJd 当前播放时间; 影点位置据此沿 centerLine 取
 */
@Composable
fun SolarPathMapCanvas(
    worldMapRad: DoubleArray,
    bordersRad: DoubleArray?,
    path: SolarEclipsePath,
    currentJd: Double,
    modifier: Modifier = Modifier
) {
    Box(
        modifier.fillMaxWidth().aspectRatio(2f)
            .clip(RoundedCornerShape(12.dp))
            .background(MapOcean)
    ) {
        Canvas(Modifier.fillMaxSize()) {
            drawRect(
                brush = Brush.verticalGradient(
                    listOf(SkyDeepNight, MapOcean, SkyDeepNight)
                ),
                size = size
            )
            drawGrid()
            drawCoastline(worldMapRad, MapLand, 1.0f)
            if (bordersRad != null && bordersRad.isNotEmpty()) {
                drawCoastline(bordersRad, MapBorder, 0.7f)
            }
            drawPathBand(path)
            drawCenterLine(path.centerLine)
            drawKeyPoints(path)
            drawMovingShadow(path.centerLine, currentJd)
            drawAxisLabels()
        }
    }
}

// ── 绘制工具 ──────────────────────────────────────────────────

private fun DrawScope.lonLatToXy(lonDeg: Double, latDeg: Double): Offset {
    // lon 范围 [-180, 180] → x [0, w]
    var ln = ((lonDeg + 540) % 360) - 180
    val x = ((ln + 180) / 360.0).toFloat() * size.width
    val y = ((90 - latDeg.coerceIn(-90.0, 90.0)) / 180.0).toFloat() * size.height
    return Offset(x, y)
}

private fun DrawScope.lonLatRadToXy(lonRad: Double, latRad: Double): Offset =
    lonLatToXy(Math.toDegrees(lonRad), Math.toDegrees(latRad))

private fun DrawScope.drawGrid() {
    val w = size.width; val h = size.height
    // 经线
    for (i in 0..12) {
        val x = i / 12f * w
        drawLine(MapGrid, Offset(x, 0f), Offset(x, h), strokeWidth = 0.6f)
    }
    // 纬线
    for (i in 0..6) {
        val y = i / 6f * h
        drawLine(MapGrid, Offset(0f, y), Offset(w, y), strokeWidth = 0.6f)
    }
    // 赤道、本初子午线加粗一点
    drawLine(Color(0x55FFFFFF), Offset(0f, h / 2), Offset(w, h / 2), strokeWidth = 1.2f)
    drawLine(Color(0x55FFFFFF), Offset(w / 2, 0f), Offset(w / 2, h), strokeWidth = 1.2f)
}

private fun DrawScope.drawCoastline(
    data: DoubleArray, color: Color, strokeWidth: Float
) {
    if (data.isEmpty()) return
    val path = Path()
    var moveNext = true
    var i = 0
    while (i < data.size) {
        if (data[i] >= 9e6) { // 1e7 segment break
            moveNext = true
            i++
            continue
        }
        if (i + 1 >= data.size) break
        val lonRad = data[i]
        val latRad = data[i + 1]
        val p = lonLatRadToXy(lonRad, latRad)
        if (moveNext) {
            path.moveTo(p.x, p.y)
            moveNext = false
        } else {
            path.lineTo(p.x, p.y)
        }
        i += 2
    }
    drawPath(path, color = color, style = Stroke(width = strokeWidth))
}

private fun DrawScope.drawPolyline(
    pts: Array<EclipseGeoPoint>, color: Color, strokeWidth: Float,
    dashed: Boolean = false
) {
    if (pts.size < 2) return
    val path = Path()
    var first = true
    var prevX = 0f
    for (p in pts) {
        val off = lonLatToXy(p.longitude, p.latitude)
        // 防止跨日期变更线时的水平大跳: 若与上一个 X 差 > 半屏, 视为新段
        if (first || abs(off.x - prevX) > size.width / 2f) {
            path.moveTo(off.x, off.y)
            first = false
        } else {
            path.lineTo(off.x, off.y)
        }
        prevX = off.x
    }
    val effect = if (dashed)
        PathEffect.dashPathEffect(floatArrayOf(8f, 5f), 0f) else null
    drawPath(path, color = color,
        style = Stroke(width = strokeWidth, pathEffect = effect))
}

private fun DrawScope.drawPathBand(path: SolarEclipsePath) {
    // 半影 (虚线, 较浅)
    drawPolyline(path.penumbraNorth, PathPenumbra.copy(alpha = 0.7f), 1.5f, dashed = true)
    drawPolyline(path.penumbraSouth, PathPenumbra.copy(alpha = 0.7f), 1.5f, dashed = true)
    // 本影 (实线)
    drawPolyline(path.umbraNorth, PathUmbra, 2f)
    drawPolyline(path.umbraSouth, PathUmbra, 2f)
}

private fun DrawScope.drawCenterLine(pts: Array<EclipseGeoPoint>) {
    drawPolyline(pts, PathCenter, 2.5f)
}

private fun DrawScope.drawKeyPoints(path: SolarEclipsePath) {
    val pts = listOf(
        path.partialStart to "P1",
        path.centralStart to "C1",
        path.greatestEclipse to "MAX",
        path.centralEnd to "C2",
        path.partialEnd to "P4"
    )
    for ((p, _) in pts) {
        if (p.longitude > 99 || p.latitude > 99) continue
        val off = lonLatToXy(p.longitude, p.latitude)
        drawCircle(Color.White, radius = 6f, center = off)
        drawCircle(PathCenter, radius = 4f, center = off)
    }
}

private fun DrawScope.drawMovingShadow(centerLine: Array<EclipseGeoPoint>, currentJd: Double) {
    if (centerLine.isEmpty()) return
    // 在中心线上找最近 jd 的两个点, 线性插值
    var idx = 0
    while (idx < centerLine.size - 1 && centerLine[idx + 1].jd < currentJd) idx++
    val p1 = centerLine[idx]
    val p2 = if (idx + 1 < centerLine.size) centerLine[idx + 1] else p1
    val span = (p2.jd - p1.jd).takeIf { it > 0 } ?: 1.0
    val t = ((currentJd - p1.jd) / span).coerceIn(0.0, 1.0).toFloat()
    val lon = p1.longitude + (p2.longitude - p1.longitude) * t
    val lat = p1.latitude + (p2.latitude - p1.latitude) * t
    val off = lonLatToXy(lon, lat)

    drawCircle(SunHalo, radius = 28f, center = off)
    drawCircle(PathMarker.copy(alpha = 0.5f), radius = 16f, center = off)
    drawCircle(MoonDark, radius = 8f, center = off)
    drawCircle(PathMarker, radius = 8f, center = off, style = Stroke(width = 2f))
}

private fun DrawScope.drawAxisLabels() {
    // 网格本身够直观, 此处保留扩展点 (后续可叠加文字标注 90W/0/90E).
}
