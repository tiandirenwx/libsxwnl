package com.sxwnl.calendar.ui.components

import androidx.compose.foundation.background
import androidx.compose.foundation.clickable
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.shape.RoundedCornerShape
import androidx.compose.material3.Text
import androidx.compose.runtime.*
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.draw.clip
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.text.font.FontWeight
import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.sp
import com.sxwnl.calendar.ui.theme.*
import kotlinx.coroutines.delay
import kotlinx.coroutines.isActive

/**
 * 动画播放控制器: 播放/暂停 + 速度选择.
 * 调用方提供 jdRange = (begin, end); 控制器输出当前 currentJd.
 *
 * 速度: 1× ≈ 实时, 60× ≈ 1 秒画完 1 分钟, 600× ≈ 10 分钟/秒.
 * 默认 600× 让常见 3 小时日食在 18 秒内播完一遍.
 */
@Composable
fun PlaybackController(
    jdBegin: Double,
    jdEnd: Double,
    currentJd: Double,
    onChange: (Double) -> Unit,
    modifier: Modifier = Modifier
) {
    var playing by remember { mutableStateOf(false) }
    var speedIdx by remember { mutableStateOf(1) } // 0=0.5×, 1=1×, 2=2×, 3=5×

    val speedMultipliers = listOf(0.5, 1.0, 2.0, 5.0)
    val speedLabels = listOf("0.5×", "1×", "2×", "5×")

    // 30 fps; 1× ≈ 0.0001 jd / 帧 ≈ 8.64 实际秒/动画秒
    LaunchedEffect(playing, speedIdx, jdBegin, jdEnd) {
        if (!playing) return@LaunchedEffect
        val totalSpan = jdEnd - jdBegin
        if (totalSpan <= 0) return@LaunchedEffect
        var jd = if (currentJd in jdBegin..jdEnd) currentJd else jdBegin
        while (isActive && playing) {
            delay(33L)
            val step = 0.0001 * speedMultipliers[speedIdx]
            jd += step
            if (jd >= jdEnd) {
                jd = jdEnd
                playing = false
            }
            onChange(jd)
        }
    }

    Row(
        modifier.fillMaxWidth(),
        verticalAlignment = Alignment.CenterVertically,
        horizontalArrangement = Arrangement.spacedBy(8.dp)
    ) {
        // 播放/暂停按钮
        Box(
            Modifier.size(40.dp)
                .clip(RoundedCornerShape(50))
                .background(if (playing) Accent else Primary)
                .clickable {
                    if (!playing && currentJd >= jdEnd - 1e-6) {
                        onChange(jdBegin)
                    }
                    playing = !playing
                },
            contentAlignment = Alignment.Center
        ) {
            Text(if (playing) "❚❚" else "▶",
                color = Color.White, fontSize = 14.sp, fontWeight = FontWeight.Bold)
        }

        // 重置
        Box(
            Modifier.size(40.dp)
                .clip(RoundedCornerShape(50))
                .background(SkyMidNight)
                .clickable {
                    playing = false
                    onChange(jdBegin)
                },
            contentAlignment = Alignment.Center
        ) {
            Text("⟲", color = Color.White, fontSize = 16.sp)
        }

        // 速度选择
        Row(
            Modifier.weight(1f)
                .clip(RoundedCornerShape(50))
                .background(SkyMidNight),
            verticalAlignment = Alignment.CenterVertically
        ) {
            speedLabels.forEachIndexed { i, label ->
                val sel = speedIdx == i
                Box(
                    Modifier.weight(1f).height(36.dp)
                        .clip(RoundedCornerShape(50))
                        .background(if (sel) Primary else Color.Transparent)
                        .clickable { speedIdx = i },
                    contentAlignment = Alignment.Center
                ) {
                    Text(label,
                        color = if (sel) Color.White else Color.White.copy(alpha = 0.7f),
                        fontSize = 11.sp,
                        fontWeight = if (sel) FontWeight.Bold else FontWeight.Normal)
                }
            }
        }
    }
}
