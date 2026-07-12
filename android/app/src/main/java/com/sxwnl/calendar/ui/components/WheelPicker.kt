package com.sxwnl.calendar.ui.components

import androidx.compose.foundation.border
import androidx.compose.foundation.gestures.snapping.rememberSnapFlingBehavior
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.lazy.LazyColumn
import androidx.compose.foundation.lazy.rememberLazyListState
import androidx.compose.runtime.*
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.platform.LocalDensity
import androidx.compose.ui.text.font.FontWeight
import androidx.compose.ui.text.style.TextAlign
import androidx.compose.ui.unit.dp
import androidx.compose.material3.Text
import com.sxwnl.calendar.ui.theme.OnSurface
import com.sxwnl.calendar.ui.theme.Primary
import com.sxwnl.calendar.ui.theme.TextSecondary
import androidx.compose.ui.unit.sp
import kotlinx.coroutines.flow.filter
import kotlinx.coroutines.launch

/**
 * 简易滚轮选择器: LazyColumn + snapFlingBehavior, 对齐鸿蒙 TextPicker 语义.
 *
 * 居中对齐的实现方式 — 前后各插 padCount 个空白占位行, **不用 contentPadding**:
 *
 *   列表结构 = [空白 × padCount] + [真实项 × size] + [空白 × padCount]
 *
 *   于是"居中项的数据索引" == "顶部可见项的列表位置(firstVisibleItemIndex)":
 *     - 视口共显示 visibleItemsCount 行, 正中那行是视口第 padCount 行;
 *     - 视口第 padCount 行的列表位置 = firstVisibleItemIndex + padCount;
 *     - 减去前置的 padCount 个空白 → 数据索引 = firstVisibleItemIndex.
 *
 *   这样 scrollToItem(dataIndex) 落位后, firstVisibleItemIndex == dataIndex,
 *   居中项恰好是 dataIndex, 与回调/高亮完全一致 —— 彻底避免 contentPadding 下
 *   "offset / scrollToItem 是否含顶部内边距"的版本歧义 (那正是之前选中错位、
 *   加粗行错行的根因).
 */
@OptIn(androidx.compose.foundation.ExperimentalFoundationApi::class)
@Composable
fun WheelPicker(
    items: List<String>,
    selectedIndex: Int,
    onSelectedChange: (Int) -> Unit,
    modifier: Modifier = Modifier,
    visibleItemsCount: Int = 5,
    itemHeight: androidx.compose.ui.unit.Dp = 36.dp
) {
    val options = items
    if (options.isEmpty()) {
        Box(modifier.height(itemHeight * visibleItemsCount))
        return
    }
    val padCount = visibleItemsCount / 2
    val totalCount = options.size + padCount * 2
    val safeIndex = selectedIndex.coerceIn(0, options.size - 1)

    // 顶部可见项列表位置 == 居中数据索引, 故初始 firstVisibleItemIndex = safeIndex.
    val listState = rememberLazyListState(initialFirstVisibleItemIndex = safeIndex)
    val flingBehavior = rememberSnapFlingBehavior(lazyListState = listState)
    val scope = rememberCoroutineScope()
    val halfItemPx = with(LocalDensity.current) { itemHeight.toPx() } / 2f

    // 居中数据索引 = firstVisibleItemIndex + (滚过半行则进位). 停顿时 offset≈0.
    val centeredIndex by remember(options.size) {
        derivedStateOf {
            val first = listState.firstVisibleItemIndex
            val off = listState.firstVisibleItemScrollOffset
            val adj = first + if (off > halfItemPx) 1 else 0
            adj.coerceIn(0, options.size - 1)
        }
    }

    // rememberUpdatedState: 长生命周期的 LaunchedEffect(listState) 始终读最新值.
    val currentSelectedIndex by rememberUpdatedState(selectedIndex)
    val currentOnSelectedChange by rememberUpdatedState(onSelectedChange)

    // 滚动停顿 → 回调"居中项"(用户真正选到的项).
    //   注意: 不能在 filter{!it} 后再 distinctUntilChanged —— snapshotFlow 只在
    //   变化时发射, 序列为 false(初始)→true(滚动)→false(停止), filter 掉 true 后
    //   剩两个相邻 false, distinctUntilChanged 会把"松手后的 false"当重复去掉,
    //   导致停顿回调永不触发, 确定时仍是初始值 (1月1日12:00). snapshotFlow 本身
    //   保证两次 false 之间必有 true, 故无需去重.
    LaunchedEffect(listState) {
        snapshotFlow { listState.isScrollInProgress }
            .filter { !it }
            .collect {
                if (centeredIndex != currentSelectedIndex) {
                    currentOnSelectedChange(centeredIndex)
                }
            }
    }

    // 外部 selectedIndex 变化(如年/月切换致日号重算) → 滚到该项居中.
    LaunchedEffect(selectedIndex, options.size) {
        val target = selectedIndex.coerceIn(0, options.size - 1)
        if (centeredIndex != target) {
            scope.launch { listState.scrollToItem(target) }
        }
    }

    Box(modifier.height(itemHeight * visibleItemsCount)) {
        LazyColumn(
            state = listState,
            flingBehavior = flingBehavior,
            modifier = Modifier.fillMaxSize(),
            horizontalAlignment = Alignment.CenterHorizontally
        ) {
            items(count = totalCount) { pos ->
                val dataIdx = pos - padCount
                Box(
                    Modifier.fillMaxWidth().height(itemHeight),
                    contentAlignment = Alignment.Center
                ) {
                    if (dataIdx in options.indices) {
                        val isSel = dataIdx == centeredIndex
                        Text(
                            options[dataIdx],
                            fontSize = if (isSel) 16.sp else 14.sp,
                            fontWeight = if (isSel) FontWeight.Bold else FontWeight.Normal,
                            color = if (isSel) OnSurface else TextSecondary,
                            textAlign = TextAlign.Center
                        )
                    }
                }
            }
        }

        // 中央高亮框 — 视口正中一行, 与居中项对齐
        Box(
            Modifier
                .align(Alignment.Center)
                .fillMaxWidth()
                .height(itemHeight)
                .border(
                    width = 0.5.dp,
                    color = Primary.copy(alpha = 0.35f)
                )
        )
    }
}
