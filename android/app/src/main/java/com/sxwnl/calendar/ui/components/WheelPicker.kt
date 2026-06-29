package com.sxwnl.calendar.ui.components

import androidx.compose.foundation.background
import androidx.compose.foundation.border
import androidx.compose.foundation.gestures.snapping.rememberSnapFlingBehavior
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.lazy.LazyColumn
import androidx.compose.foundation.lazy.rememberLazyListState
import androidx.compose.runtime.*
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.text.font.FontWeight
import androidx.compose.ui.text.style.TextAlign
import androidx.compose.ui.unit.dp
import androidx.compose.material3.Text
import com.sxwnl.calendar.ui.theme.OnSurface
import com.sxwnl.calendar.ui.theme.Primary
import com.sxwnl.calendar.ui.theme.TextSecondary
import androidx.compose.ui.unit.sp
import kotlinx.coroutines.flow.distinctUntilChanged
import kotlinx.coroutines.flow.filter
import kotlinx.coroutines.launch

/**
 * 简易滚轮选择器: LazyColumn + snapFlingBehavior.
 * 中央高亮项即为选中项, 滚动停止后通过 onSelectedChange 回调。
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
    val safeIndex = selectedIndex.coerceIn(0, options.size - 1)
    val listState = rememberLazyListState(initialFirstVisibleItemIndex = safeIndex)
    val flingBehavior = rememberSnapFlingBehavior(lazyListState = listState)
    val scope = rememberCoroutineScope()

    LaunchedEffect(listState) {
        snapshotFlow { listState.isScrollInProgress }
            .filter { !it }
            .distinctUntilChanged()
            .collect {
                val idx = listState.firstVisibleItemIndex
                val off = listState.firstVisibleItemScrollOffset
                val adj = if (off > 0) idx + 1 else idx
                val target = adj.coerceIn(0, options.size - 1)
                if (target != selectedIndex) onSelectedChange(target)
            }
    }

    LaunchedEffect(selectedIndex) {
        val cur = listState.firstVisibleItemIndex
        if (cur != selectedIndex.coerceIn(0, options.size - 1)) {
            scope.launch {
                listState.scrollToItem(selectedIndex.coerceIn(0, options.size - 1))
            }
        }
    }

    Box(modifier.height(itemHeight * visibleItemsCount)) {
        LazyColumn(
            state = listState,
            flingBehavior = flingBehavior,
            modifier = Modifier.fillMaxSize(),
            horizontalAlignment = Alignment.CenterHorizontally,
            contentPadding = PaddingValues(vertical = itemHeight * padCount)
        ) {
            items(count = options.size) { idx ->
                val isSel = idx == listState.firstVisibleItemIndex.let {
                    if (listState.firstVisibleItemScrollOffset > 0) it + 1 else it
                }
                Box(
                    Modifier.fillMaxWidth().height(itemHeight),
                    contentAlignment = Alignment.Center
                ) {
                    Text(
                        options[idx],
                        fontSize = if (isSel) 16.sp else 14.sp,
                        fontWeight = if (isSel) FontWeight.Bold else FontWeight.Normal,
                        color = if (isSel) OnSurface else TextSecondary,
                        textAlign = TextAlign.Center
                    )
                }
            }
        }

        // 中央高亮分隔线
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

        // 上下渐变蒙层占位 (可后续接入 Brush 渐变)
        Box(
            Modifier
                .align(Alignment.TopCenter)
                .fillMaxWidth()
                .height(itemHeight)
                .background(Color.Transparent)
        )
    }
}
