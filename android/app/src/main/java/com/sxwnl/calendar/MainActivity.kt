package com.sxwnl.calendar

import android.os.Bundle
import androidx.activity.ComponentActivity
import androidx.activity.compose.setContent
import androidx.activity.enableEdgeToEdge
import androidx.compose.foundation.layout.*
import androidx.compose.material.icons.Icons
import androidx.compose.material.icons.outlined.CalendarMonth
import androidx.compose.material.icons.outlined.CalendarViewMonth
import androidx.compose.material.icons.outlined.DarkMode
import androidx.compose.material.icons.outlined.GridView
import androidx.compose.material3.*
import androidx.compose.runtime.*
import androidx.compose.ui.Modifier
import androidx.compose.ui.graphics.vector.ImageVector
import com.sxwnl.calendar.ui.screens.BaziScreen
import com.sxwnl.calendar.ui.screens.CalendarScreen
import com.sxwnl.calendar.ui.screens.EclipseScreen
import com.sxwnl.calendar.ui.screens.YearCalendarScreen
import com.sxwnl.calendar.ui.theme.SxwnlTheme

class MainActivity : ComponentActivity() {
    override fun onCreate(savedInstanceState: Bundle?) {
        super.onCreate(savedInstanceState)
        enableEdgeToEdge()
        setContent {
            SxwnlTheme {
                MainScreen()
            }
        }
    }
}

@OptIn(ExperimentalMaterial3Api::class)
@Composable
fun MainScreen() {
    var selectedTab by remember { mutableIntStateOf(0) }
    // 统一使用 Material 矢量图标: 单色、同一描边风格, 由 NavigationBar 统一着色,
    // 四个 Tab 表现方式一致 (不再混用方形/圆形 emoji).
    val tabs: List<Pair<String, ImageVector>> = listOf(
        "月历" to Icons.Outlined.CalendarMonth,
        "年历" to Icons.Outlined.CalendarViewMonth,
        "八字" to Icons.Outlined.GridView,
        "日月食" to Icons.Outlined.DarkMode
    )

    Scaffold(
        bottomBar = {
            NavigationBar {
                tabs.forEachIndexed { index, (title, icon) ->
                    NavigationBarItem(
                        selected = selectedTab == index,
                        onClick = { selectedTab = index },
                        icon = { Icon(icon, contentDescription = title) },
                        label = { Text(title) }
                    )
                }
            }
        }
    ) { padding ->
        Box(Modifier.padding(padding).fillMaxSize()) {
            when (selectedTab) {
                0 -> CalendarScreen()
                1 -> YearCalendarScreen()
                2 -> BaziScreen()
                3 -> EclipseScreen()
            }
        }
    }
}
