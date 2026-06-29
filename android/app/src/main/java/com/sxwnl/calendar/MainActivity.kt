package com.sxwnl.calendar

import android.os.Bundle
import androidx.activity.ComponentActivity
import androidx.activity.compose.setContent
import androidx.activity.enableEdgeToEdge
import androidx.compose.foundation.layout.*
import androidx.compose.material3.*
import androidx.compose.runtime.*
import androidx.compose.ui.Modifier
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
    val tabs = listOf(
        "月历" to "\uD83D\uDCC5",
        "年历" to "\uD83D\uDCC6",
        "八字" to "\u262F",
        "日月食" to "\uD83C\uDF11"
    )

    Scaffold(
        bottomBar = {
            NavigationBar {
                tabs.forEachIndexed { index, (title, icon) ->
                    NavigationBarItem(
                        selected = selectedTab == index,
                        onClick = { selectedTab = index },
                        icon = { Text(icon, style = MaterialTheme.typography.titleLarge) },
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
