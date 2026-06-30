import SwiftUI
import CoreText

@main
struct SxwnlCalendarApp: App {
    init() {
        // 注册八字命书专用毛笔字体 (与鸿蒙端 BaziPage.ets 同源)
        // 走 CTFontManager 运行时注册, 无需修改 Info.plist 的 UIAppFonts
        if let url = Bundle.main.url(forResource: "WenYue", withExtension: "otf") {
            CTFontManagerRegisterFontsForURL(url as CFURL, .process, nil)
        }
    }

    var body: some Scene {
        WindowGroup {
            TabView {
                CalendarView()
                    .tabItem {
                        Label("月历", systemImage: "calendar")
                    }
                YearCalendarView()
                    .tabItem {
                        Label("年历", systemImage: "calendar.badge.clock")
                    }
                BaziView()
                    .tabItem {
                        Label("八字", systemImage: "person.text.rectangle")
                    }
                EclipseView()
                    .tabItem {
                        Label("日月食", systemImage: "moon.stars")
                    }
            }
        }
    }
}
