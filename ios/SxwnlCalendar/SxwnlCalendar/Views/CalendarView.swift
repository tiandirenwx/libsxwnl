import SwiftUI

// ════════════════════════════════════════════════════════════════
//  CalendarView — 与鸿蒙端 CalendarPage.ets 对齐的月历页面
//
//  功能:
//   - 标题栏 (渐变) + 干支年标签
//   - 导航: «/‹/›/», 年/月文本输入(支持 B212), 「今天」
//   - 年信息栏: 干支年 / 生肖 / 黄帝纪年
//   - 网格: 日号 + 朔望弦图标 + 摘要 (节日 > 节气 > 初一月名 > 日名)
//   - 点击 → Sheet 详情 (含: 农历/回历/三柱+纳音/生肖星座建除/节气/月相
//                      /节日/纪年/儒略日)
//   - 底栏: 选中日 RTS + 当月朔望/节气栅格
// ════════════════════════════════════════════════════════════════

struct CalendarView: View {
    @State private var currentYear: Int
    @State private var currentMonth: Int
    @State private var yearInput: String
    @State private var monthInput: String
    @State private var monthDays: [DayInfo] = []
    @State private var selectedDay: DayInfo?

    // 底栏: 用户当前选中日号 (默认今天)
    @State private var rtsDay: Int
    @State private var rts: DayRTS?
    @State private var moonEvents: [MonthEvent] = []
    @State private var jieQiEvents: [MonthEvent] = []

    // ── 地点选择 (对齐鸿蒙 / Android) ─────────────────────────
    //  · 国内选择 (provinces / cities) 改 location → 触发 RTS 重算
    //  · 国际选择 (continents / countries) 只影响顶栏外地实时时钟,
    //    与原 JS sxwnl indexmp.htm 行为一致 (不改观测点)
    @State private var location: GeoCity = GeoCity(
        province: "北京市", district: "天安门",
        longitude: 116.3833, latitude: 39.9, timezone: 8)
    @State private var provinces: [GeoProvince] = []
    @State private var cities: [GeoCity] = []
    @State private var provIdx: Int = 0
    @State private var cityIdx: Int = 0

    @State private var continents: [String] = []
    @State private var countries: [TimezoneGroup] = []
    @State private var contIdx: Int = 0
    @State private var countryIdx: Int = 0
    /// 不参与 SwiftUI 重渲染, 仅缓存全量分组以便切换大洲时过滤
    @State private var allTzGroups: [TimezoneGroup] = []

    // 双时钟: 本地系统时区 + 选中国际时区, 每秒刷新
    @State private var localClock: String = ""
    @State private var intlClock: String = ""

    // 老黄历静态知识 (董公总论/口诀/方位 — 全局只取一次, 懒加载常驻)
    @State private var topics: [AlmanacTopic] = []
    @State private var showTopics: Bool = false
    private let clockTimer = Timer.publish(every: 1, on: .main, in: .common).autoconnect()

    private let todayYear: Int
    private let todayMonth: Int
    private let todayDay: Int

    init() {
        let now = Calendar.current.dateComponents([.year, .month, .day], from: Date())
        let y = now.year ?? 2026
        let m = now.month ?? 1
        let d = now.day ?? 1
        _currentYear = State(initialValue: y)
        _currentMonth = State(initialValue: m)
        _yearInput = State(initialValue: YearUtil.astroYearToStr(y))
        _monthInput = State(initialValue: "\(m)")
        _rtsDay = State(initialValue: d)
        self.todayYear = y
        self.todayMonth = m
        self.todayDay = d
    }

    var body: some View {
        ZStack {
            VStack(spacing: 0) {
                headerSection
                navSection
                topLocationStrip
                yearInfoBar
                weekHeader
                calendarGrid
                bottomInfoBar
            }
            .background(AppColors.background)

            // 点击日期 → 中央深色浮层 (参考 web 版样式)
            if let day = selectedDay {
                Color.black.opacity(0.55)
                    .ignoresSafeArea()
                    .onTapGesture { dismissPopup() }
                    .transition(.opacity)
                DayDetailPopover(
                    day: day,
                    todayYear: todayYear,
                    todayMonth: todayMonth,
                    todayDay: todayDay,
                    onClose: { dismissPopup() },
                    onOpenTopics: {
                        // 切换到"董公择日 · 经典知识": 关闭日详情, 打开静态知识
                        //   切换 (而非叠加) 与鸿蒙 sheetMode='topics' 语义一致.
                        dismissPopup()
                        openTopics()
                    }
                )
                .padding(.horizontal, 32)
                .transition(.scale(scale: 0.92).combined(with: .opacity))
            }
        }
        .animation(.easeOut(duration: 0.18), value: selectedDay?.id)
        .onAppear {
            initTopBars()
            loadMonthData()
            tickClock()
        }
        .onReceive(clockTimer) { _ in tickClock() }
        .sheet(isPresented: $showTopics) {
            AlmanacTopicsSheet(topics: topics) { showTopics = false }
                .presentationDetents([.large])
                .presentationDragIndicator(.visible)
        }
    }

    /// 懒加载静态知识 (一次, 后续命中缓存) 并弹出 sheet.
    private func openTopics() {
        if topics.isEmpty {
            topics = SxwnlBridge.getAlmanacTopics()
        }
        showTopics = true
    }

    // ── 顶部标题栏 ──────────────────────────────────────────
    private var headerSection: some View {
        HStack {
            Text("\(YearUtil.astroYearToStr(currentYear))年\(currentMonth)月")
                .font(.system(size: AppDimens.fontTitle, weight: .bold))
                .foregroundColor(AppColors.onPrimary)
            Spacer()
            if let first = monthDays.first {
                HStack(spacing: 4) {
                    Text("农历")
                        .font(.system(size: AppDimens.fontCaption))
                        .foregroundColor(AppColors.onPrimary.opacity(0.7))
                    Text("\(first.yearGZ)年")
                        .font(.system(size: AppDimens.fontCaption, weight: .medium))
                        .foregroundColor(AppColors.accent)
                }
            }
        }
        .padding(.horizontal, AppDimens.paddingLg)
        .padding(.top, AppDimens.paddingXl)
        .padding(.bottom, AppDimens.paddingSm)
        .background(
            LinearGradient(
                colors: [AppColors.gradientStart, AppColors.gradientEnd],
                startPoint: .leading, endPoint: .trailing)
        )
    }

    // ── 导航条 ──────────────────────────────────────────────
    private var navSection: some View {
        HStack(spacing: 4) {
            navIconButton("«") { navigate(deltaYear: -1, deltaMonth: 0) }
            navIconButton("‹") { navigate(deltaYear: 0, deltaMonth: -1) }

            TextField("YYYY/B212", text: $yearInput)
                .textFieldStyle(.plain)
                .multilineTextAlignment(.center)
                .font(.system(size: AppDimens.fontBody))
                .foregroundColor(AppColors.onPrimary)
                .padding(.horizontal, 6)
                .frame(width: 80, height: 32)
                .background(AppColors.primaryLight.opacity(0.4))
                .cornerRadius(AppDimens.radiusSm)
                .submitLabel(.done)
                .onSubmit { applyYearMonthInput() }

            Text("年")
                .font(.system(size: AppDimens.fontCaption))
                .foregroundColor(AppColors.onPrimary)
                .padding(.horizontal, 2)

            TextField("M", text: $monthInput)
                .textFieldStyle(.plain)
                .multilineTextAlignment(.center)
                .font(.system(size: AppDimens.fontBody))
                .foregroundColor(AppColors.onPrimary)
                .frame(width: 44, height: 32)
                .background(AppColors.primaryLight.opacity(0.4))
                .cornerRadius(AppDimens.radiusSm)
                .keyboardType(.numberPad)
                .submitLabel(.done)
                .onSubmit { applyYearMonthInput() }

            Text("月")
                .font(.system(size: AppDimens.fontCaption))
                .foregroundColor(AppColors.onPrimary)
                .padding(.horizontal, 2)

            navIconButton("›") { navigate(deltaYear: 0, deltaMonth: 1) }
            navIconButton("»") { navigate(deltaYear: 1, deltaMonth: 0) }

            Spacer()

            Button("今天") { goToday() }
                .font(.system(size: AppDimens.fontCaption, weight: .medium))
                .foregroundColor(AppColors.primary)
                .padding(.horizontal, 12)
                .frame(height: 30)
                .background(AppColors.accent)
                .clipShape(Capsule())
        }
        .padding(.horizontal, AppDimens.paddingLg)
        .padding(.top, AppDimens.paddingSm)
        .padding(.bottom, AppDimens.paddingMd)
        .background(
            LinearGradient(
                colors: [AppColors.gradientStart, AppColors.gradientEnd],
                startPoint: .leading, endPoint: .trailing)
        )
    }

    private func navIconButton(_ icon: String, action: @escaping () -> Void) -> some View {
        Button(action: action) {
            Text(icon)
                .font(.system(size: 18, weight: .bold))
                .foregroundColor(AppColors.onPrimary)
                .frame(width: 30, height: 30)
                .background(AppColors.primaryLight.opacity(0.5))
                .clipShape(Circle())
        }
    }

    // ── 年信息栏 ────────────────────────────────────────────
    private var yearInfoBar: some View {
        HStack(spacing: 6) {
            Text("☰")
                .font(.system(size: 14))
                .foregroundColor(AppColors.accent)
            if let first = monthDays.first {
                Text("\(first.yearGZ)年")
                    .font(.system(size: AppDimens.fontBody, weight: .medium))
                    .foregroundColor(AppColors.primary)
                Text("|").foregroundColor(AppColors.divider)
                Text("生肖\(first.shengXiao)")
                    .font(.system(size: AppDimens.fontBody))
                    .foregroundColor(AppColors.primary)
                Spacer()
                if first.huangdiYear > 0 {
                    Text("黄帝\(first.huangdiYear)年")
                        .font(.system(size: AppDimens.fontCaption))
                        .foregroundColor(AppColors.textSecondary)
                }
            } else {
                Spacer()
            }
        }
        .padding(.horizontal, AppDimens.paddingLg)
        .padding(.vertical, AppDimens.paddingSm)
        .background(AppColors.surface)
        .overlay(Rectangle()
            .frame(height: 0.5)
            .foregroundColor(AppColors.divider),
                 alignment: .bottom)
    }

    // ── 星期头 ──────────────────────────────────────────────
    private var weekHeader: some View {
        HStack(spacing: 0) {
            ForEach(Array(AppText.weekNames.enumerated()), id: \.offset) { i, name in
                Text(name)
                    .font(.system(size: AppDimens.fontCaption, weight: .medium))
                    .foregroundColor((i == 0 || i == 6)
                                     ? AppColors.weekend : AppColors.textSecondary)
                    .frame(maxWidth: .infinity)
                    .padding(.vertical, AppDimens.paddingSm)
            }
        }
        .background(AppColors.surface)
        .padding(.horizontal, AppDimens.paddingSm)
    }

    // ── 网格 ────────────────────────────────────────────────
    private var calendarGrid: some View {
        ScrollView {
            LazyVGrid(columns: Array(repeating: GridItem(.flexible(), spacing: 0),
                                     count: 7),
                      spacing: 0) {
                ForEach(0..<getFirstDayOffset(), id: \.self) { _ in
                    Color.clear.frame(height: 64)
                }
                ForEach(monthDays, id: \.solarDay) { day in
                    dayCell(day)
                }
            }
            .padding(.horizontal, AppDimens.paddingSm)
        }
        .frame(maxHeight: .infinity)
        .background(AppColors.surface)
    }

    @ViewBuilder
    private func dayCell(_ day: DayInfo) -> some View {
        let today = isToday(day)
        let isOpen = selectedDay?.solarDay == day.solarDay &&
                     selectedDay?.solarMonth == day.solarMonth &&
                     selectedDay?.solarYear == day.solarYear
        VStack(spacing: 2) {
            HStack(alignment: .center, spacing: 2) {
                Text("\(day.solarDay)")
                    .font(.system(size: AppDimens.fontSubtitle,
                                  weight: today ? .bold : .regular))
                    .foregroundColor(daySolarColor(day))
                moonIcon(for: day)
            }
            Text(cellSubText(day))
                .font(.system(size: AppDimens.fontSmall,
                              weight: (!day.lipuJieQiName.isEmpty || !day.holiday.isEmpty)
                                      ? .medium : .regular))
                .foregroundColor(dayLunarColor(day))
                .lineLimit(1)
        }
        .frame(maxWidth: .infinity, minHeight: 64)
        .background(
            RoundedRectangle(cornerRadius: AppDimens.radiusSm)
                .fill(dayCellBg(day))
        )
        .overlay(
            RoundedRectangle(cornerRadius: AppDimens.radiusSm)
                .strokeBorder(isOpen && !today
                              ? AppColors.accent : Color.clear,
                              lineWidth: 1.5)
        )
        .onTapGesture {
            selectedDay = day
            rtsDay = day.solarDay
            loadRTS()
        }
    }

    @ViewBuilder
    private func moonIcon(for day: DayInfo) -> some View {
        switch day.yueXiangName {
        case "朔":
            Text("●").font(.system(size: 8))
                .foregroundColor(Color(hex: 0x505050))
        case "望":
            Text("●").font(.system(size: 8))
                .foregroundColor(Color(hex: 0xF0B000))
        case "上弦":
            Text("◑").font(.system(size: 8))
                .foregroundColor(Color(hex: 0xF0B000))
        case "下弦":
            Text("◐").font(.system(size: 8))
                .foregroundColor(Color(hex: 0xF0B000))
        default:
            EmptyView()
        }
    }

    private func daySolarColor(_ day: DayInfo) -> Color {
        if isToday(day) { return AppColors.onPrimary }
        if day.isOffDay || !day.holiday.isEmpty { return AppColors.weekend }
        if day.weekDay == 0 || day.weekDay == 6 { return AppColors.weekend }
        return AppColors.onSurface
    }

    private func dayLunarColor(_ day: DayInfo) -> Color {
        if !day.holiday.isEmpty { return AppColors.weekend }
        if !day.lipuJieQiName.isEmpty { return AppColors.jieQi }
        if isToday(day) { return AppColors.onPrimary.opacity(0.8) }
        return AppColors.lunarText
    }

    private func dayCellBg(_ day: DayInfo) -> Color {
        if isToday(day) { return AppColors.primary }
        if selectedDay?.solarDay == day.solarDay &&
            selectedDay?.solarMonth == day.solarMonth &&
            selectedDay?.solarYear == day.solarYear {
            return AppColors.accent.opacity(0.15)
        }
        return Color.clear
    }

    private func cellSubText(_ day: DayInfo) -> String {
        if !day.holiday.isEmpty { return shortName(day.holiday) }
        if !day.major.isEmpty { return shortName(day.major) }
        if !day.lipuJieQiName.isEmpty { return day.lipuJieQiName }
        if day.lunarDayName == "初一" { return day.lunarMonthName }
        return day.lunarDayName
    }

    private func shortName(_ s: String) -> String {
        let first = s.split(whereSeparator: { $0 == " " || $0 == "," || $0 == "，" })
            .first.map(String.init) ?? s
        if first.count > 4 { return String(first.prefix(4)) }
        return first
    }

    private func isToday(_ day: DayInfo) -> Bool {
        day.solarYear == todayYear &&
        day.solarMonth == todayMonth &&
        day.solarDay == todayDay
    }

    private func getFirstDayOffset() -> Int {
        guard let first = monthDays.first else { return 0 }
        return first.weekDay
    }

    // ── 底栏 ────────────────────────────────────────────────
    private var bottomInfoBar: some View {
        VStack(alignment: .leading, spacing: 6) {
            HStack {
                Text("\(YearUtil.astroYearToStr(currentYear))年\(currentMonth)月\(rtsDay)日 · 日月升降")
                    .font(.system(size: AppDimens.fontCaption, weight: .medium))
                    .foregroundColor(AppColors.primary)
                Spacer()
                Text(locationLabel)
                    .font(.system(size: AppDimens.fontSmall))
                    .foregroundColor(AppColors.textSecondary)
            }
            rtsTripleRow(
                l1: "日出", v1: rts?.sunRise ?? "--:--:--",
                l2: "日落", v2: rts?.sunSet ?? "--:--:--",
                l3: "中天", v3: rts?.sunMeridian ?? "--:--:--")
            rtsTripleRow(
                l1: "月出", v1: rts?.moonRise ?? "--:--:--",
                l2: "月落", v2: rts?.moonSet ?? "--:--:--",
                l3: "月中", v3: rts?.moonMeridian ?? "--:--:--")
            rtsDoubleRow(
                l1: "晨起天亮", v1: rts?.civilDawn ?? "--:--:--",
                l2: "晚上天黑", v2: rts?.civilDusk ?? "--:--:--")
            rtsDoubleRow(
                l1: "日照时间", v1: rts?.dayLength ?? "--:--:--",
                l2: "白天时间", v2: rts?.lightLength ?? "--:--:--")

            if !moonEvents.isEmpty || !jieQiEvents.isEmpty {
                Divider().background(AppColors.divider).padding(.vertical, 4)
                Text("\(currentMonth)月月相与节气")
                    .font(.system(size: AppDimens.fontSmall))
                    .foregroundColor(AppColors.textSecondary)
                eventGrid(events: moonEvents, color: AppColors.primaryLight)
                if !jieQiEvents.isEmpty {
                    eventGrid(events: jieQiEvents, color: AppColors.jieQi)
                }
            }
        }
        .padding(AppDimens.paddingMd)
        .background(AppColors.surface)
        .overlay(Rectangle()
            .frame(height: 0.5)
            .foregroundColor(AppColors.divider),
                 alignment: .top)
    }

    private func rtsTripleRow(l1: String, v1: String,
                              l2: String, v2: String,
                              l3: String, v3: String) -> some View {
        HStack(spacing: 0) {
            rtsLabeled(label: l1, value: v1)
            rtsLabeled(label: l2, value: v2)
            rtsLabeled(label: l3, value: v3)
        }
    }

    private func rtsDoubleRow(l1: String, v1: String,
                              l2: String, v2: String) -> some View {
        HStack(spacing: 0) {
            rtsLabeled(label: l1, value: v1)
            rtsLabeled(label: l2, value: v2)
        }
    }

    private func rtsLabeled(label: String, value: String) -> some View {
        HStack(spacing: 4) {
            Text(label)
                .font(.system(size: AppDimens.fontSmall))
                .foregroundColor(AppColors.textSecondary)
            Text(value)
                .font(.system(size: AppDimens.fontCaption, weight: .medium))
                .foregroundColor(value == "--:--:--"
                                 ? AppColors.textSecondary : AppColors.onSurface)
            Spacer(minLength: 0)
        }
        .frame(maxWidth: .infinity, alignment: .leading)
    }

    private func eventGrid(events: [MonthEvent], color: Color) -> some View {
        LazyVGrid(columns: Array(repeating: GridItem(.flexible(), spacing: 0),
                                 count: 2),
                  alignment: .leading, spacing: 2) {
            ForEach(events) { e in
                HStack(spacing: 4) {
                    Text(padNum(e.day, width: 2) + "日")
                        .font(.system(size: AppDimens.fontSmall))
                        .foregroundColor(AppColors.textSecondary)
                    Text(e.time)
                        .font(.system(size: AppDimens.fontSmall))
                        .foregroundColor(AppColors.onSurface)
                    Text(e.name)
                        .font(.system(size: AppDimens.fontSmall, weight: .medium))
                        .foregroundColor(color)
                    Spacer(minLength: 0)
                }
            }
        }
    }

    private func padNum(_ n: Int, width: Int) -> String {
        let s = "\(n)"
        if s.count >= width { return s }
        return String(repeating: " ", count: width - s.count) + s
    }

    // ── 数据加载 ────────────────────────────────────────────

    private func loadMonthData() {
        monthDays = SxwnlBridge.getMonthData(year: currentYear, month: currentMonth)
        yearInput = YearUtil.astroYearToStr(currentYear)
        monthInput = "\(currentMonth)"
        collectMonthEvents()
        if rtsDay > monthDays.count { rtsDay = 1 }
        loadRTS()
    }

    // 浮层关闭: 当前月若为今月 → 日月升降回今天, 否则 → 1 号
    private func dismissPopup() {
        selectedDay = nil
        let resetDay = (currentYear == todayYear && currentMonth == todayMonth) ? todayDay : 1
        if resetDay != rtsDay {
            rtsDay = resetDay
            loadRTS()
        }
    }

    private func loadRTS() {
        rts = SxwnlBridge.calcDayRTS(
            year: currentYear, month: currentMonth, day: rtsDay,
            longitude: location.longitude,
            latitude: location.latitude,
            tzHours: location.timezone
        )
    }

    private func collectMonthEvents() {
        var moon: [MonthEvent] = []
        var jq: [MonthEvent] = []
        for d in monthDays {
            if !d.yueXiangName.isEmpty {
                var label = d.yueXiangName
                if label == "朔" { label = "朔月" }
                else if label == "望" { label = "望月" }
                moon.append(MonthEvent(
                    day: d.solarDay,
                    time: extractTime(d.yueXiangTime),
                    name: label))
            }
            if !d.lipuJieQiName.isEmpty {
                jq.append(MonthEvent(
                    day: d.solarDay,
                    time: extractTime(d.jieQiTime),
                    name: d.lipuJieQiName))
            }
        }
        moonEvents = moon
        jieQiEvents = jq
    }

    private func extractTime(_ s: String) -> String {
        if let r = s.range(of: #"\d{2}:\d{2}:\d{2}"#, options: .regularExpression) {
            return String(s[r])
        }
        return s.trimmingCharacters(in: .whitespaces)
    }

    private func navigate(deltaYear: Int, deltaMonth: Int) {
        var y = currentYear + deltaYear
        var m = currentMonth + deltaMonth
        while m <= 0 { m += 12; y -= 1 }
        while m > 12 { m -= 12; y += 1 }
        currentYear = y
        currentMonth = m
        rtsDay = 1
        loadMonthData()
    }

    private func goToday() {
        currentYear = todayYear
        currentMonth = todayMonth
        rtsDay = todayDay
        loadMonthData()
    }

    private func applyYearMonthInput() {
        let y = YearUtil.yearStrToAstro(yearInput)
        let m = Int(monthInput) ?? 0
        if YearUtil.isAstroYearValid(y) && (1...12).contains(m) {
            currentYear = y
            currentMonth = m
            rtsDay = 1
            loadMonthData()
        } else {
            yearInput = YearUtil.astroYearToStr(currentYear)
            monthInput = "\(currentMonth)"
        }
    }

    // ── 顶部地点切换条 (国际行 + 国内行) ─────────────────────
    //
    //  与 JS 原版 indexmp.htm 顶部一致:
    //    第 1 行: [洲 ▼] [国家 ▼] · UTC±x · 实时时钟  (仅展示, 不改 location)
    //    第 2 行: [省 ▼] [市 ▼]   · 经X 纬Y           (改 location → 重算 RTS)
    private var topLocationStrip: some View {
        VStack(spacing: 0) {
            // 国际行
            HStack(spacing: 6) {
                if !continents.isEmpty {
                    dropdownPicker(
                        title: continents[safe: contIdx] ?? "",
                        options: continents,
                        onSelect: onContinentChange
                    )
                    dropdownPicker(
                        title: countries[safe: countryIdx]?.country ?? "",
                        options: countries.map { $0.country },
                        onSelect: onCountryChange
                    )
                }
                Spacer()
                if let tz = countries[safe: countryIdx]?.timezone {
                    Text(fmtTz(tz))
                        .font(.system(size: AppDimens.fontSmall))
                        .foregroundColor(AppColors.primary)
                }
                Text(intlClock)
                    .font(.system(size: AppDimens.fontSmall, weight: .medium))
                    .foregroundColor(AppColors.accent)
            }
            .padding(.horizontal, 8)
            .padding(.vertical, 4)
            // 国内行
            HStack(spacing: 6) {
                if !provinces.isEmpty {
                    dropdownPicker(
                        title: provinces[safe: provIdx]?.province ?? "",
                        options: provinces.map { $0.province },
                        onSelect: onProvinceChange
                    )
                    dropdownPicker(
                        title: cities[safe: cityIdx]?.district ?? "",
                        options: cities.map { $0.district },
                        onSelect: onCityChange,
                        primary: true
                    )
                }
                Spacer()
                // 国内行右侧: 经纬度 + 本地系统时钟 (与国际行 intlClock 形成对照)
                Text(fmtLonLat(location.longitude, location.latitude))
                    .font(.system(size: AppDimens.fontSmall))
                    .foregroundColor(AppColors.textSecondary)
                Text(localClock)
                    .font(.system(size: AppDimens.fontSmall, weight: .medium))
                    .foregroundColor(AppColors.primary)
            }
            .padding(.horizontal, 8)
            .padding(.vertical, 4)
        }
        .background(AppColors.surface)
        .overlay(Rectangle()
            .frame(height: 0.5)
            .foregroundColor(AppColors.divider),
                 alignment: .bottom)
    }

    /// 紧凑下拉选择: 单击展开 Menu, options 大于 ~200 时性能可接受 (省级仅 32, 单省最多约 300 市).
    @ViewBuilder
    private func dropdownPicker(title: String, options: [String],
                                onSelect: @escaping (Int) -> Void,
                                primary: Bool = false) -> some View {
        Menu {
            ForEach(options.indices, id: \.self) { idx in
                Button(options[idx]) { onSelect(idx) }
            }
        } label: {
            HStack(spacing: 2) {
                Text(title.isEmpty ? "—" : title)
                    .font(.system(size: AppDimens.fontSmall, weight: .medium))
                    .foregroundColor(primary ? AppColors.primary : AppColors.onSurface)
                    .lineLimit(1)
                Text("▼")
                    .font(.system(size: 8))
                    .foregroundColor(AppColors.textSecondary)
            }
            .padding(.horizontal, 8)
            .padding(.vertical, 4)
            .background(AppColors.primaryLight.opacity(0.1))
            .clipShape(RoundedRectangle(cornerRadius: AppDimens.radiusSm))
        }
    }

    // ── 顶栏数据加载 / 切换回调 ────────────────────────────────

    private func initTopBars() {
        // 省/市目录
        let ps = SxwnlBridge.geoListProvinces()
        if !ps.isEmpty {
            provinces = ps
            let def = SxwnlBridge.geoDefault() ?? location
            let pi = ps.firstIndex(where: { $0.province == def.province }) ?? 0
            provIdx = pi
            let cs = SxwnlBridge.geoListCities(province: ps[pi].province)
            cities = cs
            let ci = cs.firstIndex(where: { $0.district == def.district }) ?? 0
            cityIdx = ci
            if !cs.isEmpty { location = cs[ci] } else { location = def }
        }
        // 时区目录
        let tz = SxwnlBridge.geoListTimezones()
        if !tz.isEmpty {
            allTzGroups = tz
            var seen = Set<String>()
            continents = tz.compactMap { seen.insert($0.continent).inserted ? $0.continent : nil }
            let asiaIdx = continents.firstIndex(of: "亚洲") ?? 0
            contIdx = asiaIdx
            let inAsia = tz.filter { $0.continent == continents[asiaIdx] }
            countries = inAsia
            countryIdx = inAsia.firstIndex(where: { $0.country == "中国" }) ?? 0
        }
    }

    private func onContinentChange(_ idx: Int) {
        guard continents.indices.contains(idx) else { return }
        contIdx = idx
        countries = allTzGroups.filter { $0.continent == continents[idx] }
        countryIdx = 0
    }

    private func onCountryChange(_ idx: Int) {
        guard countries.indices.contains(idx) else { return }
        countryIdx = idx
    }

    private func onProvinceChange(_ idx: Int) {
        guard provinces.indices.contains(idx) else { return }
        provIdx = idx
        let cs = SxwnlBridge.geoListCities(province: provinces[idx].province)
        cities = cs
        cityIdx = 0
        if !cs.isEmpty {
            location = cs[0]
            loadRTS()
        }
    }

    private func onCityChange(_ idx: Int) {
        guard cities.indices.contains(idx) else { return }
        cityIdx = idx
        location = cities[idx]
        loadRTS()
    }

    /// 当前位置标签 (顶栏 + 底栏共用)
    private var locationLabel: String {
        let name: String
        if !location.district.isEmpty && location.district != location.province {
            name = "\(location.province)·\(location.district)"
        } else {
            name = location.district
        }
        return "\(name) (\(fmtLonLat(location.longitude, location.latitude)))"
    }

    // ── 时钟 ───────────────────────────────────────────────────

    private func tickClock() {
        let now = Date()
        localClock = formatClock(now, tzSeconds: TimeZone.current.secondsFromGMT(for: now))
        if let tz = countries[safe: countryIdx]?.timezone {
            // 选中国际时区的"当地时间"
            intlClock = formatClock(now, tzSeconds: Int((tz * 3600).rounded()))
        } else {
            intlClock = ""
        }
    }

    private func formatClock(_ d: Date, tzSeconds: Int) -> String {
        var cal = Calendar(identifier: .gregorian)
        cal.timeZone = TimeZone(secondsFromGMT: tzSeconds) ?? .current
        let c = cal.dateComponents([.month, .day, .hour, .minute, .second], from: d)
        return String(format: "%d/%d %02d:%02d:%02d",
                      c.month ?? 0, c.day ?? 0,
                      c.hour ?? 0, c.minute ?? 0, c.second ?? 0)
    }

    private func fmtTz(_ tz: Double) -> String {
        let sign = tz >= 0 ? "+" : "-"
        let absTz = abs(tz)
        let h = Int(absTz)
        let mm = Int((absTz - Double(h)) * 60)
        return mm == 0 ? "UTC\(sign)\(h)" : String(format: "UTC%@%d:%02d", sign, h, mm)
    }

    private func fmtLonLat(_ lon: Double, _ lat: Double) -> String {
        let lonDir = lon >= 0 ? "E" : "W"
        let latDir = lat >= 0 ? "N" : "S"
        return String(format: "%.1f°%@ %.1f°%@", abs(lon), lonDir, abs(lat), latDir)
    }
}

// MARK: - 安全索引
private extension Array {
    subscript(safe i: Int) -> Element? {
        indices.contains(i) ? self[i] : nil
    }
}

// MARK: - DayInfo Identifiable (for selection)

extension DayInfo: Identifiable {
    var id: String { "\(solarYear)-\(solarMonth)-\(solarDay)" }
}

// MARK: - 详情浮层 (深色紧凑卡片, 与 Android/Harmony 一致风格)
//
//  - 纯深色背景 #14163A, 不再用渐变 (避免渲染时被半透明遮罩稀释)
//  - 卡片宽 ≤ 320pt, 内边距 12-14pt, 行间距 1-2pt
//  - 字号: 标题 fontSubtitle, 主要 fontCaption, 次要 fontSmall

private enum PopupPalette {
    static let bg         = Color(red: 0.078, green: 0.086, blue: 0.227)   // #14163A
    static let border     = Color(red: 0.231, green: 0.247, blue: 0.451)   // #3B3F73
    static let text       = Color.white
    static let sub        = Color(red: 0.863, green: 0.878, blue: 0.949)   // #DCE0F2 (更亮、纯不透明)
    static let gold       = Color(red: 1.000, green: 0.835, blue: 0.310)   // #FFD54F
    static let green      = Color(red: 0.502, green: 0.878, blue: 0.655)   // #80E0A7
    static let red        = Color(red: 1.000, green: 0.616, blue: 0.616)   // #FF9D9D
    static let divider    = Color(red: 0.353, green: 0.373, blue: 0.549)   // #5A5F8C (纯色, 不再半透)
}

private struct DayDetailPopover: View {
    let day: DayInfo
    let todayYear: Int
    let todayMonth: Int
    let todayDay: Int
    let onClose: () -> Void
    let onOpenTopics: () -> Void

    // 老黄历按需加载, 切换日期时重新拉取
    @State private var almanac: DayAlmanac?

    var body: some View {
        VStack(spacing: 0) {
            ScrollView {
                VStack(spacing: 2) {
                    titleRow
                        .padding(.bottom, 4)
                    infoLine("黄帝\(day.huangdiYear)年", day.weekName, day.constellationName)
                    infoLine("\(day.yearGZ)年",
                             (day.isLeapMonth ? "闰" : "") + day.lunarMonthName,
                             "\(day.lunarDayName)日")
                    infoLine("\(day.yearGZ)年", "\(day.monthGZ)月", "\(day.dayGZ)日")
                    infoLine(day.yearNaYin, day.monthNaYin, day.dayNaYin,
                             color: PopupPalette.gold, bold: true)

                    dividerLine.padding(.vertical, 4)

                    centerLine([("生肖", day.shengXiao),
                                ("建除", day.jian12Name)])
                    centeredText(
                        "回历[\(day.moslemYear)年\(day.moslemMonth)月\(day.moslemDay)日]",
                        color: PopupPalette.sub)
                    centeredText("JD \(formatJD(day.julianDay))",
                                 color: PopupPalette.sub)

                    if hasEvent {
                        dividerLine.padding(.vertical, 4)
                        if !day.lipuJieQiName.isEmpty {
                            eventLine("🌿", "节气 \(day.lipuJieQiName)",
                                      extractTime(day.jieQiTime), PopupPalette.green)
                        }
                        if !day.yueXiangName.isEmpty {
                            eventLine("🌙", "月相 \(day.yueXiangName)",
                                      extractTime(day.yueXiangTime), PopupPalette.gold)
                        }
                        if !day.holiday.isEmpty {
                            eventLine("🎉", day.holiday, "", PopupPalette.red)
                        }
                        if !day.major.isEmpty {
                            eventLine("🎊", day.major, "", PopupPalette.gold)
                        }
                        if !day.minor.isEmpty {
                            eventLine("·", day.minor, "", PopupPalette.sub)
                        }
                        if !day.misc.isEmpty {
                            eventLine("🍂", day.misc, "", PopupPalette.green)
                        }
                    }

                    if let a = almanac {
                        dividerLine.padding(.vertical, 6)
                        almanacPanel(a)
                    }

                    // 经典知识入口 (与鸿蒙 "📜 老黄历详情 ›" 对齐)
                    Button(action: onOpenTopics) {
                        Text("📜 老黄历经典 · 董公择日要诀 ›")
                            .font(.system(size: AppDimens.fontCaption, weight: .medium))
                            .foregroundColor(PopupPalette.gold)
                            .frame(maxWidth: .infinity)
                            .padding(.vertical, 6)
                    }
                    .buttonStyle(.plain)
                    .padding(.top, 8)
                }
                .padding(.horizontal, 14)
                .padding(.vertical, 12)
            }
        }
        .frame(maxWidth: 320)
        .background(PopupPalette.bg)
        .clipShape(RoundedRectangle(cornerRadius: AppDimens.radiusMd))
        .overlay(
            RoundedRectangle(cornerRadius: AppDimens.radiusMd)
                .strokeBorder(PopupPalette.border, lineWidth: 0.5)
        )
        .shadow(color: .black.opacity(0.45), radius: 16, x: 0, y: 6)
        .onAppear { loadAlmanac() }
        .onChange(of: day.id) { loadAlmanac() }
    }

    private func loadAlmanac() {
        almanac = SxwnlBridge.getAlmanac(
            year: day.solarYear, month: day.solarMonth, day: day.solarDay)
    }

    private var isToday: Bool {
        day.solarYear == todayYear &&
        day.solarMonth == todayMonth &&
        day.solarDay == todayDay
    }

    private var hasEvent: Bool {
        !day.lipuJieQiName.isEmpty || !day.yueXiangName.isEmpty ||
        !day.holiday.isEmpty || !day.major.isEmpty ||
        !day.minor.isEmpty  || !day.misc.isEmpty
    }

    private var titleRow: some View {
        HStack(spacing: 6) {
            Text("\(YearUtil.astroYearToStr(day.solarYear))年\(day.solarMonth)月\(day.solarDay)日")
                .font(.system(size: AppDimens.fontSubtitle, weight: .bold))
                .foregroundColor(PopupPalette.text)
                .frame(maxWidth: .infinity, alignment: .center)
            if isToday {
                Text("今天")
                    .font(.system(size: AppDimens.fontSmall, weight: .medium))
                    .foregroundColor(AppColors.primary)
                    .padding(.horizontal, 5)
                    .padding(.vertical, 1)
                    .background(AppColors.accent)
                    .cornerRadius(4)
            }
            Button(action: onClose) {
                Image(systemName: "xmark")
                    .font(.system(size: 10, weight: .semibold))
                    .foregroundColor(PopupPalette.sub)
                    .frame(width: 18, height: 18)
                    .background(Circle().fill(Color.white.opacity(0.08)))
            }
        }
    }

    private func infoLine(_ a: String, _ b: String, _ c: String,
                          color: Color = PopupPalette.text,
                          bold: Bool = false) -> some View {
        // 居中紧凑三列: 不再各占 1/3 整宽, 改为水平居中 + 12pt 间隔
        HStack(spacing: 12) {
            popupCell(a, color: color, bold: bold)
            popupCell(b, color: color, bold: bold)
            popupCell(c, color: color, bold: bold)
        }
        .frame(maxWidth: .infinity)
        .padding(.vertical, 1)
    }

    private func popupCell(_ text: String, color: Color, bold: Bool = false) -> some View {
        Text(text.isEmpty ? "—" : text)
            .font(.system(size: AppDimens.fontCaption,
                          weight: bold ? .medium : .regular))
            .foregroundColor(text.isEmpty ? PopupPalette.sub.opacity(0.4) : color)
            .lineLimit(1)
            .multilineTextAlignment(.center)
    }

    private func centerLine(_ items: [(String, String)]) -> some View {
        HStack(spacing: 0) {
            ForEach(items.indices, id: \.self) { i in
                let item = items[i]
                if i > 0 {
                    Text("·").foregroundColor(PopupPalette.sub)
                        .padding(.horizontal, 6)
                        .font(.system(size: AppDimens.fontCaption))
                }
                Text(item.0).foregroundColor(PopupPalette.sub)
                    .font(.system(size: AppDimens.fontSmall))
                    .padding(.trailing, 4)
                Text(item.1).foregroundColor(PopupPalette.text)
                    .font(.system(size: AppDimens.fontCaption, weight: .medium))
            }
        }
        .frame(maxWidth: .infinity)
        .padding(.vertical, 1)
    }

    private func centeredText(_ text: String, color: Color,
                              size: CGFloat = AppDimens.fontSmall) -> some View {
        Text(text)
            .font(.system(size: size))
            .foregroundColor(color)
            .lineLimit(1)
            .frame(maxWidth: .infinity)
            .padding(.vertical, 1)
    }

    private func eventLine(_ icon: String, _ text: String,
                           _ time: String, _ color: Color) -> some View {
        HStack(spacing: 4) {
            Text(icon).font(.system(size: 11)).foregroundColor(color)
            Text(text)
                .font(.system(size: AppDimens.fontCaption, weight: .medium))
                .foregroundColor(color)
            if !time.isEmpty {
                Text(time)
                    .font(.system(size: AppDimens.fontSmall))
                    .foregroundColor(PopupPalette.sub)
            }
        }
        .frame(maxWidth: .infinity)
        .padding(.vertical, 1)
    }

    private var dividerLine: some View {
        Rectangle()
            .fill(PopupPalette.divider)
            .frame(height: 0.5)
    }

    private func extractTime(_ s: String) -> String {
        if let r = s.range(of: #"\d{2}:\d{2}:\d{2}"#, options: .regularExpression) {
            return String(s[r])
        }
        return s.trimmingCharacters(in: .whitespaces)
    }

    /// 与鸿蒙 formatJD 一致: "完整JD(J2000起算)"
    private func formatJD(_ jd: Double) -> String {
        if jd.isNaN { return "-" }
        let fullJD = Int(jd.rounded())
        let d0 = fullJD - 2451545
        return "\(fullJD)(\(d0))"
    }

    // ── 老黄历面板 (内联在 popover 滚动容器中) ─────────────────
    //
    //  紧凑展示: 二十八宿 / 黄道黑道 / 冲煞 / 五吉神方位 / 彭祖 / 神煞 /
    //  宜忌 / 吉时 / 用事择吉 / 节气提示 / 董公择日要诀语录.

    private static let zhiNames = ["子","丑","寅","卯","辰","巳",
                                    "午","未","申","酉","戌","亥"]

    @ViewBuilder
    private func almanacPanel(_ a: DayAlmanac) -> some View {
        Text("老黄历")
            .font(.system(size: AppDimens.fontCaption, weight: .medium))
            .foregroundColor(PopupPalette.gold)
            .frame(maxWidth: .infinity)
            .padding(.bottom, 2)
        infoLine("\(a.xiu)宿(\(a.xiuZheng)\(a.xiuAnimal))",
                 a.twelveGod, a.huangHei,
                 color: a.isHuangDao ? PopupPalette.green : PopupPalette.red,
                 bold: true)
        infoLine("冲\(a.chongShengXiao)", a.chongGanZhi, "煞\(a.sha)")
        // 五方位 (喜/财/福 + 阳贵/阴贵), 与鸿蒙 AlmanacComponents 对齐
        dirRow5([
            ("喜",   a.xiShenFang),
            ("财",   a.caiShenFang),
            ("福",   a.fuShenFang),
            ("阳贵", a.yangGuiFang),
            ("阴贵", a.yinGuiFang)
        ])
        if !a.pengZuGan.isEmpty {
            centeredText(a.pengZuGan, color: PopupPalette.sub)
        }
        if !a.pengZuZhi.isEmpty {
            centeredText(a.pengZuZhi, color: PopupPalette.sub)
        }

        // 神煞: 完整展示, 按权重降序 (与鸿蒙一致). 弹窗 ScrollView 承接溢出.
        let lucky = a.shenSha.filter { $0.isLucky }.sorted { $0.weight > $1.weight }
        let bad   = a.shenSha.filter { !$0.isLucky }.sorted { $0.weight > $1.weight }
        if !lucky.isEmpty { shenShaRow("吉神", lucky, PopupPalette.green) }
        if !bad.isEmpty   { shenShaRow("凶煞", bad,   PopupPalette.red) }
        if !a.yi.isEmpty { textLine("宜", a.yi, PopupPalette.green) }
        if !a.ji.isEmpty { textLine("忌", a.ji, PopupPalette.red) }

        if !a.luckyHours.isEmpty {
            HStack(alignment: .top, spacing: 6) {
                Text("吉时")
                    .font(.system(size: AppDimens.fontSmall))
                    .foregroundColor(PopupPalette.sub)
                Text(a.luckyHours.map { lh -> String in
                    let z = (0..<12).contains(lh.zhi) ? Self.zhiNames[lh.zhi] : ""
                    return "\(lh.name)(\(z))"
                }.joined(separator: "  "))
                    .font(.system(size: AppDimens.fontSmall))
                    .foregroundColor(PopupPalette.gold)
                Spacer(minLength: 0)
            }
            .frame(maxWidth: .infinity, alignment: .leading)
            .padding(.vertical, 1)
        }

        // 用事择吉 — 完整展示 (典型 4-8 条, 不再截断)
        if !a.events.isEmpty {
            centeredText("用事择吉", color: PopupPalette.sub)
            ForEach(a.events.indices, id: \.self) { i in
                eventAdviceRow(a.events[i])
            }
        }
        if !a.notes.isEmpty {
            ForEach(a.notes.indices, id: \.self) { i in
                centeredText("· \(a.notes[i])", color: PopupPalette.sub)
            }
        }
        // 董公择日要诀语录 — 完整展示 (与鸿蒙一致)
        if !a.quotes.isEmpty {
            ForEach(a.quotes.indices, id: \.self) { i in
                quoteCard(a.quotes[i])
            }
        }
    }

    /// 紧凑 5 列方位 (喜/财/福/阳贵/阴贵), 等间距均分. dir 为空时显示 "—".
    @ViewBuilder
    private func dirRow5(_ items: [(String, String)]) -> some View {
        HStack(spacing: 0) {
            ForEach(items.indices, id: \.self) { i in
                let (label, dir) = items[i]
                Text("\(label)\(dir.isEmpty ? "—" : dir)")
                    .font(.system(size: AppDimens.fontSmall, weight: .medium))
                    .foregroundColor(dir.isEmpty
                                     ? PopupPalette.sub.opacity(0.4)
                                     : PopupPalette.gold)
                    .lineLimit(1)
                    .frame(maxWidth: .infinity)
            }
        }
        .padding(.vertical, 1)
    }

    private func shenShaRow(_ label: String, _ items: [ShenSha], _ color: Color) -> some View {
        HStack(alignment: .top, spacing: 6) {
            Text(label)
                .font(.system(size: AppDimens.fontSmall))
                .foregroundColor(PopupPalette.sub)
            Text(items.map { $0.name }.joined(separator: "  "))
                .font(.system(size: AppDimens.fontSmall))
                .foregroundColor(color)
            Spacer(minLength: 0)
        }
        .frame(maxWidth: .infinity, alignment: .leading)
        .padding(.vertical, 1)
    }

    private func textLine(_ label: String, _ items: [String], _ color: Color) -> some View {
        HStack(alignment: .top, spacing: 6) {
            Text(label)
                .font(.system(size: AppDimens.fontCaption, weight: .medium))
                .foregroundColor(color)
            Text(items.joined(separator: "、"))
                .font(.system(size: AppDimens.fontCaption))
                .foregroundColor(PopupPalette.text)
            Spacer(minLength: 0)
        }
        .frame(maxWidth: .infinity, alignment: .leading)
        .padding(.vertical, 1)
    }

    private func eventAdviceRow(_ e: EventAdvice) -> some View {
        HStack(spacing: 4) {
            Text(e.suitable ? "✓" : "✗")
                .font(.system(size: AppDimens.fontCaption, weight: .bold))
                .foregroundColor(e.suitable ? PopupPalette.green : PopupPalette.red)
            Text(e.event)
                .font(.system(size: AppDimens.fontCaption))
                .foregroundColor(PopupPalette.text)
            if !e.reason.isEmpty {
                Text("(\(e.reason))")
                    .font(.system(size: AppDimens.fontSmall))
                    .foregroundColor(PopupPalette.sub)
            }
            Spacer(minLength: 0)
        }
        .frame(maxWidth: .infinity, alignment: .leading)
        .padding(.vertical, 1)
    }

    private func quoteCard(_ q: AlmanacQuote) -> some View {
        let accent: Color = {
            switch q.luck {
            case "吉": return PopupPalette.green
            case "凶": return PopupPalette.red
            case "混": return PopupPalette.gold
            default:   return PopupPalette.sub
            }
        }()
        return VStack(alignment: .leading, spacing: 2) {
            HStack {
                Text(q.title)
                    .font(.system(size: AppDimens.fontSmall, weight: .medium))
                    .foregroundColor(accent)
                Spacer()
                if !q.luck.isEmpty {
                    Text(q.luck)
                        .font(.system(size: AppDimens.fontSmall, weight: .bold))
                        .foregroundColor(accent)
                }
            }
            if !q.source.isEmpty {
                Text("— \(q.source)")
                    .font(.system(size: 9))
                    .foregroundColor(PopupPalette.sub)
            }
            if !q.body.isEmpty {
                Text(q.body)
                    .font(.system(size: AppDimens.fontSmall))
                    .foregroundColor(PopupPalette.text)
                    .lineLimit(6)
            }
        }
        .padding(8)
        .background(PopupPalette.border.opacity(0.3))
        .clipShape(RoundedRectangle(cornerRadius: AppDimens.radiusSm))
        .padding(.top, 4)
    }
}

// MARK: - 老黄历静态知识 Sheet (董公总论/择日基础理论/口诀/方位 等)
//
//  与鸿蒙 AlmanacTopicsSheet / Android AlmanacTopicsDialog 对齐:
//   · 顶部水平 Tab 选择 category, 内容区按选中分类滚动展示
//   · topics 通常 ≤ 64 条, 全量渲染无性能问题
//   · category 保序去重 — 显示顺序与 C++ 注册顺序一致

private struct AlmanacTopicsSheet: View {
    let topics: [AlmanacTopic]
    let onClose: () -> Void

    // 按 category 保序去重 (LinkedHashSet 等价语义)
    private var categories: [String] {
        var seen = Set<String>()
        var out: [String] = []
        for t in topics {
            let c = t.category.isEmpty ? "未分类" : t.category
            if seen.insert(c).inserted { out.append(c) }
        }
        return out
    }

    @State private var selectedCat: String = ""

    var body: some View {
        VStack(spacing: 0) {
            HStack {
                Text("📜 老黄历经典 · 董公择日")
                    .font(.system(size: AppDimens.fontTitle, weight: .bold))
                    .foregroundColor(AppColors.onPrimary)
                Spacer()
                Button(action: onClose) {
                    Text("✕")
                        .font(.system(size: 18, weight: .medium))
                        .foregroundColor(AppColors.onPrimary.opacity(0.8))
                }
                .buttonStyle(.plain)
            }
            .padding(.horizontal, AppDimens.paddingLg)
            .padding(.top, AppDimens.paddingLg)
            .padding(.bottom, AppDimens.paddingMd)
            .background(
                LinearGradient(colors: [AppColors.gradientStart, AppColors.gradientEnd],
                               startPoint: .topLeading, endPoint: .bottomTrailing)
            )

            if topics.isEmpty {
                Spacer()
                Text("暂无数据")
                    .font(.system(size: AppDimens.fontCaption))
                    .foregroundColor(AppColors.textSecondary)
                Spacer()
            } else {
                // 横向分类导航
                ScrollView(.horizontal, showsIndicators: false) {
                    HStack(spacing: 6) {
                        ForEach(categories, id: \.self) { cat in
                            let sel = cat == selectedCat
                            Button {
                                selectedCat = cat
                            } label: {
                                Text(cat)
                                    .font(.system(size: AppDimens.fontSmall,
                                                  weight: sel ? .bold : .regular))
                                    .foregroundColor(sel ? AppColors.primary
                                                          : AppColors.textSecondary)
                                    .padding(.horizontal, 10)
                                    .padding(.vertical, 4)
                                    .background(sel
                                        ? AppColors.primary.opacity(0.12)
                                        : Color.clear)
                                    .clipShape(RoundedRectangle(cornerRadius: AppDimens.radiusSm))
                            }
                            .buttonStyle(.plain)
                        }
                    }
                    .padding(.horizontal, AppDimens.paddingLg)
                }
                .padding(.vertical, 6)
                Divider()

                ScrollView {
                    LazyVStack(alignment: .leading, spacing: 12) {
                        let filtered = topics.filter {
                            (($0.category.isEmpty ? "未分类" : $0.category)) == selectedCat
                        }
                        ForEach(filtered.indices, id: \.self) { i in
                            let t = filtered[i]
                            VStack(alignment: .leading, spacing: 4) {
                                Text(t.title)
                                    .font(.system(size: AppDimens.fontCaption, weight: .bold))
                                    .foregroundColor(AppColors.primary)
                                Text(t.body)
                                    .font(.system(size: AppDimens.fontSmall))
                                    .foregroundColor(AppColors.onSurface)
                                    .lineSpacing(3)
                                    .fixedSize(horizontal: false, vertical: true)
                            }
                        }
                    }
                    .padding(AppDimens.paddingLg)
                }
            }
        }
        .background(AppColors.background)
        .onAppear {
            if selectedCat.isEmpty { selectedCat = categories.first ?? "" }
        }
    }
}
