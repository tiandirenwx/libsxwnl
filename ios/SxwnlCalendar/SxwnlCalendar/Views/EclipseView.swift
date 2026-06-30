import SwiftUI

// ════════════════════════════════════════════════════════════════
//  EclipseView — 日月食工具(SwiftUI 重写版)
//
//  与 Android `EclipseScreen.kt` 端 1:1 对齐, 提供:
//    · 搜索面板 (年/月/数量, 日/月食 segment, Saros 显示)
//    · 卡片列表 (emoji + 渐变徽章 + 食季/Saros 行)
//    · 点击日食 → ModalBottomSheet 3-tab: 食带地图 / 本地观测 / 数据
//    · 点击月食 → ModalBottomSheet 2-tab: 过程动画 / 数据
//    · 截图分享 + ICS 日程导出
// ════════════════════════════════════════════════════════════════

struct EclipseView: View {
    @State private var isSolar: Bool = true
    @State private var searchYear: Int
    @State private var searchMonth: Int = 1
    @State private var searchCount: Int = 8

    @State private var solarItems: [SolarEclipseInfo] = []
    @State private var lunarItems: [LunarEclipseInfo] = []
    @State private var loading: Bool = false
    @State private var searched: Bool = false

    @State private var selectedSolar: SolarEclipseInfo? = nil
    @State private var selectedLunar: LunarEclipseInfo? = nil

    init() {
        let now = Calendar.current.dateComponents([.year], from: Date())
        _searchYear = State(initialValue: now.year ?? 2026)
    }

    var body: some View {
        ScrollView {
            VStack(spacing: 0) {
                headerBanner
                searchPanel
                if searched {
                    if isSolar { solarList } else { lunarList }
                }
            }
        }
        .background(AppColors.background)
        .sheet(item: $selectedSolar) { item in
            SolarDetailSheet(item: item)
                .presentationDetents([.large])
        }
        .sheet(item: $selectedLunar) { item in
            LunarDetailSheet(item: item)
                .presentationDetents([.large])
        }
    }

    // ── 顶部 banner ────────────────────────────────────────
    private var headerBanner: some View {
        VStack(spacing: 4) {
            Text("🌑 日月食工具")
                .font(.system(size: AppDimens.fontTitle, weight: .bold))
                .foregroundColor(AppColors.onPrimary)
            Text("基于天文算法的精确预测 · 动画 · 食带地图 · Saros")
                .font(.system(size: AppDimens.fontCaption))
                .foregroundColor(AppColors.onPrimary.opacity(0.85))
        }
        .frame(maxWidth: .infinity)
        .padding(.top, AppDimens.paddingXl)
        .padding(.bottom, AppDimens.paddingLg)
        .background(
            LinearGradient(colors: [AppColors.gradientStart, AppColors.gradientEnd],
                           startPoint: .topLeading, endPoint: .bottomTrailing)
        )
    }

    // ── 搜索面板 ──────────────────────────────────────────
    private var searchPanel: some View {
        VStack(spacing: 10) {
            HStack(spacing: 0) {
                ForEach([(true, "日食"), (false, "月食")], id: \.0) { v, label in
                    let sel = isSolar == v
                    Button { isSolar = v } label: {
                        Text(label)
                            .font(.system(size: AppDimens.fontBody,
                                          weight: sel ? .bold : .regular))
                            .foregroundColor(sel ? AppColors.onPrimary : AppColors.textSecondary)
                            .frame(maxWidth: .infinity)
                            .frame(height: 36)
                            .background(sel ? AppColors.primary : Color.clear)
                            .clipShape(Capsule())
                    }
                    .buttonStyle(.plain)
                }
            }
            .background(AppColors.background)
            .clipShape(Capsule())

            HStack(spacing: 6) {
                Text("起始")
                    .font(.system(size: AppDimens.fontCaption))
                    .foregroundColor(AppColors.textSecondary)
                numField($searchYear, width: nil).frame(maxWidth: .infinity, minHeight: 38)
                Text("年").font(.system(size: AppDimens.fontBody))
                numField($searchMonth, width: 50).frame(height: 38)
                Text("月起").font(.system(size: AppDimens.fontBody))
                numField($searchCount, width: 50).frame(height: 38)
                Text("条").font(.system(size: AppDimens.fontBody))
            }

            Button(action: doSearch) {
                HStack(spacing: 6) {
                    if loading {
                        ProgressView()
                            .progressViewStyle(.circular)
                            .tint(AppColors.onPrimary)
                            .scaleEffect(0.8)
                    }
                    Text(loading ? "搜索中…" : "开始搜索")
                        .font(.system(size: AppDimens.fontBody, weight: .semibold))
                }
                .foregroundColor(AppColors.onPrimary)
                .frame(maxWidth: .infinity)
                .frame(height: 42)
                .background(AppColors.primary)
                .clipShape(Capsule())
            }
            .buttonStyle(.plain)
            .disabled(loading)
        }
        .padding(AppDimens.paddingMd)
        .background(AppColors.surface)
        .cornerRadius(AppDimens.radiusLg)
        .shadow(color: Color.black.opacity(0.06), radius: 3, y: 2)
        .padding(AppDimens.paddingMd)
    }

    // ── 列表 ────────────────────────────────────────────────
    private var solarList: some View {
        Group {
            if solarItems.isEmpty {
                emptyHint("起始时间范围内没有日食")
            } else {
                VStack(spacing: 8) {
                    ForEach(solarItems.indices, id: \.self) { idx in
                        SolarCard(item: solarItems[idx]) {
                            selectedSolar = solarItems[idx]
                        }
                    }
                }
                .padding(.horizontal, AppDimens.paddingMd)
                .padding(.bottom, 24)
            }
        }
    }

    private var lunarList: some View {
        Group {
            if lunarItems.isEmpty {
                emptyHint("起始时间范围内没有月食")
            } else {
                VStack(spacing: 8) {
                    ForEach(lunarItems.indices, id: \.self) { idx in
                        LunarCard(item: lunarItems[idx]) {
                            selectedLunar = lunarItems[idx]
                        }
                    }
                }
                .padding(.horizontal, AppDimens.paddingMd)
                .padding(.bottom, 24)
            }
        }
    }

    private func emptyHint(_ text: String) -> some View {
        Text(text)
            .font(.system(size: AppDimens.fontBody))
            .foregroundColor(AppColors.textSecondary)
            .frame(maxWidth: .infinity)
            .padding(AppDimens.paddingXl)
    }

    private func numField(_ binding: Binding<Int>, width: CGFloat?) -> some View {
        TextField("", value: binding, format: .number)
            .textFieldStyle(.plain)
            .keyboardType(.numberPad)
            .font(.system(size: AppDimens.fontBody))
            .multilineTextAlignment(.center)
            .padding(.horizontal, 6)
            .frame(width: width)
            .background(AppColors.background)
            .overlay(
                RoundedRectangle(cornerRadius: AppDimens.radiusSm)
                    .stroke(AppColors.divider, lineWidth: 1)
            )
            .clipShape(RoundedRectangle(cornerRadius: AppDimens.radiusSm))
    }

    private func doSearch() {
        loading = true
        searched = true
        let y = searchYear, m = max(1, min(12, searchMonth))
        let c = max(1, min(30, searchCount))
        let solar = isSolar
        DispatchQueue.global(qos: .userInitiated).async {
            if solar {
                let r = SxwnlEclipseMapBridge.searchSolarEclipses(year: y, month: m, count: c)
                DispatchQueue.main.async {
                    solarItems = r
                    lunarItems = []
                    loading = false
                }
            } else {
                let r = SxwnlEclipseMapBridge.searchLunarEclipses(year: y, month: m, count: c)
                DispatchQueue.main.async {
                    lunarItems = r
                    solarItems = []
                    loading = false
                }
            }
        }
    }
}

// ═══════════════════════════════════════════════════════════════
// MARK: - Cards
// ═══════════════════════════════════════════════════════════════

private struct SolarCard: View {
    let item: SolarEclipseInfo
    let onTap: () -> Void

    var body: some View {
        Button(action: onTap) {
            VStack(alignment: .leading, spacing: 8) {
                HStack {
                    Text(EclipseUtil.solarEmoji(item.type)).font(.system(size: 22))
                    VStack(alignment: .leading, spacing: 2) {
                        Text(String(format: "%04d-%02d-%02d %02d:%02d TD",
                                    item.year, item.month, item.day, item.hour, item.minute))
                            .font(.system(size: AppDimens.fontBody, weight: .bold))
                            .foregroundColor(AppColors.onSurface)
                        Text("\(item.season)  ·  Saros \(item.saros)")
                            .font(.system(size: 11))
                            .foregroundColor(AppColors.textSecondary)
                    }
                    Spacer()
                    GradientBadge(text: item.typeName,
                                  start: EclipseUtil.solarBadgeColor(item.type),
                                  end: EclipseUtil.solarBadgeColor(item.type).opacity(0.75),
                                  icon: EclipseUtil.solarEmoji(item.type))
                }

                HStack(spacing: 8) {
                    Text("食分")
                        .font(.system(size: AppDimens.fontCaption))
                        .foregroundColor(AppColors.textSecondary)
                    GeometryReader { geo in
                        ZStack(alignment: .leading) {
                            Capsule().fill(AppColors.divider).frame(height: 6)
                            Capsule().fill(EclipseUtil.solarBadgeColor(item.type))
                                .frame(width: geo.size.width *
                                       CGFloat(EclipseUtil.magnitudeProgress(item.magnitude)),
                                       height: 6)
                        }
                    }
                    .frame(height: 6)
                    Text(String(format: "%.3f", item.magnitude))
                        .font(.system(size: AppDimens.fontCaption, weight: .bold))
                        .foregroundColor(AppColors.onSurface)
                }

                HStack(spacing: 14) {
                    ParamSm(label: "γ", value: String(format: "%+.3f", item.gamma))
                    ParamSm(label: "最长", value: EclipseUtil.formatDuration(item.durationSeconds))
                    ParamSm(label: "食带宽",
                            value: item.pathWidthKm > 0
                                ? String(format: "%.0f km", item.pathWidthKm) : "—")
                    Spacer()
                }
                Text("点击查看 食带地图 / 本地观测 / 数据 →")
                    .font(.system(size: 11, weight: .semibold))
                    .foregroundColor(AppColors.primary)
            }
            .padding(AppDimens.paddingMd)
            .frame(maxWidth: .infinity, alignment: .leading)
            .background(AppColors.surface)
            .cornerRadius(AppDimens.radiusMd)
            .shadow(color: Color.black.opacity(0.05), radius: 2, y: 1)
        }
        .buttonStyle(.plain)
    }
}

private struct LunarCard: View {
    let item: LunarEclipseInfo
    let onTap: () -> Void

    var body: some View {
        Button(action: onTap) {
            VStack(alignment: .leading, spacing: 8) {
                HStack {
                    Text(EclipseUtil.lunarEmoji(item.type)).font(.system(size: 22))
                    VStack(alignment: .leading, spacing: 2) {
                        Text(String(format: "%04d-%02d-%02d %02d:%02d TD",
                                    item.year, item.month, item.day, item.hour, item.minute))
                            .font(.system(size: AppDimens.fontBody, weight: .bold))
                            .foregroundColor(AppColors.onSurface)
                        Text("\(item.season)  ·  Saros \(item.saros)")
                            .font(.system(size: 11))
                            .foregroundColor(AppColors.textSecondary)
                    }
                    Spacer()
                    GradientBadge(text: item.typeName,
                                  start: EclipseUtil.lunarBadgeColor(item.type),
                                  end: EclipseUtil.lunarBadgeColor(item.type).opacity(0.75),
                                  icon: EclipseUtil.lunarEmoji(item.type))
                }
                HStack(spacing: 8) {
                    Text("食分")
                        .font(.system(size: AppDimens.fontCaption))
                        .foregroundColor(AppColors.textSecondary)
                    GeometryReader { geo in
                        ZStack(alignment: .leading) {
                            Capsule().fill(AppColors.divider).frame(height: 6)
                            Capsule().fill(EclipseUtil.lunarBadgeColor(item.type))
                                .frame(width: geo.size.width *
                                       CGFloat(EclipseUtil.magnitudeProgress(item.magnitude)),
                                       height: 6)
                        }
                    }
                    .frame(height: 6)
                    Text(String(format: "%.3f", item.magnitude))
                        .font(.system(size: AppDimens.fontCaption, weight: .bold))
                        .foregroundColor(AppColors.onSurface)
                }
                Text("点击查看 过程动画 / 数据 →")
                    .font(.system(size: 11, weight: .semibold))
                    .foregroundColor(AppColors.primary)
            }
            .padding(AppDimens.paddingMd)
            .frame(maxWidth: .infinity, alignment: .leading)
            .background(AppColors.surface)
            .cornerRadius(AppDimens.radiusMd)
            .shadow(color: Color.black.opacity(0.05), radius: 2, y: 1)
        }
        .buttonStyle(.plain)
    }
}

// ═══════════════════════════════════════════════════════════════
// MARK: - Solar Detail Sheet (3 tabs)
// ═══════════════════════════════════════════════════════════════

private struct SolarDetailSheet: View {
    let item: SolarEclipseInfo

    @State private var tab = 0
    @State private var loading = true
    @State private var path: SolarEclipsePathData? = nil
    @State private var local: LocalSolarEclipseData? = nil
    @State private var city: ObserverCity = Cities.default
    @State private var currentJd: Double = 0
    @State private var isPlaying = false
    @State private var speed: Double = 1.0
    @State private var ditu0: [Double] = []
    @State private var ditu1: [Double] = []
    @State private var ditu2: [Double] = []
    @State private var useBigMap = true
    @State private var shareItems: [Any]? = nil
    @State private var showShare = false

    var body: some View {
        NavigationStack {
            VStack(spacing: 0) {
                tabBar
                ScrollView {
                    VStack(alignment: .leading, spacing: 12) {
                        if loading {
                            ProgressView("加载中…")
                                .frame(maxWidth: .infinity, minHeight: 140)
                        } else {
                            switch tab {
                            case 0: pathMapTab
                            case 1: localTab
                            default: dataTab
                            }
                        }
                    }
                    .padding(AppDimens.paddingMd)
                }
            }
            .navigationTitle("日食详情")
            .navigationBarTitleDisplayMode(.inline)
            .toolbar {
                ToolbarItem(placement: .topBarTrailing) {
                    Menu {
                        Button { shareSnapshot() } label: {
                            Label("分享截图", systemImage: "square.and.arrow.up")
                        }
                        Button { shareICS() } label: {
                            Label("导出日程 (ICS)", systemImage: "calendar.badge.plus")
                        }
                    } label: {
                        Image(systemName: "ellipsis.circle")
                    }
                }
            }
            .sheet(isPresented: $showShare) {
                if let items = shareItems { ShareSheet(items: items) }
            }
        }
        .task { await loadAll() }
    }

    private var tabBar: some View {
        HStack(spacing: 0) {
            ForEach(0..<3) { i in
                let title = ["食带地图", "本地观测", "数据"][i]
                Button { tab = i } label: {
                    VStack(spacing: 4) {
                        Text(title)
                            .font(.system(size: AppDimens.fontBody,
                                          weight: tab == i ? .bold : .regular))
                            .foregroundColor(tab == i ? AppColors.primary : AppColors.textSecondary)
                        Rectangle()
                            .fill(tab == i ? AppColors.primary : Color.clear)
                            .frame(height: 2)
                    }
                    .frame(maxWidth: .infinity)
                    .padding(.vertical, 8)
                }
                .buttonStyle(.plain)
            }
        }
        .background(AppColors.surface)
    }

    // ── Tab 0: 食带地图 ─────────────────────────────────────
    private var pathMapTab: some View {
        VStack(spacing: 10) {
            HStack {
                GradientBadge(text: item.typeName,
                              start: EclipseUtil.solarBadgeColor(item.type),
                              end: EclipseUtil.solarBadgeColor(item.type).opacity(0.75),
                              icon: EclipseUtil.solarEmoji(item.type))
                Spacer()
                Toggle("大图", isOn: $useBigMap)
                    .toggleStyle(.switch)
                    .labelsHidden()
                Text(useBigMap ? "大图" : "小图")
                    .font(.system(size: 11))
                    .foregroundColor(AppColors.textSecondary)
            }
            if let p = path {
                SolarPathMapCanvas(
                    worldMapRad: useBigMap && !ditu1.isEmpty ? ditu1 : ditu0,
                    bordersRad: useBigMap && !ditu2.isEmpty ? ditu2 : nil,
                    path: p,
                    currentJd: currentJd
                )
                EclipseTimelineBar(
                    jdBegin: p.partialStart.julianDay,
                    jdEnd:   p.partialEnd.julianDay,
                    currentJd: $currentJd,
                    marks: [
                        TimelineMark(label: "P₁", jd: p.partialStart.julianDay,
                                     color: AppColors.pathPenumbra),
                        TimelineMark(label: "C₁", jd: p.centralStart.julianDay,
                                     color: AppColors.pathUmbra),
                        TimelineMark(label: "MAX", jd: p.greatestEclipse.julianDay,
                                     color: AppColors.pathCenter),
                        TimelineMark(label: "C₂", jd: p.centralEnd.julianDay,
                                     color: AppColors.pathUmbra),
                        TimelineMark(label: "P₄", jd: p.partialEnd.julianDay,
                                     color: AppColors.pathPenumbra)
                    ]
                )
                PlaybackController(
                    isPlaying: $isPlaying,
                    currentJd: $currentJd,
                    jdBegin: p.partialStart.julianDay,
                    jdEnd:   p.partialEnd.julianDay,
                    speedMultiplier: $speed
                )
            } else {
                Text("此条目无路径数据 (可能是仅可见极区的偏食)")
                    .font(.system(size: AppDimens.fontBody))
                    .foregroundColor(AppColors.textSecondary)
                    .frame(maxWidth: .infinity).padding(20)
            }
        }
    }

    // ── Tab 1: 本地观测 ─────────────────────────────────────
    private var localTab: some View {
        VStack(spacing: 10) {
            LocationPicker(selected: $city)
                .onChange(of: city) { _, _ in reloadLocal() }
            if let l = local, !l.frames.isEmpty {
                SolarLocalDiscCanvas(frames: l.frames,
                                     currentJd: currentJd,
                                     city: city)
                EclipseTimelineBar(
                    jdBegin: l.firstContact > 0 ? l.firstContact : l.frames.first!.julianDay,
                    jdEnd:   l.fourthContact > 0 ? l.fourthContact : l.frames.last!.julianDay,
                    currentJd: $currentJd,
                    marks: [
                        TimelineMark(label: "P₁", jd: l.firstContact,
                                     color: AppColors.pathPenumbra),
                        TimelineMark(label: "C₁", jd: l.secondContact,
                                     color: AppColors.pathUmbra),
                        TimelineMark(label: "MAX", jd: l.maximum,
                                     color: AppColors.pathCenter),
                        TimelineMark(label: "C₂", jd: l.thirdContact,
                                     color: AppColors.pathUmbra),
                        TimelineMark(label: "P₄", jd: l.fourthContact,
                                     color: AppColors.pathPenumbra)
                    ].filter { $0.jd > 0 }
                )
                PlaybackController(
                    isPlaying: $isPlaying,
                    currentJd: $currentJd,
                    jdBegin: l.frames.first!.julianDay,
                    jdEnd:   l.frames.last!.julianDay,
                    speedMultiplier: $speed
                )
                ParamGrid(pairs: [
                    ("初亏 P₁", EclipseUtil.jdTdToLocal(l.firstContact, tzHours: city.tzHours)),
                    ("食既 C₁", EclipseUtil.jdTdToLocal(l.secondContact, tzHours: city.tzHours)),
                    ("食甚",    EclipseUtil.jdTdToLocal(l.maximum, tzHours: city.tzHours)),
                    ("生光 C₂", EclipseUtil.jdTdToLocal(l.thirdContact, tzHours: city.tzHours)),
                    ("复圆 P₄", EclipseUtil.jdTdToLocal(l.fourthContact, tzHours: city.tzHours)),
                    ("最大食分", String(format: "%.3f", l.maxMagnitude))
                ])
            } else {
                Text("该地点此次日食不可见")
                    .font(.system(size: AppDimens.fontBody))
                    .foregroundColor(AppColors.textSecondary)
                    .frame(maxWidth: .infinity).padding(20)
            }
        }
    }

    // ── Tab 2: 数据 ─────────────────────────────────────────
    private var dataTab: some View {
        VStack(alignment: .leading, spacing: 10) {
            HStack {
                Text("基础数据")
                    .font(.system(size: AppDimens.fontSubtitle, weight: .bold))
                    .foregroundColor(AppColors.primary)
                Spacer()
                GradientBadge(text: item.typeName,
                              start: EclipseUtil.solarBadgeColor(item.type),
                              end: EclipseUtil.solarBadgeColor(item.type).opacity(0.75),
                              icon: EclipseUtil.solarEmoji(item.type))
            }
            ParamGrid(pairs: [
                ("日期", String(format: "%04d-%02d-%02d", item.year, item.month, item.day)),
                ("时刻 (TD)", String(format: "%02d:%02d", item.hour, item.minute)),
                ("食分", String(format: "%.4f", item.magnitude)),
                ("γ", String(format: "%+.4f", item.gamma)),
                ("最长时长", EclipseUtil.formatDuration(item.durationSeconds)),
                ("食带宽", item.pathWidthKm > 0
                        ? String(format: "%.1f km", item.pathWidthKm) : "—"),
                ("类型代码", item.type),
                ("中心点", item.hasCenter
                        ? String(format: "%.2f°, %.2f°",
                                 item.centerLongitude, item.centerLatitude)
                        : "无中心"),
                ("Saros", "\(item.saros) (#\(item.sarosMember))"),
                ("食季", item.season)
            ])

            if let p = path {
                Text("全球路径关键点")
                    .font(.system(size: AppDimens.fontSubtitle, weight: .bold))
                    .foregroundColor(AppColors.primary)
                EventPointRow(label: "偏食始 P₁", point: p.partialStart)
                EventPointRow(label: "中心始 C₁", point: p.centralStart)
                EventPointRow(label: "食甚 Greatest", point: p.greatestEclipse)
                EventPointRow(label: "中心终 C₂", point: p.centralEnd)
                EventPointRow(label: "偏食终 P₄", point: p.partialEnd)
            }
        }
    }

    // ── 加载 ────────────────────────────────────────────────
    private func loadAll() async {
        await withCheckedContinuation { (c: CheckedContinuation<Void, Never>) in
            DispatchQueue.global(qos: .userInitiated).async {
                let p = SxwnlEclipseMapBridge.getSolarEclipsePath(
                    year: item.year, month: item.month, day: item.day)
                let d0 = SxwnlEclipseMapBridge.getWorldMapDitu0()
                let d1 = SxwnlEclipseMapBridge.getWorldMapData(idx: 1)
                let d2 = SxwnlEclipseMapBridge.getWorldMapData(idx: 2)
                let l  = SxwnlEclipseMapBridge.getLocalSolarEclipse(
                    year: item.year, month: item.month, day: item.day,
                    longitude: Cities.default.longitude,
                    latitude: Cities.default.latitude,
                    frameInterval: 60)
                DispatchQueue.main.async {
                    self.path = p
                    self.local = l
                    self.ditu0 = d0
                    self.ditu1 = d1
                    self.ditu2 = d2
                    self.currentJd = p?.greatestEclipse.julianDay ?? l?.maximum ?? 0
                    self.loading = false
                }
                c.resume()
            }
        }
    }

    private func reloadLocal() {
        DispatchQueue.global(qos: .userInitiated).async {
            let l = SxwnlEclipseMapBridge.getLocalSolarEclipse(
                year: item.year, month: item.month, day: item.day,
                longitude: city.longitude, latitude: city.latitude,
                frameInterval: 60)
            DispatchQueue.main.async {
                self.local = l
                if let m = l?.maximum, m > 0 { self.currentJd = m }
            }
        }
    }

    // ── 分享 / ICS ──────────────────────────────────────────
    private func shareSnapshot() {
        Task { @MainActor in
            let snap = EclipseShareUtil.snapshot(of:
                SolarShareSnapshot(item: item, path: path, currentJd: currentJd,
                                   ditu0: ditu0, ditu1: ditu1, ditu2: ditu2),
                fileName: "solar_\(item.year)_\(item.month)_\(item.day).png")
            if let url = snap {
                shareItems = [url]
                showShare = true
            }
        }
    }

    private func shareICS() {
        let title = "\(item.typeName) (Saros \(item.saros))"
        let desc = """
        类型: \(item.typeName) (\(item.type))
        食分: \(String(format: "%.3f", item.magnitude))
        食带宽: \(item.pathWidthKm > 0 ? String(format: "%.0f km", item.pathWidthKm) : "—")
        Saros: \(item.saros) (#\(item.sarosMember))
        \(item.season)
        """
        let begin = path?.partialStart.julianDay ?? item.julianDay - 1.0/24
        let end = path?.partialEnd.julianDay ?? item.julianDay + 2.0/24
        let ics = EclipseShareUtil.buildICS(
            uid: "solar-\(item.year)\(item.month)\(item.day)@sxwnl",
            summary: title,
            description: desc,
            startJd: begin,
            endJd: end)
        if let url = EclipseShareUtil.writeICS(ics,
            fileName: "solar_\(item.year)_\(item.month)_\(item.day).ics") {
            shareItems = [url]
            showShare = true
        }
    }
}

// 截图专用快照视图(隐藏交互, 仅展示静态画面)
private struct SolarShareSnapshot: View {
    let item: SolarEclipseInfo
    let path: SolarEclipsePathData?
    let currentJd: Double
    let ditu0: [Double]
    let ditu1: [Double]
    let ditu2: [Double]

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            HStack {
                Text("🌑 \(item.typeName)")
                    .font(.system(size: 18, weight: .bold))
                Spacer()
                Text("Saros \(item.saros)")
                    .font(.system(size: 12))
                    .foregroundColor(.gray)
            }
            Text(String(format: "%04d-%02d-%02d %02d:%02d TD",
                        item.year, item.month, item.day, item.hour, item.minute))
                .font(.system(size: 14))
                .foregroundColor(.gray)
            if let p = path {
                SolarPathMapCanvas(
                    worldMapRad: !ditu1.isEmpty ? ditu1 : ditu0,
                    bordersRad: !ditu2.isEmpty ? ditu2 : nil,
                    path: p,
                    currentJd: currentJd)
                    .frame(height: 240)
            }
            Text("食分 \(String(format: "%.3f", item.magnitude))  ·  食带宽 \(item.pathWidthKm > 0 ? String(format: "%.0f km", item.pathWidthKm) : "—")")
                .font(.system(size: 12))
                .foregroundColor(.gray)
            Text("libsxwnl · 寿星万年历")
                .font(.system(size: 10))
                .foregroundColor(.gray)
        }
        .padding(16)
        .background(Color.white)
    }
}

// ═══════════════════════════════════════════════════════════════
// MARK: - Lunar Detail Sheet (2 tabs)
// ═══════════════════════════════════════════════════════════════

private struct LunarDetailSheet: View {
    let item: LunarEclipseInfo

    @State private var tab = 0
    @State private var loading = true
    @State private var detail: LunarEclipseDetailData? = nil
    @State private var currentJd: Double = 0
    @State private var isPlaying = false
    @State private var speed: Double = 1.0
    @State private var shareItems: [Any]? = nil
    @State private var showShare = false

    var body: some View {
        NavigationStack {
            VStack(spacing: 0) {
                tabBar
                ScrollView {
                    VStack(alignment: .leading, spacing: 12) {
                        if loading {
                            ProgressView("加载中…")
                                .frame(maxWidth: .infinity, minHeight: 140)
                        } else {
                            switch tab {
                            case 0: animTab
                            default: dataTab
                            }
                        }
                    }
                    .padding(AppDimens.paddingMd)
                }
            }
            .navigationTitle("月食详情")
            .navigationBarTitleDisplayMode(.inline)
            .toolbar {
                ToolbarItem(placement: .topBarTrailing) {
                    Menu {
                        Button { shareSnapshot() } label: {
                            Label("分享截图", systemImage: "square.and.arrow.up")
                        }
                        Button { shareICS() } label: {
                            Label("导出日程 (ICS)", systemImage: "calendar.badge.plus")
                        }
                    } label: {
                        Image(systemName: "ellipsis.circle")
                    }
                }
            }
            .sheet(isPresented: $showShare) {
                if let items = shareItems { ShareSheet(items: items) }
            }
        }
        .task { await load() }
    }

    private var tabBar: some View {
        HStack(spacing: 0) {
            ForEach(0..<2) { i in
                let title = ["过程动画", "数据"][i]
                Button { tab = i } label: {
                    VStack(spacing: 4) {
                        Text(title)
                            .font(.system(size: AppDimens.fontBody,
                                          weight: tab == i ? .bold : .regular))
                            .foregroundColor(tab == i ? AppColors.primary : AppColors.textSecondary)
                        Rectangle()
                            .fill(tab == i ? AppColors.primary : Color.clear)
                            .frame(height: 2)
                    }
                    .frame(maxWidth: .infinity)
                    .padding(.vertical, 8)
                }
                .buttonStyle(.plain)
            }
        }
        .background(AppColors.surface)
    }

    private var animTab: some View {
        VStack(spacing: 10) {
            if let d = detail, !d.frames.isEmpty {
                LunarDiscCanvas(frames: d.frames, currentJd: currentJd)
                EclipseTimelineBar(
                    jdBegin: d.penumbraStart > 0 ? d.penumbraStart : d.frames.first!.julianDay,
                    jdEnd:   d.penumbraEnd   > 0 ? d.penumbraEnd   : d.frames.last!.julianDay,
                    currentJd: $currentJd,
                    marks: [
                        TimelineMark(label: "P₁", jd: d.penumbraStart,
                                     color: AppColors.pathPenumbra),
                        TimelineMark(label: "U₁", jd: d.partialStart,
                                     color: AppColors.pathUmbra),
                        TimelineMark(label: "U₂", jd: d.totalStart,
                                     color: AppColors.moonBlood),
                        TimelineMark(label: "MAX", jd: d.maximum,
                                     color: AppColors.pathCenter),
                        TimelineMark(label: "U₃", jd: d.totalEnd,
                                     color: AppColors.moonBlood),
                        TimelineMark(label: "U₄", jd: d.partialEnd,
                                     color: AppColors.pathUmbra),
                        TimelineMark(label: "P₄", jd: d.penumbraEnd,
                                     color: AppColors.pathPenumbra)
                    ].filter { $0.jd > 0 }
                )
                PlaybackController(
                    isPlaying: $isPlaying,
                    currentJd: $currentJd,
                    jdBegin: d.frames.first!.julianDay,
                    jdEnd:   d.frames.last!.julianDay,
                    speedMultiplier: $speed
                )
            }
        }
    }

    private var dataTab: some View {
        VStack(alignment: .leading, spacing: 10) {
            HStack {
                Text("基础数据")
                    .font(.system(size: AppDimens.fontSubtitle, weight: .bold))
                    .foregroundColor(AppColors.primary)
                Spacer()
                GradientBadge(text: item.typeName,
                              start: EclipseUtil.lunarBadgeColor(item.type),
                              end: EclipseUtil.lunarBadgeColor(item.type).opacity(0.75),
                              icon: EclipseUtil.lunarEmoji(item.type))
            }
            ParamGrid(pairs: [
                ("日期", String(format: "%04d-%02d-%02d", item.year, item.month, item.day)),
                ("时刻 (TD)", String(format: "%02d:%02d", item.hour, item.minute)),
                ("食分", String(format: "%.4f", item.magnitude)),
                ("类型代码", item.type),
                ("Saros", "\(item.saros) (#\(item.sarosMember))"),
                ("食季", item.season)
            ])
            if let d = detail {
                Text("关键时刻 (TD)")
                    .font(.system(size: AppDimens.fontSubtitle, weight: .bold))
                    .foregroundColor(AppColors.primary)
                ParamGrid(pairs: [
                    ("半影始 P₁", EclipseUtil.jdToTime(d.penumbraStart)),
                    ("初亏 U₁",   EclipseUtil.jdToTime(d.partialStart)),
                    ("食既 U₂",   EclipseUtil.jdToTime(d.totalStart)),
                    ("食甚",      EclipseUtil.jdToTime(d.maximum)),
                    ("生光 U₃",   EclipseUtil.jdToTime(d.totalEnd)),
                    ("复圆 U₄",   EclipseUtil.jdToTime(d.partialEnd)),
                    ("半影终 P₄", EclipseUtil.jdToTime(d.penumbraEnd))
                ])
            }
        }
    }

    private func load() async {
        await withCheckedContinuation { (c: CheckedContinuation<Void, Never>) in
            DispatchQueue.global(qos: .userInitiated).async {
                let d = SxwnlEclipseMapBridge.getLunarEclipseDetail(
                    year: item.year, month: item.month, day: item.day,
                    frameInterval: 120)
                DispatchQueue.main.async {
                    self.detail = d
                    self.currentJd = d?.maximum ?? item.julianDay
                    self.loading = false
                }
                c.resume()
            }
        }
    }

    private func shareSnapshot() {
        Task { @MainActor in
            let snap = EclipseShareUtil.snapshot(of:
                LunarShareSnapshot(item: item, detail: detail, currentJd: currentJd),
                fileName: "lunar_\(item.year)_\(item.month)_\(item.day).png")
            if let url = snap {
                shareItems = [url]
                showShare = true
            }
        }
    }

    private func shareICS() {
        let title = "\(item.typeName) (Saros \(item.saros))"
        let desc = """
        类型: \(item.typeName) (\(item.type))
        食分: \(String(format: "%.3f", item.magnitude))
        Saros: \(item.saros) (#\(item.sarosMember))
        \(item.season)
        """
        let begin = detail?.penumbraStart ?? item.julianDay - 1.5/24
        let end   = detail?.penumbraEnd   ?? item.julianDay + 1.5/24
        let ics = EclipseShareUtil.buildICS(
            uid: "lunar-\(item.year)\(item.month)\(item.day)@sxwnl",
            summary: title,
            description: desc,
            startJd: begin,
            endJd: end)
        if let url = EclipseShareUtil.writeICS(ics,
            fileName: "lunar_\(item.year)_\(item.month)_\(item.day).ics") {
            shareItems = [url]
            showShare = true
        }
    }
}

private struct LunarShareSnapshot: View {
    let item: LunarEclipseInfo
    let detail: LunarEclipseDetailData?
    let currentJd: Double

    var body: some View {
        VStack(alignment: .leading, spacing: 10) {
            HStack {
                Text("🌕 \(item.typeName)")
                    .font(.system(size: 18, weight: .bold))
                Spacer()
                Text("Saros \(item.saros)")
                    .font(.system(size: 12))
                    .foregroundColor(.gray)
            }
            Text(String(format: "%04d-%02d-%02d %02d:%02d TD",
                        item.year, item.month, item.day, item.hour, item.minute))
                .font(.system(size: 14))
                .foregroundColor(.gray)
            if let d = detail, !d.frames.isEmpty {
                LunarDiscCanvas(frames: d.frames, currentJd: currentJd)
                    .frame(height: 240)
            }
            Text("食分 \(String(format: "%.3f", item.magnitude))  ·  \(item.season)")
                .font(.system(size: 12))
                .foregroundColor(.gray)
            Text("libsxwnl · 寿星万年历")
                .font(.system(size: 10))
                .foregroundColor(.gray)
        }
        .padding(16)
        .background(Color.white)
    }
}

// ═══════════════════════════════════════════════════════════════
// MARK: - Common atoms
// ═══════════════════════════════════════════════════════════════

private struct ParamSm: View {
    let label: String; let value: String
    var body: some View {
        VStack(alignment: .leading, spacing: 1) {
            Text(label).font(.system(size: 10)).foregroundColor(AppColors.textTertiary)
            Text(value).font(.system(size: 12, weight: .medium))
                .foregroundColor(AppColors.onSurface)
        }
    }
}

private struct ParamGrid: View {
    let pairs: [(String, String)]
    var body: some View {
        VStack(spacing: 4) {
            ForEach(0..<((pairs.count + 1) / 2), id: \.self) { row in
                HStack(spacing: 10) {
                    cell(pairs[row * 2])
                    if row * 2 + 1 < pairs.count { cell(pairs[row * 2 + 1]) }
                    else { Spacer().frame(maxWidth: .infinity) }
                }
            }
        }
    }
    private func cell(_ kv: (String, String)) -> some View {
        HStack(spacing: 6) {
            Text(kv.0).font(.system(size: 11))
                .foregroundColor(AppColors.textSecondary)
                .frame(width: 76, alignment: .leading)
            Text(kv.1).font(.system(size: 12, weight: .medium).monospaced())
                .foregroundColor(AppColors.onSurface)
            Spacer()
        }
        .padding(.horizontal, 8).padding(.vertical, 5)
        .background(AppColors.background)
        .clipShape(RoundedRectangle(cornerRadius: 6))
        .frame(maxWidth: .infinity)
    }
}

private struct EventPointRow: View {
    let label: String; let point: GeoPoint
    var body: some View {
        HStack(alignment: .top) {
            Text(label).font(.system(size: AppDimens.fontCaption))
                .foregroundColor(AppColors.textSecondary)
                .frame(width: 110, alignment: .leading)
            VStack(alignment: .leading, spacing: 1) {
                Text(point.longitude < 99
                     ? String(format: "%.2f°,  %.2f°", point.longitude, point.latitude)
                     : "—")
                    .font(.system(size: AppDimens.fontCaption).monospaced())
                    .foregroundColor(AppColors.onSurface)
                Text(EclipseUtil.jdToDateTime(point.julianDay))
                    .font(.system(size: 11).monospaced())
                    .foregroundColor(AppColors.textTertiary)
            }
            Spacer()
        }
        .padding(.horizontal, 10).padding(.vertical, 6)
        .background(AppColors.background)
        .clipShape(RoundedRectangle(cornerRadius: 6))
    }
}

// ─── Identifiable for sheets ────────────────────────────────

extension SolarEclipseInfo: Identifiable {
    var id: Double { julianDay }
}

extension LunarEclipseInfo: Identifiable {
    var id: Double { julianDay }
}
