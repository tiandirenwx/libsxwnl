import SwiftUI

// ════════════════════════════════════════════════════════════════
//  BaziDateTimePicker — 与鸿蒙 components/DateTimePicker.ets 对齐
//  三 tab: 公历 / 农历 / 四柱反推
//  选择完成后通过 onConfirm 回传一个 DateSelection (反推统一为公历)
// ════════════════════════════════════════════════════════════════

struct BaziDateTimePicker: View {
    @Binding var isPresented: Bool
    let initialSelection: DateSelection
    let onConfirm: (DateSelection) -> Void

    @State private var inputMode: BirthInputMode
    @State private var year: Int
    @State private var yearText: String
    @State private var month: Int
    @State private var day: Int
    @State private var hour: Int
    @State private var minute: Int
    @State private var isLeap: Bool
    @State private var isSpec: Bool
    @State private var lunarMonths: [LunarMonth] = []
    @State private var validDays: [Int] = []

    // 四柱反推
    @State private var revIdx: [Int] = [0, 2, 0, 0]
    @State private var revHourUnknown: Bool = false
    @State private var revStartText: String = "1900"
    @State private var revEndText: String = "2100"
    @State private var revResults: [ReverseItem] = []
    @State private var revSearched: Bool = false
    @State private var revPickerIdx: Int? = nil

    private let tabs = ["公历", "农历", "四柱反推"]
    private let pillarNames = ["年柱", "月柱", "日柱", "时柱"]

    init(isPresented: Binding<Bool>,
         initialSelection: DateSelection,
         onConfirm: @escaping (DateSelection) -> Void) {
        _isPresented = isPresented
        self.initialSelection = initialSelection
        self.onConfirm = onConfirm
        let init0 = initialSelection
        _inputMode = State(initialValue:
            init0.inputMode == .lunar ? .lunar : .solar)
        _year = State(initialValue: init0.year)
        _yearText = State(initialValue: "\(init0.year)")
        _month = State(initialValue: init0.month)
        _day = State(initialValue: init0.day)
        _hour = State(initialValue: init0.hour)
        _minute = State(initialValue: init0.minute)
        _isLeap = State(initialValue: init0.isLeap)
        _isSpec = State(initialValue: init0.isSpec)
    }

    var body: some View {
        VStack(spacing: 0) {
            topBar
            Divider().background(AppColors.divider)
            if inputMode == .reverse {
                reverseBody
            } else {
                dateWheels
            }
        }
        .background(AppColors.surface)
        .clipShape(RoundedCorner(radius: AppDimens.radiusLg, corners: [.topLeft, .topRight]))
        .onAppear { refreshAll() }
    }

    // ── 顶部条 ────────────────────────────────────────────
    private var topBar: some View {
        HStack {
            Button("取消") { isPresented = false }
                .font(.system(size: AppDimens.fontBody))
                .foregroundColor(AppColors.textSecondary)

            Spacer()

            HStack(spacing: 4) {
                ForEach(Array(tabs.enumerated()), id: \.offset) { idx, t in
                    Button {
                        switchTab(idx)
                    } label: {
                        Text(t)
                            .font(.system(size: AppDimens.fontCaption,
                                          weight: inputMode.rawValue == idx
                                          ? .bold : .regular))
                            .foregroundColor(inputMode.rawValue == idx
                                             ? AppColors.onPrimary
                                             : AppColors.textSecondary)
                            .padding(.horizontal, 12)
                            .padding(.vertical, 6)
                            .background(inputMode.rawValue == idx
                                        ? AppColors.accent : .clear)
                            .cornerRadius(AppDimens.radiusSm)
                    }
                }
            }

            Spacer()

            if inputMode != .reverse {
                Button("确认") { confirm() }
                    .font(.system(size: AppDimens.fontBody, weight: .bold))
                    .foregroundColor(AppColors.accent)
            } else {
                Text("确认").hidden()
            }
        }
        .padding(.horizontal, AppDimens.paddingLg)
        .padding(.vertical, AppDimens.paddingMd)
    }

    // ── 公历/农历: 年输入 + 月/日/时/分滚轮 ─────────────
    private var dateWheels: some View {
        VStack(spacing: 0) {
            HStack(spacing: 8) {
                Text("年份")
                    .font(.system(size: AppDimens.fontCaption))
                    .foregroundColor(AppColors.textSecondary)
                TextField("如 1990 / B212 / -211", text: $yearText)
                    .textFieldStyle(.plain)
                    .font(.system(size: AppDimens.fontBody))
                    .padding(.horizontal, 10)
                    .frame(maxWidth: .infinity, maxHeight: 38)
                    .background(AppColors.background)
                    .cornerRadius(AppDimens.radiusSm)
                    .onChange(of: yearText) { _ in onYearInput() }
                    .submitLabel(.done)
                Text(BaziCalc.formatYear(year))
                    .font(.system(size: AppDimens.fontCaption))
                    .foregroundColor(AppColors.primary)
            }
            .padding(.horizontal, AppDimens.paddingLg)
            .padding(.top, AppDimens.paddingMd)

            HStack(spacing: 0) {
                if inputMode == .lunar {
                    Picker("", selection: $month) {
                        ForEach(lunarMonths.indices, id: \.self) { i in
                            Text(lunarMonths[i].name).tag(lunarMonths[i].month)
                        }
                    }
                    .pickerStyle(.wheel)
                    .frame(maxWidth: .infinity)
                    .onChange(of: month) { _ in
                        if let m = lunarMonths.first(where: { $0.month == month }) {
                            isLeap = m.isLeap
                            isSpec = m.isSpec
                        }
                        refreshValidDays()
                    }
                } else {
                    Picker("", selection: $month) {
                        ForEach(1...12, id: \.self) { Text("\($0)月").tag($0) }
                    }
                    .pickerStyle(.wheel)
                    .frame(maxWidth: .infinity)
                    .onChange(of: month) { _ in
                        isLeap = false; isSpec = false
                        refreshValidDays()
                    }
                }
                Picker("", selection: $day) {
                    ForEach(validDays, id: \.self) { Text("\($0)日").tag($0) }
                }
                .pickerStyle(.wheel)
                .frame(maxWidth: .infinity)

                Picker("", selection: $hour) {
                    ForEach(0..<24, id: \.self) { Text("\($0)时").tag($0) }
                }
                .pickerStyle(.wheel)
                .frame(maxWidth: .infinity)

                Picker("", selection: $minute) {
                    ForEach(0..<60, id: \.self) { Text("\($0)分").tag($0) }
                }
                .pickerStyle(.wheel)
                .frame(maxWidth: .infinity)
            }
            .frame(height: 180)
            .padding(.horizontal, AppDimens.paddingSm)
        }
        .padding(.bottom, AppDimens.paddingLg)
    }

    // ── 四柱反推 ──────────────────────────────────────────
    private var reverseBody: some View {
        VStack(spacing: 0) {
            HStack(spacing: 6) {
                ForEach(Array(pillarNames.enumerated()), id: \.offset) { i, name in
                    pillarChip(i, name)
                }
            }
            .padding(.top, AppDimens.paddingMd)

            HStack {
                Toggle("", isOn: $revHourUnknown)
                    .labelsHidden()
                    .toggleStyle(SwitchToggleStyle(tint: AppColors.accent))
                Text("时辰未知")
                    .font(.system(size: AppDimens.fontCaption))
                    .foregroundColor(AppColors.textSecondary)
                Spacer()
            }
            .padding(.top, AppDimens.paddingMd)

            HStack(spacing: 6) {
                Text("搜索范围")
                    .font(.system(size: AppDimens.fontCaption))
                    .foregroundColor(AppColors.textSecondary)
                TextField("", text: $revStartText)
                    .textFieldStyle(.plain)
                    .font(.system(size: AppDimens.fontCaption))
                    .frame(width: 80, height: 36)
                    .padding(.horizontal, 6)
                    .background(AppColors.background)
                    .cornerRadius(AppDimens.radiusSm)
                Text("—")
                    .font(.system(size: AppDimens.fontCaption))
                    .foregroundColor(AppColors.textSecondary)
                TextField("", text: $revEndText)
                    .textFieldStyle(.plain)
                    .font(.system(size: AppDimens.fontCaption))
                    .frame(width: 80, height: 36)
                    .padding(.horizontal, 6)
                    .background(AppColors.background)
                    .cornerRadius(AppDimens.radiusSm)
                Spacer()
                Button("反推") { doReverse() }
                    .font(.system(size: AppDimens.fontCaption))
                    .foregroundColor(AppColors.primaryDark)
                    .padding(.horizontal, 12)
                    .frame(height: 36)
                    .background(AppColors.accent)
                    .cornerRadius(AppDimens.radiusSm)
            }
            .padding(.top, AppDimens.paddingMd)

            if revSearched {
                Text("匹配结果（\(revResults.count) 条\(revResults.count >= 60 ? "，已截断" : "")），点击选用")
                    .font(.system(size: AppDimens.fontSmall))
                    .foregroundColor(AppColors.textSecondary)
                    .frame(maxWidth: .infinity, alignment: .leading)
                    .padding(.top, AppDimens.paddingSm)

                ScrollView {
                    VStack(spacing: 0) {
                        ForEach(Array(revResults.enumerated()), id: \.offset) { _, item in
                            Button { applyReverse(item) } label: {
                                HStack {
                                    Text(formatReverseItem(item))
                                        .font(.system(size: AppDimens.fontCaption))
                                        .foregroundColor(AppColors.onSurface)
                                    Spacer()
                                    Text("选用 ›")
                                        .font(.system(size: AppDimens.fontCaption))
                                        .foregroundColor(AppColors.accent)
                                }
                                .padding(.vertical, 8)
                                .overlay(Rectangle()
                                    .frame(height: 0.5)
                                    .foregroundColor(AppColors.divider),
                                         alignment: .bottom)
                            }
                            .buttonStyle(.plain)
                        }
                    }
                }
                .frame(height: 160)
            }
        }
        .padding(.horizontal, AppDimens.paddingLg)
        .padding(.bottom, AppDimens.paddingLg)
        .sheet(item: Binding(
            get: { revPickerIdx.map { PillarIndex(idx: $0) } },
            set: { revPickerIdx = $0?.idx }
        )) { wrap in
            pillarChooser(wrap.idx)
                .presentationDetents([.medium])
        }
    }

    private struct PillarIndex: Identifiable {
        let idx: Int
        var id: Int { idx }
    }

    private func pillarChip(_ i: Int, _ name: String) -> some View {
        Button {
            if !(i == 3 && revHourUnknown) {
                revPickerIdx = i
            }
        } label: {
            VStack(spacing: 4) {
                Text(name)
                    .font(.system(size: AppDimens.fontSmall))
                    .foregroundColor(AppColors.textSecondary)
                Text(revHourUnknown && i == 3
                     ? "未知" : BaziCalc.jiaZiName(revIdx[i]))
                    .font(.system(size: 20, weight: .bold))
                    .foregroundColor(revHourUnknown && i == 3
                                     ? AppColors.textSecondary : AppColors.primary)
                Text("▾")
                    .font(.system(size: AppDimens.fontSmall))
                    .foregroundColor(AppColors.textSecondary)
            }
            .frame(maxWidth: .infinity)
            .frame(height: 70)
            .background(AppColors.background)
            .cornerRadius(AppDimens.radiusSm)
        }
        .buttonStyle(.plain)
    }

    private func pillarChooser(_ i: Int) -> some View {
        let list: [Int] = {
            if i == 1 { return BaziCalc.validMonths(yearIdx: revIdx[0]) }
            if i == 3 { return BaziCalc.validHours(dayIdx: revIdx[2]) }
            return Array(0..<60)
        }()
        let names = list.map { BaziCalc.jiaZiName($0) }
        let sel = list.firstIndex(of: ((revIdx[i] % 60) + 60) % 60) ?? 0
        return PickerSheet(title: pillarNames[i],
                           options: names,
                           selectedIndex: sel) { idx in
            var arr = revIdx
            arr[i] = list[idx]
            if i == 0 {
                arr[1] = BaziCalc.remapByZhi(
                    validList: BaziCalc.validMonths(yearIdx: arr[0]),
                    oldIdx: arr[1])
            }
            if i == 2 {
                arr[3] = BaziCalc.remapByZhi(
                    validList: BaziCalc.validHours(dayIdx: arr[2]),
                    oldIdx: arr[3])
            }
            revIdx = arr
            revPickerIdx = nil
        } onCancel: {
            revPickerIdx = nil
        }
    }

    // ── 数据刷新 ──────────────────────────────────────────
    private func refreshAll() {
        if inputMode == .lunar {
            refreshLunarMonths()
        }
        refreshValidDays()
    }

    private func refreshLunarMonths() {
        lunarMonths = SxwnlBridge.getLunarMonths(year: year)
    }

    private func refreshValidDays() {
        var days: [Int]
        if inputMode == .lunar {
            let n = BaziCalc.lunarDaysInMonth(year: year, month: month,
                                              isLeap: isLeap, isSpec: isSpec)
            days = Array(1...max(n, 1))
        } else {
            let d = SxwnlBridge.getSolarMonthValidDays(year: year, month: month)
            days = d.isEmpty
                ? Array(1...BaziCalc.solarDaysInMonth(year: year, month: month))
                : d
        }
        validDays = days
        if !days.contains(day) {
            var pick = days.first ?? 1
            for d in days where d <= day { pick = d }
            day = pick
        }
    }

    private func switchTab(_ i: Int) {
        guard let m = BirthInputMode(rawValue: i) else { return }
        inputMode = m
        if m == .lunar {
            refreshLunarMonths()
        }
        if m != .reverse {
            refreshValidDays()
        }
    }

    private func onYearInput() {
        if let y = BaziCalc.parseYear(yearText) {
            year = y
            if inputMode == .lunar { refreshLunarMonths() }
            refreshValidDays()
        }
    }

    private func confirm() {
        var sel = DateSelection()
        sel.inputMode = inputMode
        sel.year = year; sel.month = month; sel.day = day
        sel.hour = hour; sel.minute = minute
        sel.isLeap = isLeap; sel.isSpec = isSpec
        onConfirm(sel)
        isPresented = false
    }

    private func doReverse() {
        let sy = BaziCalc.parseYear(revStartText) ?? Int.min
        let ey = BaziCalc.parseYear(revEndText) ?? Int.min
        if sy == Int.min || ey == Int.min || sy > ey {
            revResults = []; revSearched = true; return
        }
        var g = [Int](), z = [Int]()
        for i in 0..<4 {
            let idx = ((revIdx[i] % 60) + 60) % 60
            g.append(idx % 10)
            z.append(idx % 12)
        }
        let hg = revHourUnknown ? -1 : g[3]
        let hz = revHourUnknown ? -1 : z[3]
        revResults = SxwnlBridge.baziReverse(
            yg: g[0], yz: z[0], mg: g[1], mz: z[1],
            dg: g[2], dz: z[2], hg: hg, hz: hz,
            startYear: sy, endYear: ey)
        revSearched = true
    }

    private func applyReverse(_ item: ReverseItem) {
        var sel = DateSelection()
        sel.inputMode = .solar
        sel.year = item.year; sel.month = item.month; sel.day = item.day
        sel.hour = item.hour >= 0 ? item.hour : 12
        sel.minute = 0
        onConfirm(sel)
        isPresented = false
    }

    private func formatReverseItem(_ item: ReverseItem) -> String {
        var s = "\(BaziCalc.formatYear(item.year))\(item.month)月\(item.day)日"
        if item.hour >= 0 { s += " \(item.hour)时" }
        return s
    }
}

// MARK: - 通用列表选项 sheet
private struct PickerSheet: View {
    let title: String
    let options: [String]
    let selectedIndex: Int
    let onSelect: (Int) -> Void
    let onCancel: () -> Void

    @State private var pick: Int

    init(title: String, options: [String],
         selectedIndex: Int,
         onSelect: @escaping (Int) -> Void,
         onCancel: @escaping () -> Void) {
        self.title = title
        self.options = options
        self.selectedIndex = selectedIndex
        self.onSelect = onSelect
        self.onCancel = onCancel
        _pick = State(initialValue: selectedIndex)
    }

    var body: some View {
        VStack(spacing: 0) {
            HStack {
                Button("取消") { onCancel() }
                    .foregroundColor(AppColors.textSecondary)
                Spacer()
                Text(title)
                    .font(.system(size: AppDimens.fontBody, weight: .bold))
                Spacer()
                Button("确定") { onSelect(pick) }
                    .foregroundColor(AppColors.accent)
            }
            .padding(AppDimens.paddingLg)
            Divider().background(AppColors.divider)
            Picker("", selection: $pick) {
                ForEach(options.indices, id: \.self) { i in
                    Text(options[i]).tag(i)
                }
            }
            .pickerStyle(.wheel)
            .frame(maxWidth: .infinity)
        }
        .background(AppColors.surface)
    }
}

// MARK: - 圆角形状 (顶部圆角)
private struct RoundedCorner: Shape {
    var radius: CGFloat = .infinity
    var corners: UIRectCorner = .allCorners

    func path(in rect: CGRect) -> Path {
        let path = UIBezierPath(
            roundedRect: rect,
            byRoundingCorners: corners,
            cornerRadii: CGSize(width: radius, height: radius))
        return Path(path.cgPath)
    }
}
