import SwiftUI

// ════════════════════════════════════════════════════════════════
//  BaziView — 与鸿蒙 BaziPage.ets 对齐的八字输入主页
//
//  顶部 banner + 输入卡片(姓名/性别/出生时间/经纬度/历法/真太阳时)
//  「排盘」→ 弹出 BaziResultView (sheet 展示简洁/专业版)
// ════════════════════════════════════════════════════════════════

struct BaziView: View {
    @State private var name: String = ""
    @State private var gender: Int = 0    // 0男 1女
    @State private var lifa: LiFaModern = .dingQi
    @State private var astEnabled: Bool = false

    // 出生地: 默认「北京市 / 天安门」, 手动开关关闭=省/区县两级菜单, 打开=手动输入经纬度
    @State private var manualGeo: Bool = false
    @State private var regionIdx: Int = 0
    @State private var cityIdx: Int = 0
    @State private var longitude: Double = Cities.default.longitude
    @State private var latitude: Double = Cities.default.latitude
    @State private var lonText: String = String(Cities.default.longitude)
    @State private var latText: String = String(Cities.default.latitude)

    @State private var sel: DateSelection = {
        var s = DateSelection()
        s.inputMode = .solar
        s.year = 1990; s.month = 1; s.day = 1
        s.hour = 12; s.minute = 0
        return s
    }()
    @State private var picked: Bool = false
    @State private var lunarMonths: [LunarMonth] = []

    @State private var showPicker = false
    @State private var showCityPicker = false
    @State private var resultArg: BaziResultArg?
    @State private var errorMsg: String?

    var body: some View {
        ScrollView {
            VStack(spacing: 0) {
                headerBanner
                inputCard
                    .padding(AppDimens.paddingMd)
            }
            .padding(.bottom, 32)
        }
        .background(AppColors.background)
        .onAppear { lunarMonths = SxwnlBridge.getLunarMonths(year: sel.year) }
        .sheet(isPresented: $showPicker) {
            BaziDateTimePicker(
                isPresented: $showPicker,
                initialSelection: sel,
                onConfirm: { newSel in
                    sel = newSel
                    lunarMonths = SxwnlBridge.getLunarMonths(year: newSel.year)
                    picked = true
                }
            )
            .presentationDetents([.medium, .large])
            .presentationDragIndicator(.visible)
        }
        .sheet(item: $resultArg) { a in
            BaziResultView(arg: a)
        }
        .sheet(isPresented: $showCityPicker) {
            CityPickerSheet(
                initRegionIdx: regionIdx,
                initCityIdx: cityIdx,
                onConfirm: { r, c, city in
                    regionIdx = r
                    cityIdx = c
                    longitude = city.longitude
                    latitude  = city.latitude
                    lonText   = String(city.longitude)
                    latText   = String(city.latitude)
                }
            )
            .presentationDetents([.medium, .large])
            .presentationDragIndicator(.visible)
        }
        .alert("排盘失败",
               isPresented: Binding(get: { errorMsg != nil },
                                    set: { if !$0 { errorMsg = nil } })) {
            Button("好") { errorMsg = nil }
        } message: {
            Text(errorMsg ?? "")
        }
    }

    // ── 顶部 banner ────────────────────────────────────────
    private var headerBanner: some View {
        VStack {
            Text("☯ 八字排盘")
                .font(.system(size: AppDimens.fontTitle, weight: .bold))
                .foregroundColor(AppColors.onPrimary)
        }
        .frame(maxWidth: .infinity)
        .padding(.top, AppDimens.paddingXl)
        .padding(.bottom, AppDimens.paddingLg)
        .background(
            LinearGradient(
                colors: [AppColors.gradientStart, AppColors.gradientEnd],
                startPoint: .topLeading, endPoint: .bottomTrailing)
        )
    }

    // ── 输入卡片 ────────────────────────────────────────────
    private var inputCard: some View {
        VStack(spacing: AppDimens.paddingMd) {
            // 姓名 + 性别
            HStack {
                TextField("姓名(选填)", text: $name)
                    .textFieldStyle(.plain)
                    .font(.system(size: AppDimens.fontBody))
                    .padding(.horizontal, 12)
                    .frame(height: 44)
                    .background(AppColors.background)
                    .cornerRadius(AppDimens.radiusSm)
                HStack(spacing: 0) {
                    genderChip("男", 0)
                    genderChip("女", 1)
                }
                .padding(.leading, 10)
                .overlay(RoundedRectangle(cornerRadius: AppDimens.radiusSm)
                    .stroke(AppColors.divider, lineWidth: 1)
                    .padding(.leading, 10))
                .clipShape(RoundedRectangle(cornerRadius: AppDimens.radiusSm)
                    .offset(x: 10))
            }

            // 出生时间(整条) — 点击弹出选择器
            Button { showPicker = true } label: {
                HStack {
                    Text(recordText())
                        .font(.system(size: AppDimens.fontSubtitle))
                        .foregroundColor(picked ? AppColors.onSurface : AppColors.textSecondary)
                    Spacer()
                    Text("▾")
                        .font(.system(size: AppDimens.fontBody))
                        .foregroundColor(AppColors.textSecondary)
                }
                .padding(.horizontal, 14)
                .frame(height: 48)
                .background(AppColors.background)
                .cornerRadius(AppDimens.radiusSm)
            }
            .buttonStyle(.plain)

            // 出生地: 省/区县级联(默认) 或 手动输入经纬度(开关切换)
            birthPlaceSection

            // 历法
            HStack(spacing: 0) {
                Text("历法")
                    .font(.system(size: AppDimens.fontCaption))
                    .foregroundColor(AppColors.textSecondary)
                    .frame(width: 48, alignment: .leading)
                HStack(spacing: 6) {
                    ForEach(LiFaModern.allCases) { opt in
                        lifaChip(opt)
                    }
                }
                .frame(maxWidth: .infinity)
            }

            // 真太阳时
            HStack {
                Toggle(isOn: $astEnabled) { EmptyView() }
                    .labelsHidden()
                    .toggleStyle(CheckboxToggleStyle())
                Text("使用真太阳时")
                    .font(.system(size: AppDimens.fontBody))
                    .foregroundColor(AppColors.onSurface)
                    .padding(.leading, 6)
                Spacer()
            }
            .padding(.bottom, AppDimens.paddingSm)

            Button { paiPan() } label: {
                Text("排  盘")
                    .font(.system(size: AppDimens.fontSubtitle, weight: .bold))
                    .foregroundColor(AppColors.primaryDark)
                    .frame(maxWidth: .infinity)
                    .frame(height: 48)
                    .background(AppColors.accent)
                    .cornerRadius(AppDimens.radiusLg)
            }
            .buttonStyle(.plain)
        }
        .padding(AppDimens.paddingLg)
        .background(AppColors.surface)
        .cornerRadius(AppDimens.radiusLg)
        .shadow(color: Color.black.opacity(0.08), radius: 6, y: 2)
    }

    private func genderChip(_ label: String, _ value: Int) -> some View {
        Button { gender = value } label: {
            Text(label)
                .font(.system(size: AppDimens.fontBody,
                              weight: gender == value ? .bold : .regular))
                .foregroundColor(gender == value
                                 ? AppColors.onPrimary : AppColors.textSecondary)
                .frame(width: 48, height: 44)
                .background(gender == value ? AppColors.primary : AppColors.surface)
        }
        .buttonStyle(.plain)
    }

    private func lifaChip(_ opt: LiFaModern) -> some View {
        Button { lifa = opt } label: {
            Text(opt.label)
                .font(.system(size: AppDimens.fontCaption,
                              weight: lifa == opt ? .bold : .regular))
                .foregroundColor(lifa == opt
                                 ? AppColors.onPrimary : AppColors.textSecondary)
                .frame(maxWidth: .infinity)
                .frame(height: 34)
                .background(lifa == opt ? AppColors.primary : AppColors.background)
                .cornerRadius(AppDimens.radiusSm)
        }
        .buttonStyle(.plain)
    }

    private func recordText() -> String {
        if !picked { return "请选择出生时间" }
        return BaziCalc.formatRecord(sel, lunarMonths: lunarMonths)
    }

    // ── 出生地选择 ──────────────────────────────────────────
    //   manualGeo=false: 两级 Menu [省/直辖市] → [区/县] (自动填入经纬度)
    //   manualGeo=true:  直接输入经纬度
    //   时间一律按北京时间处理, 经纬度仅用于真太阳时方程修正
    private var birthPlaceSection: some View {
        VStack(alignment: .leading, spacing: 6) {
            HStack {
                Text("出生地")
                    .font(.system(size: AppDimens.fontCaption))
                    .foregroundColor(AppColors.textSecondary)
                    .frame(width: 48, alignment: .leading)
                Text(birthPlaceSummary)
                    .font(.system(size: AppDimens.fontSmall))
                    .foregroundColor(AppColors.textSecondary)
                    .lineLimit(1)
                Spacer(minLength: 4)
                Text("手动")
                    .font(.system(size: AppDimens.fontSmall))
                    .foregroundColor(AppColors.textSecondary)
                Toggle("", isOn: $manualGeo)
                    .labelsHidden()
                    .toggleStyle(.switch)
                    .tint(AppColors.primary)
                    .scaleEffect(0.8)
            }

            if manualGeo {
                HStack(spacing: 4) {
                    Text("经度")
                        .font(.system(size: AppDimens.fontSmall))
                        .foregroundColor(AppColors.textSecondary)
                    TextField("", text: $lonText)
                        .textFieldStyle(.plain)
                        .keyboardType(.numbersAndPunctuation)
                        .font(.system(size: AppDimens.fontCaption))
                        .padding(.horizontal, 6)
                        .frame(maxWidth: .infinity)
                        .frame(height: 38)
                        .background(AppColors.background)
                        .cornerRadius(AppDimens.radiusSm)
                        .onChange(of: lonText) { v in
                            if let n = Double(v) { longitude = n }
                        }
                    Text("纬度")
                        .font(.system(size: AppDimens.fontSmall))
                        .foregroundColor(AppColors.textSecondary)
                        .padding(.leading, 4)
                    TextField("", text: $latText)
                        .textFieldStyle(.plain)
                        .keyboardType(.numbersAndPunctuation)
                        .font(.system(size: AppDimens.fontCaption))
                        .padding(.horizontal, 6)
                        .frame(maxWidth: .infinity)
                        .frame(height: 38)
                        .background(AppColors.background)
                        .cornerRadius(AppDimens.radiusSm)
                        .onChange(of: latText) { v in
                            if let n = Double(v) { latitude = n }
                        }
                }
            } else {
                // 触发器: 点击整行打开两栏弹窗 CityPickerSheet
                Button {
                    showCityPicker = true
                } label: {
                    HStack {
                        Text(regionCityLabel)
                            .font(.system(size: AppDimens.fontCaption))
                            .foregroundColor(AppColors.onSurface)
                            .lineLimit(1)
                        Spacer()
                        Text("▾")
                            .font(.system(size: AppDimens.fontSmall))
                            .foregroundColor(AppColors.textSecondary)
                    }
                    .padding(.horizontal, 14)
                    .frame(maxWidth: .infinity)
                    .frame(height: 38)
                    .background(AppColors.background)
                    .cornerRadius(AppDimens.radiusSm)
                }
                .buttonStyle(.plain)
            }
        }
    }

    private var birthPlaceSummary: String {
        let ew = longitude >= 0 ? "E" : "W"
        let ns = latitude  >= 0 ? "N" : "S"
        let lon = String(format: "%.2f", abs(longitude))
        let lat = String(format: "%.2f", abs(latitude))
        if manualGeo {
            return "\(lon)°\(ew), \(lat)°\(ns)"
        }
        let r = Cities.regions[regionIdx]
        let c = r.cities[cityIdx]
        return "\(r.name) · \(c.name) · \(lon)°\(ew), \(lat)°\(ns)"
    }

    private var regionCityLabel: String {
        let r = Cities.regions[regionIdx]
        let c = r.cities[cityIdx]
        return "\(r.name)  ▸  \(c.name)"
    }

    // ── 排盘 ────────────────────────────────────────────────
    private func paiPan() {
        // 八字时间一律按北京时间处理, 经纬度仅用于真太阳时方程修正
        let result = SxwnlBridge.calcBazi(
            name: name.isEmpty ? "匿名" : name,
            gender: gender == 1,
            year: sel.year, month: sel.month, day: sel.day,
            hour: sel.hour, minute: sel.minute, second: 0,
            isLunar: sel.inputMode == .lunar,
            isLeap: sel.inputMode == .lunar ? sel.isLeap : false,
            isSpec: sel.inputMode == .lunar ? sel.isSpec : false,
            astEnabled: astEnabled,
            longitude: longitude, latitude: latitude,
            lifa: lifa.rawValue
        )
        guard let r = result else {
            errorMsg = "排盘计算失败, 请检查输入"
            return
        }
        resultArg = BaziResultArg(
            result: r, birthYear: sel.year,
            astEnabled: astEnabled,
            longitude: longitude, latitude: latitude,
            lifaLabel: lifa.label
        )
    }
}

// 让 BaziResultArg 可作为 sheet item
extension BaziResultArg: Identifiable {
    var id: String {
        "\(result.userName)-\(birthYear)-\(astEnabled)-\(longitude)-\(latitude)-\(result.solarBirth)"
    }
}

// 自定义 Checkbox Toggle (鸿蒙 Checkbox 风格)
private struct CheckboxToggleStyle: ToggleStyle {
    func makeBody(configuration: Configuration) -> some View {
        Button {
            configuration.isOn.toggle()
        } label: {
            ZStack {
                RoundedRectangle(cornerRadius: 4)
                    .stroke(configuration.isOn
                            ? AppColors.primary : AppColors.divider,
                            lineWidth: 1.5)
                    .background(
                        RoundedRectangle(cornerRadius: 4)
                            .fill(configuration.isOn ? AppColors.primary : .clear)
                    )
                    .frame(width: 20, height: 20)
                if configuration.isOn {
                    Image(systemName: "checkmark")
                        .font(.system(size: 12, weight: .bold))
                        .foregroundColor(.white)
                }
            }
        }
        .buttonStyle(.plain)
    }
}
