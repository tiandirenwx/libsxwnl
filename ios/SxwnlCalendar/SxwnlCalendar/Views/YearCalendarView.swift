import SwiftUI

// ════════════════════════════════════════════════════════════════
//  YearCalendarView — 与鸿蒙端 YearCalendarPage.ets 对齐的年历页面
//
//  按农历月聚合的一年数据 (对应 sxwnl/lunar.js 的 "时刻/干支" 模式):
//  每行展示一个农历月: [闰?][月名][大/小][朔日干支][朔日公历日期]
//                    该月内的节气: [日序][干支][节气名][日期]
// ════════════════════════════════════════════════════════════════

struct YearCalendarView: View {
    @State private var currentYear: Int
    @State private var yearInput: String
    @State private var months: [YearCalMonth] = []
    @State private var loading: Bool = false

    init() {
        let y = Calendar.current.dateComponents([.year], from: Date()).year ?? 2026
        _currentYear = State(initialValue: y)
        _yearInput = State(initialValue: YearUtil.astroYearToStr(y))
    }

    var body: some View {
        VStack(spacing: 0) {
            headerSection
            navSection
            summaryBar

            if loading {
                VStack {
                    ProgressView()
                        .progressViewStyle(CircularProgressViewStyle(tint: AppColors.primary))
                        .scaleEffect(1.2)
                    Text("正在计算 \(YearUtil.astroYearToStr(currentYear)) 年历...")
                        .font(.system(size: AppDimens.fontCaption))
                        .foregroundColor(AppColors.textSecondary)
                        .padding(.top, 8)
                }
                .frame(maxWidth: .infinity, maxHeight: .infinity)
            } else {
                ScrollView {
                    VStack(spacing: AppDimens.paddingSm) {
                        ForEach(months, id: \.monthIdx) { m in
                            monthRow(m)
                        }
                        Text("数据来源: 寿星万年历内核(C++) · 节气精确到秒")
                            .font(.system(size: AppDimens.fontSmall))
                            .foregroundColor(AppColors.textSecondary)
                            .frame(maxWidth: .infinity)
                            .padding(.vertical, AppDimens.paddingLg)
                    }
                    .padding(.horizontal, AppDimens.paddingMd)
                    .padding(.top, AppDimens.paddingMd)
                }
                .frame(maxHeight: .infinity)
            }
        }
        .background(AppColors.background)
        .onAppear { reload() }
    }

    // ── 顶部标题 ────────────────────────────────────────────
    private var headerSection: some View {
        HStack(alignment: .center) {
            VStack(alignment: .leading, spacing: 2) {
                Text("\(YearUtil.astroYearToStr(currentYear)) 年")
                    .font(.system(size: AppDimens.fontTitle, weight: .bold))
                    .foregroundColor(AppColors.onPrimary)
                Text("农历年视图 (按朔月排列)")
                    .font(.system(size: AppDimens.fontSmall))
                    .foregroundColor(AppColors.onPrimary.opacity(0.7))
            }
            Spacer()
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
            navIconButton("«") { shiftYear(-10) }
            navIconButton("‹") { shiftYear(-1) }

            TextField("YYYY/B212", text: $yearInput)
                .textFieldStyle(.plain)
                .multilineTextAlignment(.center)
                .font(.system(size: AppDimens.fontBody))
                .foregroundColor(AppColors.onPrimary)
                .frame(width: 80, height: 32)
                .background(AppColors.primaryLight.opacity(0.4))
                .cornerRadius(AppDimens.radiusSm)
                .submitLabel(.done)
                .onSubmit { applyYearInput() }

            Text("年")
                .font(.system(size: AppDimens.fontCaption))
                .foregroundColor(AppColors.onPrimary)
                .padding(.leading, 2)
                .padding(.trailing, 6)

            navIconButton("›") { shiftYear(1) }
            navIconButton("»") { shiftYear(10) }

            Spacer()

            Button("本年") { goCurYear() }
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

    // ── 摘要栏 ──────────────────────────────────────────────
    private var summaryBar: some View {
        HStack {
            Text("本年共 \(months.count) 个农历月")
                .font(.system(size: AppDimens.fontCaption))
                .foregroundColor(AppColors.textSecondary)
            if hasLeapMonth() {
                Text("  · 含 \(leapMonthName())")
                    .font(.system(size: AppDimens.fontCaption, weight: .medium))
                    .foregroundColor(AppColors.jieQi)
            }
            Spacer()
            Text("24 节气 / \(totalJieQi()) 项")
                .font(.system(size: AppDimens.fontSmall))
                .foregroundColor(AppColors.textSecondary)
        }
        .padding(.horizontal, AppDimens.paddingLg)
        .padding(.vertical, AppDimens.paddingSm)
        .background(AppColors.surface)
        .overlay(Rectangle()
            .frame(height: 0.5)
            .foregroundColor(AppColors.divider),
                 alignment: .bottom)
    }

    // ── 单月行 ──────────────────────────────────────────────
    private func monthRow(_ m: YearCalMonth) -> some View {
        VStack(alignment: .leading, spacing: 0) {
            HStack(alignment: .center, spacing: AppDimens.paddingMd) {
                Text(m.monthName)
                    .font(.system(size: AppDimens.fontSubtitle, weight: .bold))
                    .foregroundColor(m.isLeap ? AppColors.jieQi : AppColors.primary)
                    .frame(width: 70, alignment: .leading)

                VStack(alignment: .leading, spacing: 0) {
                    Text("朔 \(m.shuoGZ)")
                        .font(.system(size: AppDimens.fontBody, weight: .medium))
                        .foregroundColor(AppColors.onSurface)
                    Text(formatShuoDate(m))
                        .font(.system(size: AppDimens.fontSmall))
                        .foregroundColor(AppColors.textSecondary)
                }

                Spacer()

                Text(m.dayCount > 29 ? "大" : "小")
                    .font(.system(size: AppDimens.fontCaption))
                    .foregroundColor(AppColors.onPrimary)
                    .padding(.horizontal, 6)
                    .padding(.vertical, 2)
                    .background(m.dayCount > 29
                                ? AppColors.primary : AppColors.textSecondary)
                    .cornerRadius(4)
                    .padding(.trailing, 4)

                Text("\(m.dayCount) 天")
                    .font(.system(size: AppDimens.fontSmall))
                    .foregroundColor(AppColors.textSecondary)
            }
            .padding(.top, AppDimens.paddingMd)
            .padding(.bottom, AppDimens.paddingSm)

            if !m.jieQi.isEmpty {
                VStack(spacing: 0) {
                    ForEach(Array(m.jieQi.enumerated()), id: \.offset) { _, jq in
                        jieQiRow(jq, m)
                    }
                }
                .padding(.horizontal, AppDimens.paddingMd)
                .padding(.vertical, AppDimens.paddingSm)
                .background(AppColors.background)
                .cornerRadius(AppDimens.radiusSm)
                .padding(.top, AppDimens.paddingXs)
            }
        }
        .padding(.horizontal, AppDimens.paddingLg)
        .padding(.vertical, AppDimens.paddingSm)
        .frame(maxWidth: .infinity, alignment: .leading)
        .background(AppColors.surface)
        .cornerRadius(AppDimens.radiusMd)
    }

    private func jieQiRow(_ jq: YearCalJieQi, _ m: YearCalMonth) -> some View {
        HStack(spacing: 0) {
            Text(jq.dayName)
                .font(.system(size: AppDimens.fontCaption))
                .foregroundColor(AppColors.textSecondary)
                .frame(width: 36, alignment: .leading)
            Text(jq.gz)
                .font(.system(size: AppDimens.fontCaption, weight: .medium))
                .foregroundColor(AppColors.onSurface)
                .frame(width: 40, alignment: .leading)
            Text(jq.name)
                .font(.system(size: AppDimens.fontBody, weight: .medium))
                .foregroundColor(AppColors.jieQi)
                .frame(width: 48, alignment: .leading)
            Text(formatJieQiDate(jq, m))
                .font(.system(size: AppDimens.fontCaption))
                .foregroundColor(AppColors.textSecondary)
                .frame(maxWidth: .infinity, alignment: .leading)
        }
        .padding(.vertical, 2)
    }

    // ── 辅助方法 ────────────────────────────────────────────
    private func hasLeapMonth() -> Bool {
        months.contains(where: { $0.isLeap })
    }

    private func leapMonthName() -> String {
        months.first(where: { $0.isLeap })?.monthName ?? ""
    }

    private func totalJieQi() -> Int {
        months.reduce(0) { $0 + $1.jieQi.count }
    }

    private func formatShuoDate(_ m: YearCalMonth) -> String {
        let y = YearUtil.astroYearToStr(m.solarYear)
        return "\(y)-\(pad2(m.solarMonth))-\(pad2(m.solarDay))"
    }

    private func formatJieQiDate(_ jq: YearCalJieQi, _ m: YearCalMonth) -> String {
        let date = "\(pad2(jq.solarMonth))-\(pad2(jq.solarDay))"
        if !jq.time.isEmpty {
            return "\(date) \(jq.time)"
        }
        return "\(YearUtil.astroYearToStr(m.solarYear))-\(date)"
    }

    private func pad2(_ n: Int) -> String {
        n < 10 ? "0\(n)" : "\(n)"
    }

    // ── 数据加载 ────────────────────────────────────────────
    private func reload() {
        loading = true
        DispatchQueue.global(qos: .userInitiated).async {
            let data = SxwnlBridge.getYearCalendar(year: currentYear)
            DispatchQueue.main.async {
                months = data
                yearInput = YearUtil.astroYearToStr(currentYear)
                loading = false
            }
        }
    }

    private func shiftYear(_ delta: Int) {
        currentYear += delta
        reload()
    }

    private func goCurYear() {
        let y = Calendar.current.dateComponents([.year], from: Date()).year ?? 2026
        currentYear = y
        reload()
    }

    private func applyYearInput() {
        let y = YearUtil.yearStrToAstro(yearInput)
        if YearUtil.isAstroYearValid(y) {
            currentYear = y
            reload()
        } else {
            yearInput = YearUtil.astroYearToStr(currentYear)
        }
    }
}
