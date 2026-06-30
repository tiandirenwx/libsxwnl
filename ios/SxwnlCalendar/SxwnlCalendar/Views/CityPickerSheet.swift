import SwiftUI

// ════════════════════════════════════════════════════════════════
//  CityPickerSheet — 出生地选择 (与鸿蒙 CityPickerDialog / Android CityPickerSheet 对齐)
//
//  sheet + HStack { 省 ScrollView | 区县 ScrollView } 联动
//  左栏选省 → 自动重置右栏到第 0 项, 顶部「确定」回调最终结果.
// ════════════════════════════════════════════════════════════════

struct CityPickerSheet: View {
    let initRegionIdx: Int
    let initCityIdx: Int
    let onConfirm: (Int, Int, ObserverCity) -> Void

    @Environment(\.dismiss) private var dismiss
    @State private var tmpRegionIdx: Int = 0
    @State private var tmpCityIdx: Int = 0

    var body: some View {
        VStack(spacing: 0) {
            // 标题栏
            HStack {
                Button("取消") { dismiss() }
                    .foregroundColor(AppColors.textSecondary)
                    .font(.system(size: AppDimens.fontBody))
                Spacer()
                Text("选择出生地")
                    .font(.system(size: AppDimens.fontSubtitle, weight: .bold))
                    .foregroundColor(AppColors.onSurface)
                Spacer()
                Button("确定") {
                    let r = Cities.regions[tmpRegionIdx]
                    let c = r.cities[tmpCityIdx]
                    onConfirm(tmpRegionIdx, tmpCityIdx, c)
                    dismiss()
                }
                .foregroundColor(AppColors.primary)
                .font(.system(size: AppDimens.fontBody, weight: .bold))
            }
            .padding(.horizontal, 16)
            .padding(.vertical, 12)

            Divider().background(AppColors.divider)

            // 当前选中预览
            Text(currentSummary)
                .font(.system(size: AppDimens.fontSmall))
                .foregroundColor(AppColors.textSecondary)
                .frame(maxWidth: .infinity, alignment: .leading)
                .lineLimit(1)
                .padding(.horizontal, 16)
                .padding(.vertical, 8)

            // 左省 + 右市
            HStack(spacing: 0) {
                // 左栏: 省/直辖市
                ScrollViewReader { proxy in
                    ScrollView {
                        LazyVStack(spacing: 0) {
                            ForEach(0..<Cities.regions.count, id: \.self) { idx in
                                regionRow(idx: idx).id(idx)
                            }
                        }
                    }
                    .onAppear {
                        if tmpRegionIdx > 0 {
                            DispatchQueue.main.async { proxy.scrollTo(tmpRegionIdx, anchor: .center) }
                        }
                    }
                }
                .frame(maxWidth: .infinity, maxHeight: .infinity)
                .background(AppColors.surface)

                // 右栏: 区/县 (随左栏切换 → ForEach key 含 tmpRegionIdx 触发重建)
                ScrollViewReader { proxy in
                    ScrollView {
                        LazyVStack(spacing: 0) {
                            let cities = Cities.regions[tmpRegionIdx].cities
                            ForEach(0..<cities.count, id: \.self) { idx in
                                cityRow(idx: idx, city: cities[idx]).id(idx)
                            }
                        }
                    }
                    .id("city_\(tmpRegionIdx)")        // 切省时强制重建 ScrollView, 自动回滚顶部
                    .onChange(of: tmpRegionIdx) { _, _ in
                        // 切省后右栏从顶开始
                        DispatchQueue.main.async { proxy.scrollTo(0, anchor: .top) }
                    }
                }
                .frame(maxWidth: .infinity, maxHeight: .infinity)
                .background(AppColors.background)
            }
            .frame(maxHeight: 360)
        }
        .background(AppColors.surface)
        .onAppear {
            tmpRegionIdx = max(0, min(initRegionIdx, Cities.regions.count - 1))
            tmpCityIdx = max(0, min(initCityIdx, Cities.regions[tmpRegionIdx].cities.count - 1))
        }
    }

    @ViewBuilder
    private func regionRow(idx: Int) -> some View {
        let r = Cities.regions[idx]
        let active = idx == tmpRegionIdx
        Button {
            if idx != tmpRegionIdx {
                tmpRegionIdx = idx
                tmpCityIdx = 0
            }
        } label: {
            HStack {
                Text(r.name)
                    .font(.system(size: AppDimens.fontCaption,
                                  weight: active ? .bold : .regular))
                    .foregroundColor(active ? AppColors.primary : AppColors.onSurface)
                Spacer()
            }
            .padding(.leading, 14).padding(.trailing, 8)
            .frame(height: 44)
            .frame(maxWidth: .infinity, alignment: .leading)
            .background(active ? AppColors.background : AppColors.surface)
        }
        .buttonStyle(.plain)
    }

    @ViewBuilder
    private func cityRow(idx: Int, city: ObserverCity) -> some View {
        let active = idx == tmpCityIdx
        Button {
            tmpCityIdx = idx
        } label: {
            HStack {
                Text(city.name)
                    .font(.system(size: AppDimens.fontCaption,
                                  weight: active ? .bold : .regular))
                    .foregroundColor(active ? AppColors.primary : AppColors.onSurface)
                Spacer()
            }
            .padding(.leading, 14).padding(.trailing, 8)
            .frame(height: 44)
            .frame(maxWidth: .infinity, alignment: .leading)
            .background(active ? AppColors.surface : AppColors.background)
        }
        .buttonStyle(.plain)
    }

    private var currentSummary: String {
        let r = Cities.regions[tmpRegionIdx]
        let c = r.cities[tmpCityIdx]
        let ew = c.longitude >= 0 ? "E" : "W"
        let ns = c.latitude  >= 0 ? "N" : "S"
        let lon = String(format: "%.2f", abs(c.longitude))
        let lat = String(format: "%.2f", abs(c.latitude))
        return "\(r.name) · \(c.name) · \(lon)°\(ew), \(lat)°\(ns)"
    }
}
