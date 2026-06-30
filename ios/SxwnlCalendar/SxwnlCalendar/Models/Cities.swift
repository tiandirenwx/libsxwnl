import Foundation

// ════════════════════════════════════════════════════════════════
//  Cities — 出生地经纬度目录 (薄包装 SxwnlBridge.geoXxx)
//
//  数据来源: libsxwnl C++ 层的 `GeoPosition` 单例 (src/geo.cpp 中的 `JWv`).
//
//  设计:
//  ─────────────────────────────────────────────────────────────
//  · 本文件**不再持有任何硬编码城市表**, 跨端唯一数据源 = C++ 核心.
//  · `regions / `default` / `preset` 使用 Swift `static let` 懒求值, 整个
//    App 共享同一份结果, 一次拉取 ~50ms.
//  · 任何 C++ 桥接异常 (无目录/空目录等) 都返回单条 FALLBACK, 防止 UI 崩溃.
//
//  与 Android `Cities.kt` / 鸿蒙 `Cities.ets` 字段结构对齐.
//  tzHours 字段供日食/RTS 等需要观测点时区的模块使用, 八字模块不读.
// ════════════════════════════════════════════════════════════════

struct ObserverCity: Identifiable, Hashable {
    let id = UUID()
    let name: String
    let longitude: Double
    let latitude: Double
    let tzHours: Double
}

struct Region: Identifiable, Hashable {
    let id = UUID()
    let name: String
    let cities: [ObserverCity]
}

enum Cities {

    // 兜底默认: native 数据完全不可用时使用 (避免 UI 空数组崩溃).
    private static let fallback = ObserverCity(
        name: "天安门", longitude: 116.3833, latitude: 39.9, tzHours: 8
    )

    static let regions: [Region] = loadRegions()

    static let `default`: ObserverCity = regions.first?.cities.first ?? fallback

    /// 旧版兼容: 扁平预设(供日食观测点快速选择条使用), 取首批 ~30 个.
    static let preset: [ObserverCity] = {
        guard let first = regions.first else { return [fallback] }
        return Array(first.cities.prefix(30))
    }()

    // MARK: - Loader

    private static func loadRegions() -> [Region] {
        let provinces = SxwnlBridge.geoListProvinces()
        guard !provinces.isEmpty else {
            return [Region(name: fallback.name, cities: [fallback])]
        }
        var out: [Region] = []
        out.reserveCapacity(provinces.count)
        for p in provinces {
            let raw = SxwnlBridge.geoListCities(province: p.province)
            guard !raw.isEmpty else { continue }
            let cities = raw.map { c in
                ObserverCity(name: c.district,
                             longitude: c.longitude,
                             latitude: c.latitude,
                             tzHours: c.timezone)
            }
            out.append(Region(name: p.province, cities: cities))
        }
        return out.isEmpty
            ? [Region(name: fallback.name, cities: [fallback])]
            : out
    }
}
