import Foundation
import CoreLocation

// MARK: - Data Models

struct GeoPoint {
    let coordinate: CLLocationCoordinate2D
    let julianDay: Double

    var longitude: Double { coordinate.longitude }
    var latitude: Double { coordinate.latitude }
}

struct SolarEclipseInfo {
    let julianDay: Double
    let year: Int
    let month: Int
    let day: Int
    let hour: Int
    let minute: Int
    let type: String
    let typeName: String
    let gamma: Double
    let magnitude: Double
    let centerLongitude: Double
    let centerLatitude: Double
    let pathWidthKm: Double
    let durationSeconds: Double
    let saros: Int
    let sarosMember: Int
    let season: String

    var hasCenter: Bool { centerLongitude < 99 && centerLatitude < 99 }
}

struct SolarEclipsePathData {
    let type: String
    let typeName: String
    let gamma: Double
    let magnitude: Double
    let pathWidthKm: Double
    let durationSeconds: Double
    let maxEclipseJD: Double

    let centerLine: [GeoPoint]
    let umbraNorth: [GeoPoint]
    let umbraSouth: [GeoPoint]
    let penumbraNorth: [GeoPoint]
    let penumbraSouth: [GeoPoint]

    let partialStart: GeoPoint
    let centralStart: GeoPoint
    let greatestEclipse: GeoPoint
    let centralEnd: GeoPoint
    let partialEnd: GeoPoint
}

struct EclipseFrame {
    let sunRadius: Double
    let moonX: Double
    let moonY: Double
    let moonRadius: Double
    let julianDay: Double
    let magnitude: Double
}

struct LocalSolarEclipseData {
    let type: String
    let maxMagnitude: Double
    let moonSunRatio: Double

    let firstContact: Double   // JD
    let secondContact: Double
    let maximum: Double
    let thirdContact: Double
    let fourthContact: Double
    let sunrise: Double
    let sunset: Double

    let frames: [EclipseFrame]
}

struct LunarEclipseInfo {
    let julianDay: Double
    let year: Int
    let month: Int
    let day: Int
    let hour: Int
    let minute: Int
    let type: String
    let typeName: String
    let magnitude: Double
    let saros: Int
    let sarosMember: Int
    let season: String
}

struct LunarEclipseFrame {
    let moonX: Double
    let moonY: Double
    let moonRadius: Double
    let umbraRadius: Double
    let penumbraRadius: Double
    let julianDay: Double
    let coverage: Double
}

struct LunarEclipseDetailData {
    let type: String
    let typeName: String
    let maxMagnitude: Double

    let penumbraStart: Double  // JD (0 = N/A)
    let partialStart: Double
    let totalStart: Double
    let maximum: Double
    let totalEnd: Double
    let partialEnd: Double
    let penumbraEnd: Double

    let frames: [LunarEclipseFrame]
}

// MARK: - Bridge

final class SxwnlEclipseMapBridge {

    // MARK: - Solar Eclipse

    static func searchSolarEclipses(year: Int, month: Int, count: Int) -> [SolarEclipseInfo] {
        guard let list = sxwnl_solar_eclipse_search(Int32(year), Int32(month), Int32(count)) else {
            return []
        }
        defer { sxwnl_solar_eclipse_list_free(list) }

        var results: [SolarEclipseInfo] = []
        let n = Int(list.pointee.count)
        guard let items = list.pointee.items else { return [] }

        for i in 0..<n {
            let item = items[i]
            results.append(SolarEclipseInfo(
                julianDay: item.jd,
                year: Int(item.year),
                month: Int(item.month),
                day: Int(item.day),
                hour: Int(item.hour),
                minute: Int(item.minute),
                type: cStr(item.type),
                typeName: cStr(item.type_name),
                gamma: item.gamma,
                magnitude: item.magnitude,
                centerLongitude: item.center_lon,
                centerLatitude: item.center_lat,
                pathWidthKm: item.path_width,
                durationSeconds: item.duration,
                saros: Int(item.saros),
                sarosMember: Int(item.saros_member),
                season: cStr(item.season)
            ))
        }
        return results
    }

    static func getSolarEclipsePath(year: Int, month: Int, day: Int) -> SolarEclipsePathData? {
        guard let path = sxwnl_solar_eclipse_path(Int32(year), Int32(month), Int32(day)) else {
            return nil
        }
        defer { sxwnl_solar_eclipse_path_free(path) }

        let p = path.pointee
        return SolarEclipsePathData(
            type: cStr(p.type),
            typeName: cStr(p.type_name),
            gamma: p.gamma,
            magnitude: p.magnitude,
            pathWidthKm: p.path_width_km,
            durationSeconds: p.duration_seconds,
            maxEclipseJD: p.max_eclipse_jd,
            centerLine: convertPoints(p.center_line, count: Int(p.center_count)),
            umbraNorth: convertPoints(p.umbra_north, count: Int(p.umbra_north_count)),
            umbraSouth: convertPoints(p.umbra_south, count: Int(p.umbra_south_count)),
            penumbraNorth: convertPoints(p.penumbra_north, count: Int(p.penumbra_north_count)),
            penumbraSouth: convertPoints(p.penumbra_south, count: Int(p.penumbra_south_count)),
            partialStart: convertGeoPoint(p.partial_start),
            centralStart: convertGeoPoint(p.central_start),
            greatestEclipse: convertGeoPoint(p.greatest_eclipse),
            centralEnd: convertGeoPoint(p.central_end),
            partialEnd: convertGeoPoint(p.partial_end)
        )
    }

    static func getLocalSolarEclipse(
        year: Int, month: Int, day: Int,
        longitude: Double, latitude: Double,
        frameInterval: Int = 60
    ) -> LocalSolarEclipseData? {
        guard let data = sxwnl_solar_eclipse_local(
            Int32(year), Int32(month), Int32(day),
            longitude, latitude, Int32(frameInterval)
        ) else { return nil }
        defer { sxwnl_solar_eclipse_local_free(data) }

        let d = data.pointee
        var frames: [EclipseFrame] = []
        if let framesPtr = d.frames {
            for i in 0..<Int(d.count) {
                let f = framesPtr[i]
                frames.append(EclipseFrame(
                    sunRadius: f.sun_radius,
                    moonX: f.moon_x,
                    moonY: f.moon_y,
                    moonRadius: f.moon_radius,
                    julianDay: f.jd,
                    magnitude: f.magnitude
                ))
            }
        }

        return LocalSolarEclipseData(
            type: cStr(d.type),
            maxMagnitude: d.max_magnitude,
            moonSunRatio: d.moon_sun_ratio,
            firstContact: d.t_c1,
            secondContact: d.t_c2,
            maximum: d.t_max,
            thirdContact: d.t_c3,
            fourthContact: d.t_c4,
            sunrise: d.t_sunrise,
            sunset: d.t_sunset,
            frames: frames
        )
    }

    // MARK: - Lunar Eclipse

    static func searchLunarEclipses(year: Int, month: Int, count: Int) -> [LunarEclipseInfo] {
        guard let list = sxwnl_lunar_eclipse_search(Int32(year), Int32(month), Int32(count)) else {
            return []
        }
        defer { sxwnl_lunar_eclipse_list_free(list) }

        var results: [LunarEclipseInfo] = []
        let n = Int(list.pointee.count)
        guard let items = list.pointee.items else { return [] }

        for i in 0..<n {
            let item = items[i]
            results.append(LunarEclipseInfo(
                julianDay: item.jd,
                year: Int(item.year),
                month: Int(item.month),
                day: Int(item.day),
                hour: Int(item.hour),
                minute: Int(item.minute),
                type: cStr(item.type),
                typeName: cStr(item.type_name),
                magnitude: item.magnitude,
                saros: Int(item.saros),
                sarosMember: Int(item.saros_member),
                season: cStr(item.season)
            ))
        }
        return results
    }

    static func getLunarEclipseDetail(
        year: Int, month: Int, day: Int,
        frameInterval: Int = 120
    ) -> LunarEclipseDetailData? {
        guard let data = sxwnl_lunar_eclipse_detail(
            Int32(year), Int32(month), Int32(day), Int32(frameInterval)
        ) else { return nil }
        defer { sxwnl_lunar_eclipse_detail_free(data) }

        let d = data.pointee
        var frames: [LunarEclipseFrame] = []
        if let framesPtr = d.frames {
            for i in 0..<Int(d.count) {
                let f = framesPtr[i]
                frames.append(LunarEclipseFrame(
                    moonX: f.moon_x,
                    moonY: f.moon_y,
                    moonRadius: f.moon_radius,
                    umbraRadius: f.umbra_radius,
                    penumbraRadius: f.penumbra_radius,
                    julianDay: f.jd,
                    coverage: f.coverage
                ))
            }
        }

        return LunarEclipseDetailData(
            type: cStr(d.type),
            typeName: cStr(d.type_name),
            maxMagnitude: d.max_magnitude,
            penumbraStart: d.t_p1,
            partialStart: d.t_u1,
            totalStart: d.t_u2,
            maximum: d.t_max,
            totalEnd: d.t_u3,
            partialEnd: d.t_u4,
            penumbraEnd: d.t_p4,
            frames: frames
        )
    }

    // MARK: - World Map (海岸线轮廓)
    //
    // 返回 [Double], 经纬度(弧度) 交替: [lon0, lat0, lon1, lat1, ...]
    // 段间分隔点用 1e7 标记 (Move-To).
    //
    // 防御性: C 协议保证 got <= need, 但仍主动 clamp 以防 C 层升级后
    // 出现破坏性变更导致 removeSubrange 越界 crash.
    static func getWorldMapDitu0() -> [Double] {
        // 协议: 缓冲不足时返回 -need; 0 表示没数据.
        let probe = Int(sxwnl_world_map_get_ditu0(nil, 0))
        let need = probe < 0 ? -probe : probe
        guard need > 0 else { return [] }

        var buf = [Double](repeating: 0, count: need)
        var got = Int(sxwnl_world_map_get_ditu0(&buf, Int32(need)))
        guard got > 0 else { return [] }
        if got > need { got = need }
        if got < need { buf.removeSubrange(got..<need) }
        return buf
    }

    /// idx=1 → ditu1 (大图海岸); idx=2 → ditu2 (国界).
    /// 首次调用时自动从内建表加载.
    static func getWorldMapData(idx: Int) -> [Double] {
        // 协议: 缓冲不足时返回 -need; 0 表示没数据.
        let probe = Int(sxwnl_world_map_get_data(Int32(idx), nil, 0))
        let need = probe < 0 ? -probe : probe
        guard need > 0 else { return [] }
        var buf = [Double](repeating: 0, count: need)
        var got = Int(sxwnl_world_map_get_data(Int32(idx), &buf, Int32(need)))
        guard got > 0 else { return [] }
        if got > need { got = need }
        if got < need { buf.removeSubrange(got..<need) }
        return buf
    }

    // MARK: - Private Helpers

    private static func convertPoints(_ ptr: UnsafeMutablePointer<SxwnlGeoPoint>?, count: Int) -> [GeoPoint] {
        guard let ptr = ptr, count > 0 else { return [] }
        return (0..<count).map { i in
            convertGeoPoint(ptr[i])
        }
    }

    private static func convertGeoPoint(_ p: SxwnlGeoPoint) -> GeoPoint {
        GeoPoint(
            coordinate: CLLocationCoordinate2D(latitude: p.latitude, longitude: p.longitude),
            julianDay: p.jd
        )
    }

    private static func cStr<T>(_ tuple: T) -> String {
        withUnsafePointer(to: tuple) { ptr in
            ptr.withMemoryRebound(to: CChar.self, capacity: MemoryLayout<T>.size) {
                String(cString: $0)
            }
        }
    }
}
