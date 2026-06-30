import SwiftUI
import UIKit
import CoreLocation

// ════════════════════════════════════════════════════════════════
//  EclipseComponents — 日月食可视化组件集
//  与 Android `ui/components/*` 端 1:1 对齐
//
//  组件:
//    SolarPathMapCanvas      日食食带世界地图
//    SolarLocalDiscCanvas    本地日食 日/月圆盘动画
//    LunarDiscCanvas         月食 月+本影+半影 动画
//    EclipseTimelineBar      关键事件时间轴 + scrubber
//    PlaybackController      播放/暂停/速度
//    LocationPicker          预设城市 + 自定义观测点
//    GradientBadge           渐变徽章
// ════════════════════════════════════════════════════════════════

// MARK: - 通用模型

struct TimelineMark: Identifiable, Hashable {
    let id = UUID()
    let label: String
    let jd: Double
    let color: Color
}

// MARK: - GradientBadge

struct GradientBadge: View {
    let text: String
    let start: Color
    let end: Color
    var icon: String? = nil

    var body: some View {
        HStack(spacing: 4) {
            if let icon = icon { Text(icon).font(.system(size: 11)) }
            Text(text)
                .font(.system(size: 11, weight: .bold))
                .foregroundColor(.white)
        }
        .padding(.horizontal, 10).padding(.vertical, 4)
        .background(
            LinearGradient(colors: [start, end],
                           startPoint: .topLeading, endPoint: .bottomTrailing)
        )
        .clipShape(Capsule())
    }
}

// MARK: - EclipseTimelineBar

struct EclipseTimelineBar: View {
    let jdBegin: Double
    let jdEnd: Double
    @Binding var currentJd: Double
    let marks: [TimelineMark]

    var body: some View {
        GeometryReader { geo in
            ZStack(alignment: .leading) {
                Capsule().fill(AppColors.divider).frame(height: 6)

                // 进度填充
                Capsule().fill(LinearGradient(colors: [AppColors.gradientStart,
                                                       AppColors.gradientEnd],
                                              startPoint: .leading, endPoint: .trailing))
                    .frame(width: progressW(geo.size.width), height: 6)

                // 事件标记
                ForEach(marks) { m in
                    let x = mapX(m.jd, width: geo.size.width)
                    VStack(spacing: 2) {
                        Circle().fill(m.color).frame(width: 8, height: 8)
                        Text(m.label)
                            .font(.system(size: 9, weight: .semibold))
                            .foregroundColor(AppColors.onSurface)
                            .fixedSize()
                    }
                    .offset(x: x - 14, y: -22)
                }

                // 拖拽指示器
                Circle()
                    .fill(AppColors.onPrimary)
                    .overlay(Circle().stroke(AppColors.primary, lineWidth: 2))
                    .frame(width: 14, height: 14)
                    .offset(x: progressW(geo.size.width) - 7, y: 0)
                    .gesture(
                        DragGesture()
                            .onChanged { value in
                                let frac = max(0, min(1, value.location.x / geo.size.width))
                                currentJd = jdBegin + (jdEnd - jdBegin) * Double(frac)
                            }
                    )
            }
            .contentShape(Rectangle())
            .gesture(
                DragGesture(minimumDistance: 0)
                    .onChanged { value in
                        let frac = max(0, min(1, value.location.x / geo.size.width))
                        currentJd = jdBegin + (jdEnd - jdBegin) * Double(frac)
                    }
            )
        }
        .frame(height: 36)
    }

    private func progressW(_ width: CGFloat) -> CGFloat {
        guard jdEnd > jdBegin else { return 0 }
        let frac = (currentJd - jdBegin) / (jdEnd - jdBegin)
        return width * CGFloat(max(0, min(1, frac)))
    }

    private func mapX(_ jd: Double, width: CGFloat) -> CGFloat {
        guard jdEnd > jdBegin else { return 0 }
        let frac = (jd - jdBegin) / (jdEnd - jdBegin)
        return width * CGFloat(max(0, min(1, frac)))
    }
}

// MARK: - PlaybackController

struct PlaybackController: View {
    @Binding var isPlaying: Bool
    @Binding var currentJd: Double
    let jdBegin: Double
    let jdEnd: Double
    @Binding var speedMultiplier: Double

    /// 每秒推进多少 JD (×speed). 默认 60 倍速.
    private let baseStep: Double = 0.0001  // ≈ 8.6 秒/秒 (×speed)

    @State private var timer: Timer? = nil

    var body: some View {
        HStack(spacing: 14) {
            Button { reset() } label: {
                Image(systemName: "backward.end.fill")
                    .font(.system(size: 14))
                    .foregroundColor(AppColors.primary)
            }
            Button { togglePlay() } label: {
                Image(systemName: isPlaying ? "pause.circle.fill" : "play.circle.fill")
                    .font(.system(size: 32))
                    .foregroundColor(AppColors.primary)
            }
            Spacer()
            ForEach([0.5, 1.0, 2.0, 5.0], id: \.self) { sp in
                Button {
                    speedMultiplier = sp
                } label: {
                    Text(sp == 0.5 ? "0.5×" : (sp == 1 ? "1×" : (sp == 2 ? "2×" : "5×")))
                        .font(.system(size: 11, weight: speedMultiplier == sp ? .bold : .regular))
                        .foregroundColor(speedMultiplier == sp ? AppColors.onPrimary : AppColors.primary)
                        .padding(.horizontal, 8).padding(.vertical, 4)
                        .background(speedMultiplier == sp ? AppColors.primary : Color.clear)
                        .clipShape(Capsule())
                }
                .buttonStyle(.plain)
            }
        }
        .onDisappear { timer?.invalidate(); timer = nil; isPlaying = false }
    }

    private func togglePlay() {
        isPlaying.toggle()
        if isPlaying {
            timer = Timer.scheduledTimer(withTimeInterval: 1.0/30.0, repeats: true) { _ in
                let step = baseStep * speedMultiplier
                let next = currentJd + step
                if next >= jdEnd { currentJd = jdBegin; isPlaying = false; timer?.invalidate(); timer = nil }
                else { currentJd = next }
            }
        } else {
            timer?.invalidate(); timer = nil
        }
    }

    private func reset() {
        currentJd = jdBegin
        if isPlaying { togglePlay() }
    }
}

// MARK: - LocationPicker

struct LocationPicker: View {
    @Binding var selected: ObserverCity
    @State private var showCustom = false
    @State private var customLon: String = "116.4"
    @State private var customLat: String = "39.9"
    @State private var customTz: String  = "8"

    var body: some View {
        ScrollView(.horizontal, showsIndicators: false) {
            HStack(spacing: 8) {
                ForEach(Cities.preset) { c in
                    let active = selected.id == c.id
                    Button { selected = c } label: {
                        Text(c.name)
                            .font(.system(size: 12, weight: active ? .bold : .regular))
                            .foregroundColor(active ? AppColors.onPrimary : AppColors.primary)
                            .padding(.horizontal, 10).padding(.vertical, 6)
                            .background(active ? AppColors.primary : AppColors.background)
                            .clipShape(Capsule())
                            .overlay(Capsule().stroke(AppColors.primary.opacity(0.4),
                                                      lineWidth: active ? 0 : 1))
                    }
                    .buttonStyle(.plain)
                }
                Button { showCustom = true } label: {
                    Text("自定义…")
                        .font(.system(size: 12, weight: .semibold))
                        .foregroundColor(AppColors.primary)
                        .padding(.horizontal, 10).padding(.vertical, 6)
                        .overlay(Capsule().stroke(AppColors.primary, lineWidth: 1))
                }
                .buttonStyle(.plain)
            }
            .padding(.horizontal, 2)
        }
        .alert("自定义观测点", isPresented: $showCustom) {
            TextField("经度 (东+/西-)", text: $customLon).keyboardType(.numbersAndPunctuation)
            TextField("纬度 (北+/南-)", text: $customLat).keyboardType(.numbersAndPunctuation)
            TextField("时区 (h, 东+/西-)", text: $customTz).keyboardType(.numbersAndPunctuation)
            Button("确定") {
                if let lon = Double(customLon), let lat = Double(customLat), let tz = Double(customTz) {
                    selected = ObserverCity(name: "自定义",
                                             longitude: lon, latitude: lat, tzHours: tz)
                }
            }
            Button("取消", role: .cancel) {}
        } message: {
            Text("输入观测点的经纬度与时区")
        }
    }
}

// MARK: - SolarPathMapCanvas

struct SolarPathMapCanvas: View {
    let worldMapRad: [Double]      // 海岸 (ditu0 / ditu1)
    let bordersRad: [Double]?      // 国界 (ditu2), 可选
    let path: SolarEclipsePathData
    let currentJd: Double

    var body: some View {
        GeometryReader { geo in
            ZStack(alignment: .topLeading) {
                AppColors.mapOcean

                Canvas { ctx, size in
                    drawGraticule(ctx: ctx, size: size)
                    drawCoast(ctx: ctx, size: size, data: worldMapRad,
                              color: AppColors.mapLand, lineWidth: 0.6)
                    if let b = bordersRad, !b.isEmpty {
                        drawCoast(ctx: ctx, size: size, data: b,
                                  color: AppColors.mapBorder, lineWidth: 0.4)
                    }
                    drawBoundary(ctx: ctx, size: size, points: path.penumbraNorth,
                                 color: AppColors.pathPenumbra)
                    drawBoundary(ctx: ctx, size: size, points: path.penumbraSouth,
                                 color: AppColors.pathPenumbra)
                    drawBoundary(ctx: ctx, size: size, points: path.umbraNorth,
                                 color: AppColors.pathUmbra)
                    drawBoundary(ctx: ctx, size: size, points: path.umbraSouth,
                                 color: AppColors.pathUmbra)
                    drawCenter(ctx: ctx, size: size, points: path.centerLine)

                    drawKeyPoint(ctx: ctx, size: size, p: path.partialStart, label: "P₁",
                                 color: AppColors.pathPenumbra)
                    drawKeyPoint(ctx: ctx, size: size, p: path.centralStart, label: "C₁",
                                 color: AppColors.pathUmbra)
                    drawKeyPoint(ctx: ctx, size: size, p: path.greatestEclipse, label: "MAX",
                                 color: AppColors.pathCenter)
                    drawKeyPoint(ctx: ctx, size: size, p: path.centralEnd, label: "C₂",
                                 color: AppColors.pathUmbra)
                    drawKeyPoint(ctx: ctx, size: size, p: path.partialEnd, label: "P₄",
                                 color: AppColors.pathPenumbra)

                    if let pt = interpolateAtJd(path.centerLine, jd: currentJd) {
                        drawShadowMarker(ctx: ctx, size: size, point: pt)
                    }
                }
            }
        }
        .aspectRatio(2, contentMode: .fit)
        .clipShape(RoundedRectangle(cornerRadius: AppDimens.radiusMd))
    }

    // ─── 绘制工具 ───────────────────────────────────────────

    private func drawGraticule(ctx: GraphicsContext, size: CGSize) {
        for lon in stride(from: -180, through: 180, by: 30) {
            let x = size.width * CGFloat(EclipseUtil.projectLonX(Double(lon) * .pi / 180))
            var p = Path(); p.move(to: CGPoint(x: x, y: 0)); p.addLine(to: CGPoint(x: x, y: size.height))
            ctx.stroke(p, with: .color(AppColors.mapGrid), lineWidth: 0.4)
        }
        for lat in stride(from: -60, through: 60, by: 30) {
            let y = size.height * CGFloat(EclipseUtil.projectLatY(Double(lat) * .pi / 180))
            var p = Path(); p.move(to: CGPoint(x: 0, y: y)); p.addLine(to: CGPoint(x: size.width, y: y))
            ctx.stroke(p, with: .color(AppColors.mapGrid), lineWidth: 0.4)
        }
    }

    private func drawCoast(ctx: GraphicsContext, size: CGSize, data: [Double],
                           color: Color, lineWidth: CGFloat) {
        guard data.count >= 4 else { return }
        var p = Path()
        var isFirst = true
        var i = 0
        while i + 1 < data.count {
            let lon = data[i]; let lat = data[i+1]
            if lon >= 1e6 {
                isFirst = true; i += 2; continue
            }
            let x = size.width  * CGFloat(EclipseUtil.projectLonX(lon))
            let y = size.height * CGFloat(EclipseUtil.projectLatY(lat))
            if isFirst { p.move(to: CGPoint(x: x, y: y)); isFirst = false }
            else       { p.addLine(to: CGPoint(x: x, y: y)) }
            i += 2
        }
        ctx.stroke(p, with: .color(color), lineWidth: lineWidth)
    }

    private func drawBoundary(ctx: GraphicsContext, size: CGSize,
                              points: [GeoPoint], color: Color) {
        guard points.count > 1 else { return }
        var p = Path(); var first = true
        for pt in points {
            if pt.longitude > 99 || pt.latitude > 99 { first = true; continue }
            let x = size.width  * CGFloat(EclipseUtil.projectLonX(pt.longitude * .pi / 180))
            let y = size.height * CGFloat(EclipseUtil.projectLatY(pt.latitude  * .pi / 180))
            if first { p.move(to: CGPoint(x: x, y: y)); first = false }
            else     { p.addLine(to: CGPoint(x: x, y: y)) }
        }
        ctx.stroke(p, with: .color(color.opacity(0.85)), lineWidth: 1.4)
    }

    private func drawCenter(ctx: GraphicsContext, size: CGSize, points: [GeoPoint]) {
        guard points.count > 1 else { return }
        var p = Path(); var first = true
        for pt in points {
            if pt.longitude > 99 || pt.latitude > 99 { first = true; continue }
            let x = size.width  * CGFloat(EclipseUtil.projectLonX(pt.longitude * .pi / 180))
            let y = size.height * CGFloat(EclipseUtil.projectLatY(pt.latitude  * .pi / 180))
            if first { p.move(to: CGPoint(x: x, y: y)); first = false }
            else     { p.addLine(to: CGPoint(x: x, y: y)) }
        }
        ctx.stroke(p, with: .color(AppColors.pathCenter), lineWidth: 1.6)
    }

    private func drawKeyPoint(ctx: GraphicsContext, size: CGSize,
                              p pt: GeoPoint, label: String, color: Color) {
        if pt.longitude > 99 || pt.latitude > 99 { return }
        let x = size.width  * CGFloat(EclipseUtil.projectLonX(pt.longitude * .pi / 180))
        let y = size.height * CGFloat(EclipseUtil.projectLatY(pt.latitude  * .pi / 180))
        ctx.fill(Path(ellipseIn: CGRect(x: x - 4, y: y - 4, width: 8, height: 8)),
                 with: .color(color))
        ctx.stroke(Path(ellipseIn: CGRect(x: x - 6, y: y - 6, width: 12, height: 12)),
                   with: .color(color.opacity(0.5)), lineWidth: 1)
        ctx.draw(Text(label).font(.system(size: 9, weight: .bold))
                    .foregroundColor(.white), at: CGPoint(x: x, y: y - 14))
    }

    private func drawShadowMarker(ctx: GraphicsContext, size: CGSize, point: GeoPoint) {
        let x = size.width  * CGFloat(EclipseUtil.projectLonX(point.longitude * .pi / 180))
        let y = size.height * CGFloat(EclipseUtil.projectLatY(point.latitude  * .pi / 180))
        for r in [16.0, 12.0, 8.0] {
            ctx.fill(Path(ellipseIn: CGRect(x: x - r/2, y: y - r/2, width: r, height: r)),
                     with: .color(AppColors.pathMarker.opacity(0.20)))
        }
        ctx.fill(Path(ellipseIn: CGRect(x: x - 3, y: y - 3, width: 6, height: 6)),
                 with: .color(AppColors.pathMarker))
    }

    private func interpolateAtJd(_ pts: [GeoPoint], jd: Double) -> GeoPoint? {
        guard pts.count >= 2 else { return nil }
        if jd <= pts.first!.julianDay { return pts.first }
        if jd >= pts.last!.julianDay  { return pts.last  }
        for i in 0..<(pts.count - 1) {
            let a = pts[i], b = pts[i + 1]
            if jd >= a.julianDay && jd <= b.julianDay {
                let t = (jd - a.julianDay) / (b.julianDay - a.julianDay)
                return GeoPoint(coordinate: CLLocationCoordinate2D(
                    latitude: a.latitude + (b.latitude - a.latitude) * t,
                    longitude: a.longitude + (b.longitude - a.longitude) * t),
                                julianDay: jd)
            }
        }
        return nil
    }
}

// MARK: - SolarLocalDiscCanvas

struct SolarLocalDiscCanvas: View {
    let frames: [EclipseFrame]
    let currentJd: Double
    let city: ObserverCity

    var body: some View {
        GeometryReader { geo in
            ZStack {
                LinearGradient(colors: [AppColors.skyDeepNight, AppColors.skyMidNight,
                                        AppColors.skyDawn],
                               startPoint: .top, endPoint: .bottom)

                Canvas { ctx, size in
                    drawStars(ctx: ctx, size: size)

                    let frame = lerpFrame(currentJd) ?? frames.first
                    guard let f = frame, f.sunRadius > 0 else { return }
                    let cx = size.width  / 2
                    let cy = size.height / 2
                    let scale = min(size.width, size.height) * 0.35 / CGFloat(f.sunRadius)

                    // Halo
                    for r in [1.7, 1.4, 1.15] {
                        let rad = CGFloat(f.sunRadius) * scale * r
                        ctx.fill(Path(ellipseIn: CGRect(x: cx - rad, y: cy - rad,
                                                         width: rad * 2, height: rad * 2)),
                                 with: .color(AppColors.sunHalo.opacity(0.18)))
                    }

                    // Sun core gradient
                    let sunR = CGFloat(f.sunRadius) * scale
                    ctx.fill(Path(ellipseIn: CGRect(x: cx - sunR, y: cy - sunR,
                                                     width: sunR * 2, height: sunR * 2)),
                             with: .radialGradient(
                                Gradient(colors: [AppColors.sunCore, AppColors.sunGlow]),
                                center: CGPoint(x: cx, y: cy),
                                startRadius: 0, endRadius: sunR))

                    // Moon disc
                    let mx = cx + CGFloat(f.moonX) * scale
                    let my = cy - CGFloat(f.moonY) * scale  // Y 反向
                    let mr = CGFloat(f.moonRadius) * scale
                    ctx.fill(Path(ellipseIn: CGRect(x: mx - mr, y: my - mr,
                                                     width: mr * 2, height: mr * 2)),
                             with: .color(AppColors.moonDark))
                    ctx.stroke(Path(ellipseIn: CGRect(x: mx - mr, y: my - mr,
                                                      width: mr * 2, height: mr * 2)),
                               with: .color(AppColors.moonRim), lineWidth: 1)
                }

                VStack {
                    HStack {
                        DataChip(label: "食分",
                                 value: String(format: "%.3f",
                                               lerpFrame(currentJd)?.magnitude ?? 0))
                        DataChip(label: city.name,
                                 value: String(format: "%.1f°, %.1f°",
                                               city.longitude, city.latitude))
                        Spacer()
                        DataChip(label: EclipseUtil.tzLabel(city.tzHours),
                                 value: EclipseUtil.jdTdToLocal(currentJd,
                                                                tzHours: city.tzHours))
                    }
                    .padding(10)
                    Spacer()
                }
            }
        }
        .aspectRatio(1.6, contentMode: .fit)
        .clipShape(RoundedRectangle(cornerRadius: AppDimens.radiusMd))
    }

    private func drawStars(ctx: GraphicsContext, size: CGSize) {
        var rng = SystemRandomNumberGenerator()
        // 用固定种子的伪随机以减少抖动: 我们用 Hash
        for i in 0..<60 {
            let h = Double((i &* 2654435761) & 0xFFFF)
            let x = CGFloat((h * 0.0123).truncatingRemainder(dividingBy: 1)) * size.width
            let y = CGFloat((h * 0.04567).truncatingRemainder(dividingBy: 1)) * size.height
            let r = CGFloat((h * 0.0789).truncatingRemainder(dividingBy: 1)) * 1.2 + 0.4
            ctx.fill(Path(ellipseIn: CGRect(x: x, y: y, width: r, height: r)),
                     with: .color(.white.opacity(0.6)))
            _ = rng.next() // suppress unused warning
        }
    }

    private func lerpFrame(_ jd: Double) -> EclipseFrame? {
        guard frames.count >= 2 else { return frames.first }
        if jd <= frames.first!.julianDay { return frames.first }
        if jd >= frames.last!.julianDay  { return frames.last  }
        for i in 0..<(frames.count - 1) {
            let a = frames[i], b = frames[i + 1]
            if jd >= a.julianDay && jd <= b.julianDay {
                let t = (jd - a.julianDay) / (b.julianDay - a.julianDay)
                return EclipseFrame(
                    sunRadius: EclipseUtil.lerp(a.sunRadius, b.sunRadius, t),
                    moonX: EclipseUtil.lerp(a.moonX, b.moonX, t),
                    moonY: EclipseUtil.lerp(a.moonY, b.moonY, t),
                    moonRadius: EclipseUtil.lerp(a.moonRadius, b.moonRadius, t),
                    julianDay: jd,
                    magnitude: EclipseUtil.lerp(a.magnitude, b.magnitude, t))
            }
        }
        return frames.last
    }
}

// MARK: - LunarDiscCanvas

struct LunarDiscCanvas: View {
    let frames: [LunarEclipseFrame]
    let currentJd: Double

    var body: some View {
        GeometryReader { geo in
            ZStack {
                LinearGradient(colors: [AppColors.skyDeepNight, AppColors.skyMidNight],
                               startPoint: .top, endPoint: .bottom)

                Canvas { ctx, size in
                    let f = lerpFrame(currentJd) ?? frames.first
                    guard let f = f, f.penumbraRadius > 0 else { return }
                    let cx = size.width  / 2
                    let cy = size.height / 2
                    let scale = min(size.width, size.height) * 0.40 / CGFloat(f.penumbraRadius)

                    // Penumbra
                    let pr = CGFloat(f.penumbraRadius) * scale
                    ctx.fill(Path(ellipseIn: CGRect(x: cx - pr, y: cy - pr,
                                                     width: pr * 2, height: pr * 2)),
                             with: .radialGradient(
                                Gradient(colors: [AppColors.earthPenumbra.opacity(0.55),
                                                  AppColors.earthPenumbra.opacity(0.0)]),
                                center: CGPoint(x: cx, y: cy),
                                startRadius: pr * 0.6, endRadius: pr))

                    // Umbra
                    let ur = CGFloat(f.umbraRadius) * scale
                    ctx.fill(Path(ellipseIn: CGRect(x: cx - ur, y: cy - ur,
                                                     width: ur * 2, height: ur * 2)),
                             with: .color(AppColors.earthUmbra))

                    // Moon — colour shifts to blood moon as coverage grows
                    let mx = cx + CGFloat(f.moonX) * scale
                    let my = cy - CGFloat(f.moonY) * scale
                    let mr = CGFloat(f.moonRadius) * scale
                    let t = max(0, min(1, f.coverage))
                    let mixed = mix(AppColors.moonBright, AppColors.moonBlood, t)
                    ctx.fill(Path(ellipseIn: CGRect(x: mx - mr, y: my - mr,
                                                     width: mr * 2, height: mr * 2)),
                             with: .color(mixed))
                    ctx.stroke(Path(ellipseIn: CGRect(x: mx - mr, y: my - mr,
                                                      width: mr * 2, height: mr * 2)),
                               with: .color(AppColors.moonRim), lineWidth: 0.8)
                }

                VStack {
                    HStack {
                        DataChip(label: "覆盖度",
                                 value: String(format: "%.2f", lerpFrame(currentJd)?.coverage ?? 0))
                        Spacer()
                        DataChip(label: "时刻",
                                 value: EclipseUtil.jdToTime(currentJd))
                    }
                    .padding(10)
                    Spacer()
                }
            }
        }
        .aspectRatio(1.4, contentMode: .fit)
        .clipShape(RoundedRectangle(cornerRadius: AppDimens.radiusMd))
    }

    private func lerpFrame(_ jd: Double) -> LunarEclipseFrame? {
        guard frames.count >= 2 else { return frames.first }
        if jd <= frames.first!.julianDay { return frames.first }
        if jd >= frames.last!.julianDay  { return frames.last  }
        for i in 0..<(frames.count - 1) {
            let a = frames[i], b = frames[i + 1]
            if jd >= a.julianDay && jd <= b.julianDay {
                let t = (jd - a.julianDay) / (b.julianDay - a.julianDay)
                return LunarEclipseFrame(
                    moonX: EclipseUtil.lerp(a.moonX, b.moonX, t),
                    moonY: EclipseUtil.lerp(a.moonY, b.moonY, t),
                    moonRadius: EclipseUtil.lerp(a.moonRadius, b.moonRadius, t),
                    umbraRadius: EclipseUtil.lerp(a.umbraRadius, b.umbraRadius, t),
                    penumbraRadius: EclipseUtil.lerp(a.penumbraRadius, b.penumbraRadius, t),
                    julianDay: jd,
                    coverage: EclipseUtil.lerp(a.coverage, b.coverage, t))
            }
        }
        return frames.last
    }

    private func mix(_ a: Color, _ b: Color, _ t: Double) -> Color {
        // SwiftUI Color 不直接暴露 RGB; 用 UIColor 中转
        let ua = UIColor(a), ub = UIColor(b)
        var ar: CGFloat = 0, ag: CGFloat = 0, ab: CGFloat = 0, aa: CGFloat = 0
        var br: CGFloat = 0, bg: CGFloat = 0, bb: CGFloat = 0, ba: CGFloat = 0
        ua.getRed(&ar, green: &ag, blue: &ab, alpha: &aa)
        ub.getRed(&br, green: &bg, blue: &bb, alpha: &ba)
        return Color(red: Double(ar + (br - ar) * CGFloat(t)),
                     green: Double(ag + (bg - ag) * CGFloat(t)),
                     blue: Double(ab + (bb - ab) * CGFloat(t)))
    }
}

// MARK: - DataChip

struct DataChip: View {
    let label: String
    let value: String
    var body: some View {
        VStack(alignment: .leading, spacing: 1) {
            Text(label)
                .font(.system(size: 9, weight: .semibold))
                .foregroundColor(.white.opacity(0.7))
            Text(value)
                .font(.system(size: 12, weight: .bold).monospaced())
                .foregroundColor(.white)
        }
        .padding(.horizontal, 8).padding(.vertical, 4)
        .background(.ultraThinMaterial)
        .clipShape(RoundedRectangle(cornerRadius: 6))
    }
}
