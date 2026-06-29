import { AppColors, AppDimens } from '../common/Constants'
import { EclipseUtil } from '../common/EclipseUtil'
import { Cities, ObserverCity } from '../common/Cities'
import {
  EclipseFrame, LunarEclipseFrame, EclipseGeoPoint, SolarEclipsePath
} from '../services/NativeBridge'

// ════════════════════════════════════════════════════════════════
//  EclipseComponents — 日月食可视化组件集 (ArkTS)
//
//  · TimelineMark / EclipseTimelineBar
//  · PlaybackController
//  · LocationStrip
//  · SolarPathMapCanvas
//  · SolarLocalDiscCanvas
//  · LunarDiscCanvas
//  · GradientBadge
//
//  与 Android `ui/components/*` / iOS `EclipseComponents.swift` 对齐
// ════════════════════════════════════════════════════════════════

// 改为 class: ArkTS 严格模式下 interface 不能用 {} as Interface 安全构造,
// 用 class + 构造函数避免 arkts-no-untyped-obj-literals 在 6.1+ 引发的运行时崩溃
export class TimelineMark {
  label: string
  jd: number
  color: string
  constructor(label: string, jd: number, color: string) {
    this.label = label
    this.jd = jd
    this.color = color
  }
}

@Component
export struct GradientBadge {
  text: string = ''
  color: string = '#1A237E'
  icon: string = ''

  build() {
    Row({ space: 4 }) {
      if (this.icon.length > 0) Text(this.icon).fontSize(11)
      Text(this.text).fontSize(11).fontWeight(FontWeight.Bold).fontColor('#FFFFFF')
    }
    .padding({ left: 10, right: 10, top: 4, bottom: 4 })
    .backgroundColor(this.color)
    .borderRadius(20)
  }
}

@Component
export struct EclipseTimelineBar {
  jdBegin: number = 0
  jdEnd: number = 0
  @Link currentJd: number
  marks: TimelineMark[] = []
  // 实测时间轴的真实宽度 (vp), 由 onAreaChange 写入, 用于把 jd → 像素偏移
  @State private barWidthVp: number = 0

  build() {
    // 改成 Column: 上半放标记, 下半放进度条; 不再依赖 position(%) 这种
    // 在 HarmonyOS 6.1 + bindSheet 重渲染时可能 native crash 的 API.
    Column() {
      // ─── 顶部 22vp: 事件标记 (Circle + 文字), 用 offset(vp) 定位 ───
      Stack({ alignContent: Alignment.TopStart }) {
        if (this.barWidthVp > 0) {
          ForEach(this.marks, (m: TimelineMark) => {
            Column({ space: 2 }) {
              Circle().width(8).height(8).fill(m.color)
              Text(m.label).fontSize(9).fontWeight(FontWeight.Bold).fontColor(AppColors.ON_SURFACE)
            }
            .width(20)
            .offset({ x: this.posVp(m.jd) - 10, y: 0 })
          }, (m: TimelineMark) => m.label + m.jd)
        }
      }.width('100%').height(22)

      // ─── 下半: 进度条本体, 实测宽度后用于 marks 偏移换算 ───
      Stack({ alignContent: Alignment.Start }) {
        Row().width('100%').height(6)
          .backgroundColor(AppColors.DIVIDER).borderRadius(3)
        Row().width(this.progress() + '%').height(6)
          .borderRadius(3)
          .linearGradient({ angle: 90, colors: [[AppColors.GRADIENT_START, 0], [AppColors.GRADIENT_END, 1]] })
      }
      .width('100%').height(6)
      .onAreaChange((_o: Area, n: Area) => {
        const w = Number(n.width) || 0
        if (w > 0 && Math.abs(w - this.barWidthVp) > 0.5) this.barWidthVp = w
      })
      .onTouch((e: TouchEvent) => {
        if (!e.touches || e.touches.length === 0) return
        let w = this.barWidthVp
        if (w <= 0) {
          try { w = Number(e.target.area.width) || 0 } catch (_) { w = 0 }
        }
        if (w <= 0) return
        const ratio = Math.max(0, Math.min(1, e.touches[0].x / w))
        this.currentJd = this.jdBegin + (this.jdEnd - this.jdBegin) * ratio
      })
    }.width('100%')
  }

  private progress(): number {
    if (this.jdEnd <= this.jdBegin) return 0
    return Math.max(0, Math.min(100,
      ((this.currentJd - this.jdBegin) / (this.jdEnd - this.jdBegin)) * 100))
  }

  private posVp(jd: number): number {
    if (this.jdEnd <= this.jdBegin || this.barWidthVp <= 0) return 0
    const ratio = Math.max(0, Math.min(1, (jd - this.jdBegin) / (this.jdEnd - this.jdBegin)))
    return ratio * this.barWidthVp
  }
}

@Component
export struct PlaybackController {
  jdBegin: number = 0
  jdEnd: number = 0
  @Link currentJd: number
  @State isPlaying: boolean = false
  @State speed: number = 1.0
  private intervalId: number = -1

  build() {
    Row({ space: 12 }) {
      Text('⏮')
        .fontSize(20).fontColor(AppColors.PRIMARY)
        .onClick(() => { this.currentJd = this.jdBegin })
      Text(this.isPlaying ? '❚❚' : '▶')
        .fontSize(28).fontColor(AppColors.PRIMARY).fontWeight(FontWeight.Bold)
        .onClick(() => this.togglePlay())
      Blank().layoutWeight(1)
      ForEach([0.5, 1, 2, 5], (sp: number) => {
        Text(sp === 0.5 ? '0.5×' : (sp === 1 ? '1×' : (sp === 2 ? '2×' : '5×')))
          .fontSize(11)
          .fontWeight(this.speed === sp ? FontWeight.Bold : FontWeight.Normal)
          .fontColor(this.speed === sp ? AppColors.ON_PRIMARY : AppColors.PRIMARY)
          .backgroundColor(this.speed === sp ? AppColors.PRIMARY : '#00000000')
          .padding({ left: 8, right: 8, top: 3, bottom: 3 })
          .borderRadius(12)
          .onClick(() => { this.speed = sp })
      })
    }
    .width('100%')
    .onDisAppear(() => {
      if (this.intervalId >= 0) clearInterval(this.intervalId)
    })
  }

  private togglePlay() {
    if (this.isPlaying) {
      if (this.intervalId >= 0) clearInterval(this.intervalId)
      this.intervalId = -1
      this.isPlaying = false
    } else {
      this.isPlaying = true
      this.intervalId = setInterval(() => {
        const step = 0.0001 * this.speed
        const next = this.currentJd + step
        if (next >= this.jdEnd) {
          this.currentJd = this.jdBegin
          this.isPlaying = false
          if (this.intervalId >= 0) clearInterval(this.intervalId)
          this.intervalId = -1
        } else {
          this.currentJd = next
        }
      }, 33)
    }
  }
}

@Component
export struct LocationStrip {
  @Link selected: ObserverCity
  @State private showCustom: boolean = false
  @State private customLon: string = '116.4'
  @State private customLat: string = '39.9'
  @State private customTz: string  = '8'

  build() {
    Column({ space: 6 }) {
      Scroll() {
        Row({ space: 8 }) {
          ForEach(Cities.PRESET, (c: ObserverCity) => {
            Text(c.name)
              .fontSize(12)
              .fontWeight(this.selected.name === c.name ? FontWeight.Bold : FontWeight.Normal)
              .fontColor(this.selected.name === c.name ? AppColors.ON_PRIMARY : AppColors.PRIMARY)
              .backgroundColor(this.selected.name === c.name ? AppColors.PRIMARY : AppColors.BACKGROUND)
              .padding({ left: 10, right: 10, top: 5, bottom: 5 })
              .borderRadius(16)
              .border({ width: 1, color: AppColors.PRIMARY, style: BorderStyle.Solid })
              .onClick(() => { this.selected = c })
          }, (c: ObserverCity) => c.name)
          Text(this.showCustom ? '收起' : '自定义…')
            .fontSize(12).fontColor(AppColors.PRIMARY)
            .backgroundColor(AppColors.BACKGROUND)
            .padding({ left: 10, right: 10, top: 5, bottom: 5 })
            .borderRadius(16)
            .border({ width: 1, color: AppColors.PRIMARY, style: BorderStyle.Dashed })
            .onClick(() => { this.showCustom = !this.showCustom })
        }.padding({ left: 4, right: 4 })
      }
      .width('100%').scrollable(ScrollDirection.Horizontal).scrollBar(BarState.Off)

      if (this.showCustom) {
        Row({ space: 6 }) {
          TextInput({ text: this.customLon, placeholder: '经度' })
            .layoutWeight(1).height(34).fontSize(12)
            .onChange((v: string) => { this.customLon = v })
          TextInput({ text: this.customLat, placeholder: '纬度' })
            .layoutWeight(1).height(34).fontSize(12)
            .onChange((v: string) => { this.customLat = v })
          TextInput({ text: this.customTz, placeholder: '时区' })
            .width(60).height(34).fontSize(12)
            .onChange((v: string) => { this.customTz = v })
          Text('应用')
            .fontSize(12).fontColor(AppColors.ON_PRIMARY)
            .backgroundColor(AppColors.PRIMARY)
            .padding({ left: 10, right: 10, top: 7, bottom: 7 })
            .borderRadius(8)
            .onClick(() => {
              const lon = Number(this.customLon)
              const lat = Number(this.customLat)
              const tz  = Number(this.customTz)
              if (!isNaN(lon) && !isNaN(lat) && !isNaN(tz)) {
                const c: ObserverCity = {
                  name: '自定义', longitude: lon, latitude: lat, tzHours: tz
                }
                this.selected = c
              }
            })
        }.width('100%').padding({ left: 4, right: 4 })
      }
    }
    .width('100%')
  }
}

// ─── SolarPathMapCanvas ──────────────────────────────────────

// 局部插值结果类 — 避免 ArkTS 严格模式禁止 {} as Interface 的运行时崩溃
class GeoMarker {
  longitude: number
  latitude: number
  jd: number
  constructor(longitude: number, latitude: number, jd: number) {
    this.longitude = longitude
    this.latitude = latitude
    this.jd = jd
  }
}

@Component
export struct SolarPathMapCanvas {
  // 用 number[] 代替 Float64Array - HarmonyOS ArkTS 严格模式不认
  // NAPI 返回的 Float64Array, 访问 .length 会抛
  // "This value does not have a [[TypedArrayName]] internal slot.".
  worldMapRad: number[] = []
  bordersRad: number[] = []
  path: SolarEclipsePath | null = null
  // dataVersion 由父组件每次数据更新时 +1, 触发 @Watch 后重绘
  @Prop @Watch('onDataChange') dataVersion: number = 0
  @Prop @Watch('onJdChange') currentJd: number = 0
  private settings: RenderingContextSettings = new RenderingContextSettings(true)
  private ctx: CanvasRenderingContext2D = new CanvasRenderingContext2D(this.settings)
  private ready: boolean = false

  build() {
    Stack() {
      Canvas(this.ctx)
        .width('100%').aspectRatio(2)
        .backgroundColor(AppColors.MAP_OCEAN)
        .borderRadius(AppDimens.BORDER_RADIUS_MD)
        .onReady(() => { this.ready = true; this.safeDraw() })
    }
  }

  onJdChange(): void { if (this.ready) this.safeDraw() }
  onDataChange(): void { if (this.ready) this.safeDraw() }

  // 整体 try/catch — 单帧绘制异常不应该杀进程
  private safeDraw(): void {
    try { this.draw() } catch (e) {
      console.error('[SolarPathMapCanvas.draw] ' + (e as Error).message)
    }
  }

  draw() {
    const w = this.ctx.width, h = this.ctx.height
    if (w <= 0 || h <= 0) return
    this.ctx.fillStyle = AppColors.MAP_OCEAN
    this.ctx.fillRect(0, 0, w, h)
    // 网格
    this.ctx.lineWidth = 0.4
    this.ctx.strokeStyle = AppColors.MAP_GRID
    for (let lon = -180; lon <= 180; lon += 30) {
      const x = w * EclipseUtil.projectLonX(lon * Math.PI / 180)
      this.ctx.beginPath(); this.ctx.moveTo(x, 0); this.ctx.lineTo(x, h); this.ctx.stroke()
    }
    for (let lat = -60; lat <= 60; lat += 30) {
      const y = h * EclipseUtil.projectLatY(lat * Math.PI / 180)
      this.ctx.beginPath(); this.ctx.moveTo(0, y); this.ctx.lineTo(w, y); this.ctx.stroke()
    }
    // 海岸
    this.drawCoast(this.worldMapRad, AppColors.MAP_LAND, 0.6)
    if (this.bordersRad && this.bordersRad.length > 0) {
      this.drawCoast(this.bordersRad, AppColors.MAP_BORDER, 0.4)
    }
    // 食带边界
    const p = this.path
    if (p) {
      this.drawPts(p.penumbraNorth, AppColors.PATH_PENUMBRA, 1.4)
      this.drawPts(p.penumbraSouth, AppColors.PATH_PENUMBRA, 1.4)
      this.drawPts(p.umbraNorth, AppColors.PATH_UMBRA, 1.4)
      this.drawPts(p.umbraSouth, AppColors.PATH_UMBRA, 1.4)
      this.drawPts(p.centerLine, AppColors.PATH_CENTER, 1.6)
      this.drawKey(p.partialStart, 'P₁', AppColors.PATH_PENUMBRA)
      this.drawKey(p.centralStart, 'C₁', AppColors.PATH_UMBRA)
      this.drawKey(p.greatestEclipse, 'MAX', AppColors.PATH_CENTER)
      this.drawKey(p.centralEnd, 'C₂', AppColors.PATH_UMBRA)
      this.drawKey(p.partialEnd, 'P₄', AppColors.PATH_PENUMBRA)
      const marker = this.interpolate(p.centerLine, this.currentJd)
      if (marker) {
        const mx = w * EclipseUtil.projectLonX(marker.longitude * Math.PI / 180)
        const my = h * EclipseUtil.projectLatY(marker.latitude * Math.PI / 180)
        this.ctx.fillStyle = AppColors.PATH_MARKER
        for (const r of [16, 12, 8]) {
          this.ctx.globalAlpha = 0.2
          this.ctx.beginPath(); this.ctx.arc(mx, my, r, 0, 2 * Math.PI); this.ctx.fill()
        }
        this.ctx.globalAlpha = 1
        this.ctx.beginPath(); this.ctx.arc(mx, my, 3, 0, 2 * Math.PI); this.ctx.fill()
      }
    }
  }

  private drawCoast(data: number[], color: string, lineWidth: number) {
    if (!data || data.length < 4) return
    const w = this.ctx.width, h = this.ctx.height
    this.ctx.lineWidth = lineWidth
    this.ctx.strokeStyle = color
    this.ctx.beginPath()
    let first = true
    for (let i = 0; i + 1 < data.length; i += 2) {
      const lon = data[i], lat = data[i + 1]
      if (lon >= 1e6) { first = true; continue }
      const x = w * EclipseUtil.projectLonX(lon)
      const y = h * EclipseUtil.projectLatY(lat)
      if (first) { this.ctx.moveTo(x, y); first = false }
      else this.ctx.lineTo(x, y)
    }
    this.ctx.stroke()
  }

  private drawPts(pts: EclipseGeoPoint[], color: string, lw: number) {
    if (!pts || pts.length < 2) return
    const w = this.ctx.width, h = this.ctx.height
    this.ctx.strokeStyle = color; this.ctx.lineWidth = lw
    this.ctx.beginPath()
    let first = true
    for (const p of pts) {
      if (p.longitude > 99 || p.latitude > 99) { first = true; continue }
      const x = w * EclipseUtil.projectLonX(p.longitude * Math.PI / 180)
      const y = h * EclipseUtil.projectLatY(p.latitude * Math.PI / 180)
      if (first) { this.ctx.moveTo(x, y); first = false }
      else this.ctx.lineTo(x, y)
    }
    this.ctx.stroke()
  }

  private drawKey(p: EclipseGeoPoint, label: string, color: string) {
    if (!p || p.longitude > 99 || p.latitude > 99) return
    const w = this.ctx.width, h = this.ctx.height
    const x = w * EclipseUtil.projectLonX(p.longitude * Math.PI / 180)
    const y = h * EclipseUtil.projectLatY(p.latitude * Math.PI / 180)
    this.ctx.fillStyle = color
    this.ctx.beginPath(); this.ctx.arc(x, y, 4, 0, 2 * Math.PI); this.ctx.fill()
    this.ctx.fillStyle = '#FFFFFF'
    this.ctx.font = 'bold 9px sans-serif'
    this.ctx.textAlign = 'center'
    this.ctx.fillText(label, x, y - 6)
  }

  // 注意: 返回值用 GeoMarker class 而不是 {...} as EclipseGeoPoint;
  // ArkTS 严格模式下 interface 的对象字面量断言会触发 arkts-no-untyped-obj-literals
  private interpolate(pts: EclipseGeoPoint[], jd: number): GeoMarker | null {
    if (!pts || pts.length < 2) return null
    if (jd <= pts[0].jd) return new GeoMarker(pts[0].longitude, pts[0].latitude, pts[0].jd)
    const last = pts[pts.length - 1]
    if (jd >= last.jd) return new GeoMarker(last.longitude, last.latitude, last.jd)
    for (let i = 0; i < pts.length - 1; i++) {
      const a = pts[i], b = pts[i + 1]
      if (jd >= a.jd && jd <= b.jd) {
        const t = (jd - a.jd) / (b.jd - a.jd)
        return new GeoMarker(
          a.longitude + (b.longitude - a.longitude) * t,
          a.latitude + (b.latitude - a.latitude) * t,
          jd
        )
      }
    }
    return null
  }
}

// ─── SolarLocalDiscCanvas ────────────────────────────────────

// 局部插值帧类 — 避免 ArkTS 对接口对象字面量断言的运行时崩溃
class SolarFrameLerp {
  sunRadius: number
  moonX: number
  moonY: number
  moonRadius: number
  jd: number
  magnitude: number
  constructor(sunRadius: number, moonX: number, moonY: number,
              moonRadius: number, jd: number, magnitude: number) {
    this.sunRadius = sunRadius
    this.moonX = moonX
    this.moonY = moonY
    this.moonRadius = moonRadius
    this.jd = jd
    this.magnitude = magnitude
  }
}

@Component
export struct SolarLocalDiscCanvas {
  frames: EclipseFrame[] = []
  @Prop @Watch('onJdChange') currentJd: number = 0
  @Prop @Watch('onJdChange') dataVersion: number = 0
  tzHours: number = 8
  private settings: RenderingContextSettings = new RenderingContextSettings(true)
  private ctx: CanvasRenderingContext2D = new CanvasRenderingContext2D(this.settings)
  private ready: boolean = false

  build() {
    Canvas(this.ctx)
      .width('100%').aspectRatio(1.6)
      .borderRadius(AppDimens.BORDER_RADIUS_MD)
      .linearGradient({ angle: 180, colors: [
        [AppColors.SKY_DEEP_NIGHT, 0], [AppColors.SKY_MID_NIGHT, 0.6], [AppColors.SKY_DAWN, 1]
      ]})
      .onReady(() => { this.ready = true; this.safeDraw() })
  }

  onJdChange(): void { if (this.ready) this.safeDraw() }

  private safeDraw(): void {
    try { this.draw() } catch (e) {
      console.error('[SolarLocalDiscCanvas.draw] ' + (e as Error).message)
    }
  }

  draw() {
    const w = this.ctx.width, h = this.ctx.height
    this.ctx.clearRect(0, 0, w, h)
    // 星空背景
    this.ctx.fillStyle = '#FFFFFF'
    this.ctx.globalAlpha = 0.6
    for (let i = 0; i < 60; i++) {
      const h2 = (i * 2654435761) & 0xFFFF
      const x = ((h2 * 0.0123) % 1) * w
      const y = ((h2 * 0.04567) % 1) * h
      const r = ((h2 * 0.0789) % 1) * 1.2 + 0.4
      this.ctx.beginPath(); this.ctx.arc(x, y, r, 0, 2 * Math.PI); this.ctx.fill()
    }
    this.ctx.globalAlpha = 1

    const f = this.lerpFrame(this.currentJd)
    if (!f || f.sunRadius <= 0) return
    const cx = w / 2, cy = h / 2
    const scale = Math.min(w, h) * 0.35 / f.sunRadius
    // 光晕
    this.ctx.fillStyle = AppColors.SUN_HALO
    this.ctx.globalAlpha = 0.18
    for (const r of [1.7, 1.4, 1.15]) {
      const rad = f.sunRadius * scale * r
      this.ctx.beginPath(); this.ctx.arc(cx, cy, rad, 0, 2 * Math.PI); this.ctx.fill()
    }
    this.ctx.globalAlpha = 1
    // 日盘渐变
    const sunR = f.sunRadius * scale
    const grad = this.ctx.createRadialGradient(cx, cy, 0, cx, cy, sunR)
    grad.addColorStop(0, AppColors.SUN_CORE)
    grad.addColorStop(1, AppColors.SUN_GLOW)
    this.ctx.fillStyle = grad
    this.ctx.beginPath(); this.ctx.arc(cx, cy, sunR, 0, 2 * Math.PI); this.ctx.fill()
    // 月盘
    const mx = cx + f.moonX * scale
    const my = cy - f.moonY * scale
    const mr = f.moonRadius * scale
    this.ctx.fillStyle = AppColors.MOON_DARK
    this.ctx.beginPath(); this.ctx.arc(mx, my, mr, 0, 2 * Math.PI); this.ctx.fill()
    this.ctx.strokeStyle = AppColors.MOON_RIM
    this.ctx.lineWidth = 1
    this.ctx.beginPath(); this.ctx.arc(mx, my, mr, 0, 2 * Math.PI); this.ctx.stroke()
  }

  private lerpFrame(jd: number): SolarFrameLerp | null {
    if (!this.frames || this.frames.length === 0) return null
    const wrap = (f: EclipseFrame): SolarFrameLerp =>
      new SolarFrameLerp(f.sunRadius, f.moonX, f.moonY, f.moonRadius, f.jd, f.magnitude)
    if (this.frames.length === 1) return wrap(this.frames[0])
    if (jd <= this.frames[0].jd) return wrap(this.frames[0])
    if (jd >= this.frames[this.frames.length - 1].jd) return wrap(this.frames[this.frames.length - 1])
    for (let i = 0; i < this.frames.length - 1; i++) {
      const a = this.frames[i], b = this.frames[i + 1]
      if (jd >= a.jd && jd <= b.jd) {
        const t = (jd - a.jd) / (b.jd - a.jd)
        return new SolarFrameLerp(
          EclipseUtil.lerp(a.sunRadius, b.sunRadius, t),
          EclipseUtil.lerp(a.moonX, b.moonX, t),
          EclipseUtil.lerp(a.moonY, b.moonY, t),
          EclipseUtil.lerp(a.moonRadius, b.moonRadius, t),
          jd,
          EclipseUtil.lerp(a.magnitude, b.magnitude, t)
        )
      }
    }
    return wrap(this.frames[this.frames.length - 1])
  }
}

// ─── LunarDiscCanvas ────────────────────────────────────────

class LunarFrameLerp {
  moonX: number
  moonY: number
  moonRadius: number
  umbraRadius: number
  penumbraRadius: number
  jd: number
  coverage: number
  constructor(moonX: number, moonY: number, moonRadius: number,
              umbraRadius: number, penumbraRadius: number, jd: number, coverage: number) {
    this.moonX = moonX
    this.moonY = moonY
    this.moonRadius = moonRadius
    this.umbraRadius = umbraRadius
    this.penumbraRadius = penumbraRadius
    this.jd = jd
    this.coverage = coverage
  }
}

@Component
export struct LunarDiscCanvas {
  frames: LunarEclipseFrame[] = []
  @Prop @Watch('onJdChange') currentJd: number = 0
  @Prop @Watch('onJdChange') dataVersion: number = 0
  private settings: RenderingContextSettings = new RenderingContextSettings(true)
  private ctx: CanvasRenderingContext2D = new CanvasRenderingContext2D(this.settings)
  private ready: boolean = false

  build() {
    Canvas(this.ctx)
      .width('100%').aspectRatio(1.4)
      .borderRadius(AppDimens.BORDER_RADIUS_MD)
      .linearGradient({ angle: 180, colors: [
        [AppColors.SKY_DEEP_NIGHT, 0], [AppColors.SKY_MID_NIGHT, 1]
      ]})
      .onReady(() => { this.ready = true; this.safeDraw() })
  }

  onJdChange(): void { if (this.ready) this.safeDraw() }

  private safeDraw(): void {
    try { this.draw() } catch (e) {
      console.error('[LunarDiscCanvas.draw] ' + (e as Error).message)
    }
  }

  draw() {
    const w = this.ctx.width, h = this.ctx.height
    this.ctx.clearRect(0, 0, w, h)
    const f = this.lerpFrame(this.currentJd)
    if (!f || f.penumbraRadius <= 0) return
    const cx = w / 2, cy = h / 2
    const scale = Math.min(w, h) * 0.4 / f.penumbraRadius

    // 半影
    const pr = f.penumbraRadius * scale
    const pg = this.ctx.createRadialGradient(cx, cy, pr * 0.6, cx, cy, pr)
    pg.addColorStop(0, '#1E2840AA')
    pg.addColorStop(1, '#1E284000')
    this.ctx.fillStyle = pg
    this.ctx.beginPath(); this.ctx.arc(cx, cy, pr, 0, 2 * Math.PI); this.ctx.fill()
    // 本影
    const ur = f.umbraRadius * scale
    this.ctx.fillStyle = AppColors.EARTH_UMBRA
    this.ctx.beginPath(); this.ctx.arc(cx, cy, ur, 0, 2 * Math.PI); this.ctx.fill()
    // 月亮 — 颜色按覆盖度从亮到血月
    const mx = cx + f.moonX * scale
    const my = cy - f.moonY * scale
    const mr = f.moonRadius * scale
    this.ctx.fillStyle = this.mix(AppColors.MOON_BRIGHT, AppColors.MOON_BLOOD, Math.max(0, Math.min(1, f.coverage)))
    this.ctx.beginPath(); this.ctx.arc(mx, my, mr, 0, 2 * Math.PI); this.ctx.fill()
    this.ctx.strokeStyle = AppColors.MOON_RIM; this.ctx.lineWidth = 0.8
    this.ctx.beginPath(); this.ctx.arc(mx, my, mr, 0, 2 * Math.PI); this.ctx.stroke()
  }

  private lerpFrame(jd: number): LunarFrameLerp | null {
    if (!this.frames || this.frames.length === 0) return null
    const wrap = (f: LunarEclipseFrame): LunarFrameLerp =>
      new LunarFrameLerp(f.moonX, f.moonY, f.moonRadius, f.umbraRadius, f.penumbraRadius, f.jd, f.coverage)
    if (this.frames.length === 1) return wrap(this.frames[0])
    if (jd <= this.frames[0].jd) return wrap(this.frames[0])
    if (jd >= this.frames[this.frames.length - 1].jd) return wrap(this.frames[this.frames.length - 1])
    for (let i = 0; i < this.frames.length - 1; i++) {
      const a = this.frames[i], b = this.frames[i + 1]
      if (jd >= a.jd && jd <= b.jd) {
        const t = (jd - a.jd) / (b.jd - a.jd)
        return new LunarFrameLerp(
          EclipseUtil.lerp(a.moonX, b.moonX, t),
          EclipseUtil.lerp(a.moonY, b.moonY, t),
          EclipseUtil.lerp(a.moonRadius, b.moonRadius, t),
          EclipseUtil.lerp(a.umbraRadius, b.umbraRadius, t),
          EclipseUtil.lerp(a.penumbraRadius, b.penumbraRadius, t),
          jd,
          EclipseUtil.lerp(a.coverage, b.coverage, t)
        )
      }
    }
    return wrap(this.frames[this.frames.length - 1])
  }

  private mix(a: string, b: string, t: number): string {
    const pa = this.parseHex(a), pb = this.parseHex(b)
    const r = Math.round(pa[0] + (pb[0] - pa[0]) * t)
    const g = Math.round(pa[1] + (pb[1] - pa[1]) * t)
    const bl = Math.round(pa[2] + (pb[2] - pa[2]) * t)
    return '#' + this.h2(r) + this.h2(g) + this.h2(bl)
  }
  private parseHex(s: string): number[] {
    if (s.charAt(0) === '#') s = s.substring(1)
    const v = parseInt(s, 16)
    return [(v >> 16) & 0xFF, (v >> 8) & 0xFF, v & 0xFF]
  }
  private h2(v: number): string {
    const s = v.toString(16)
    return s.length >= 2 ? s : '0' + s
  }
}
