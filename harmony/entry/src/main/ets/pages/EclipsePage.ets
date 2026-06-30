import { AppColors, AppDimens } from '../common/Constants'
import { EclipseUtil } from '../common/EclipseUtil'
import { EclipseShareUtil } from '../common/EclipseShareUtil'
import { Cities, ObserverCity } from '../common/Cities'
import {
  NativeBridge, SolarEclipseItem, LunarEclipseItem,
  SolarEclipsePath, LocalSolarEclipse, LunarEclipseDetail, EclipseGeoPoint
} from '../services/NativeBridge'
import {
  EclipseTimelineBar, PlaybackController, LocationStrip,
  SolarPathMapCanvas, SolarLocalDiscCanvas, LunarDiscCanvas,
  GradientBadge, TimelineMark
} from '../components/EclipseComponents'
import common from '@ohos.app.ability.common'

// ════════════════════════════════════════════════════════════════
//  EclipsePage — 与 Android/iOS 端 1:1 对齐
//
//  · 卡片含 emoji + 渐变徽章 + Saros/食季
//  · 日食详情 3 tab: 食带地图 / 本地观测 / 数据
//  · 月食详情 2 tab: 过程动画 / 数据
//  · ICS 导出 (HarmonyOS 通过 cacheDir + 弹窗提示)
// ════════════════════════════════════════════════════════════════

@Component
export struct EclipsePage {
  @State isSolar: boolean = true
  @State yearStr: string = `${new Date().getFullYear()}`
  @State monthStr: string = '1'
  @State countStr: string = '8'

  @State solarItems: SolarEclipseItem[] = []
  @State lunarItems: LunarEclipseItem[] = []
  @State loading: boolean = false
  @State searched: boolean = false

  // ArkTS 规则下同一组件只能挂一个 bindSheet, 用 sheetMode 复用同一容器.
  //   sheetMode='solar' 显示日食详情; 'lunar' 显示月食详情.
  @State showSheet: boolean = false
  @State sheetMode: string = 'solar'
  @State solarSelected: SolarEclipseItem | null = null
  @State lunarSelected: LunarEclipseItem | null = null

  search(): void {
    this.loading = true
    this.searched = true
    const y = parseInt(this.yearStr) || new Date().getFullYear()
    const m = Math.max(1, Math.min(12, parseInt(this.monthStr) || 1))
    const c = Math.max(1, Math.min(30, parseInt(this.countStr) || 8))
    setTimeout(() => {
      if (this.isSolar) {
        this.solarItems = NativeBridge.searchSolarEclipses(y, m, c)
        this.lunarItems = []
      } else {
        this.lunarItems = NativeBridge.searchLunarEclipses(y, m, c)
        this.solarItems = []
      }
      this.loading = false
    }, 16)
  }

  build() {
    Scroll() {
      Column() {
        this.headerBanner()
        this.searchPanel()
        if (this.searched) {
          if (this.isSolar) this.solarList()
          else this.lunarList()
        }
        Blank().height(24)
      }.width('100%')
    }
    .width('100%').height('100%').backgroundColor(AppColors.BACKGROUND)
    .scrollBar(BarState.Auto)
    .bindSheet($$this.showSheet, this.detailSheetBuilder(), {
      height: SheetSize.LARGE, dragBar: true,
      title: { title: this.sheetMode === 'lunar' ? '月食详情' : '日食详情' },
      onDisappear: (): void => {
        this.solarSelected = null
        this.lunarSelected = null
      }
    })
  }

  @Builder
  detailSheetBuilder() {
    if (this.sheetMode === 'lunar') {
      this.lunarSheet()
    } else {
      this.solarSheet()
    }
  }

  @Builder headerBanner() {
    Column({ space: 4 }) {
      Text('🌑 日月食工具').fontSize(AppDimens.FONT_TITLE)
        .fontWeight(FontWeight.Bold).fontColor(AppColors.ON_PRIMARY)
      Text('天文算法 · 动画 · 食带地图 · Saros')
        .fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.ON_PRIMARY).opacity(0.85)
    }
    .width('100%').padding({ top: AppDimens.PADDING_XL, bottom: AppDimens.PADDING_LG })
    .linearGradient({ angle: 135,
      colors: [[AppColors.GRADIENT_START, 0], [AppColors.GRADIENT_END, 1]] })
  }

  @Builder searchPanel() {
    Column({ space: 10 }) {
      Row() {
        this.segItem('日食', true); this.segItem('月食', false)
      }.width('100%').height(36).backgroundColor(AppColors.BACKGROUND).borderRadius(20)
      Row({ space: 6 }) {
        Text('起始').fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.TEXT_SECONDARY)
        TextInput({ text: this.yearStr }).type(InputType.Number)
          .height(38).layoutWeight(1.4).backgroundColor(AppColors.BACKGROUND)
          .borderRadius(AppDimens.BORDER_RADIUS_SM)
          .border({ width: 1, color: AppColors.DIVIDER }).textAlign(TextAlign.Center)
          .onChange((v: string) => { this.yearStr = v })
        Text('年').fontSize(AppDimens.FONT_BODY)
        TextInput({ text: this.monthStr }).type(InputType.Number)
          .height(38).width(50).backgroundColor(AppColors.BACKGROUND)
          .borderRadius(AppDimens.BORDER_RADIUS_SM)
          .border({ width: 1, color: AppColors.DIVIDER }).textAlign(TextAlign.Center)
          .onChange((v: string) => { this.monthStr = v })
        Text('月起').fontSize(AppDimens.FONT_BODY)
        TextInput({ text: this.countStr }).type(InputType.Number)
          .height(38).width(50).backgroundColor(AppColors.BACKGROUND)
          .borderRadius(AppDimens.BORDER_RADIUS_SM)
          .border({ width: 1, color: AppColors.DIVIDER }).textAlign(TextAlign.Center)
          .onChange((v: string) => { this.countStr = v })
        Text('条').fontSize(AppDimens.FONT_BODY)
      }.width('100%')
      Button(this.loading ? '搜索中…' : '开始搜索')
        .width('100%').height(42).fontSize(AppDimens.FONT_BODY)
        .fontColor(AppColors.ON_PRIMARY).backgroundColor(AppColors.PRIMARY)
        .borderRadius(21).enabled(!this.loading)
        .onClick(() => this.search())
    }
    .width('90%').padding(AppDimens.PADDING_MD)
    .backgroundColor(AppColors.SURFACE).borderRadius(AppDimens.BORDER_RADIUS_LG)
    .margin({ top: AppDimens.PADDING_MD, bottom: AppDimens.PADDING_MD })
    .shadow({ radius: 6, color: AppColors.CARD_SHADOW_COLOR, offsetY: 2 })
  }

  @Builder segItem(label: string, v: boolean) {
    Text(label)
      .layoutWeight(1).height(36).textAlign(TextAlign.Center)
      .fontSize(AppDimens.FONT_BODY)
      .fontColor(this.isSolar === v ? AppColors.ON_PRIMARY : AppColors.TEXT_SECONDARY)
      .fontWeight(this.isSolar === v ? FontWeight.Bold : FontWeight.Normal)
      .backgroundColor(this.isSolar === v ? AppColors.PRIMARY : '#00000000')
      .borderRadius(20)
      .onClick(() => { this.isSolar = v })
  }

  @Builder solarList() {
    if (this.solarItems.length === 0) {
      this.emptyHint('起始时间范围内没有日食')
    } else {
      Column({ space: 8 }) {
        ForEach(this.solarItems, (item: SolarEclipseItem) => {
          this.solarCard(item)
        }, (item: SolarEclipseItem) => `${item.jd}`)
      }.width('90%')
    }
  }

  @Builder lunarList() {
    if (this.lunarItems.length === 0) {
      this.emptyHint('起始时间范围内没有月食')
    } else {
      Column({ space: 8 }) {
        ForEach(this.lunarItems, (item: LunarEclipseItem) => {
          this.lunarCard(item)
        }, (item: LunarEclipseItem) => `${item.jd}`)
      }.width('90%')
    }
  }

  @Builder emptyHint(text: string) {
    Text(text).fontSize(AppDimens.FONT_BODY).fontColor(AppColors.TEXT_SECONDARY)
      .width('90%').textAlign(TextAlign.Center).padding(AppDimens.PADDING_XL)
  }

  @Builder solarCard(item: SolarEclipseItem) {
    Column({ space: 8 }) {
      Row({ space: 8 }) {
        Text(EclipseUtil.solarEmoji(item.type)).fontSize(22)
        Column({ space: 2 }) {
          Text(this.fmtDateTime(item) + ' TD')
            .fontSize(AppDimens.FONT_BODY).fontWeight(FontWeight.Bold)
            .fontColor(AppColors.ON_SURFACE)
          Text(item.season + '  ·  Saros ' + item.saros)
            .fontSize(11).fontColor(AppColors.TEXT_SECONDARY)
        }.alignItems(HorizontalAlign.Start).layoutWeight(1)
        GradientBadge({
          text: item.typeName,
          color: EclipseUtil.solarBadgeColor(item.type),
          icon: EclipseUtil.solarEmoji(item.type)
        })
      }.width('100%')
      Row({ space: 8 }) {
        Text('食分').fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.TEXT_SECONDARY)
        Progress({ value: EclipseUtil.magnitudeProgress(item.magnitude) * 100, total: 100, type: ProgressType.Linear })
          .layoutWeight(1).height(6)
          .color(EclipseUtil.solarBadgeColor(item.type))
          .backgroundColor(AppColors.DIVIDER)
        Text(item.magnitude.toFixed(3))
          .fontSize(AppDimens.FONT_CAPTION).fontWeight(FontWeight.Bold).fontColor(AppColors.ON_SURFACE)
      }.width('100%')
      Text('点击 → 食带地图 · 本地观测 · 数据')
        .fontSize(11).fontColor(AppColors.PRIMARY).fontWeight(FontWeight.Medium)
    }
    .width('100%').padding(AppDimens.PADDING_MD)
    .backgroundColor(AppColors.SURFACE).borderRadius(AppDimens.BORDER_RADIUS_MD)
    .shadow({ radius: 4, color: AppColors.CARD_SHADOW_COLOR, offsetY: 1 })
    .onClick(() => {
      this.solarSelected = item
      this.lunarSelected = null
      this.sheetMode = 'solar'
      this.showSheet = true
    })
  }

  @Builder lunarCard(item: LunarEclipseItem) {
    Column({ space: 8 }) {
      Row({ space: 8 }) {
        Text(EclipseUtil.lunarEmoji(item.type)).fontSize(22)
        Column({ space: 2 }) {
          Text(this.fmtDateTime(item) + ' TD')
            .fontSize(AppDimens.FONT_BODY).fontWeight(FontWeight.Bold)
            .fontColor(AppColors.ON_SURFACE)
          Text(item.season + '  ·  Saros ' + item.saros)
            .fontSize(11).fontColor(AppColors.TEXT_SECONDARY)
        }.alignItems(HorizontalAlign.Start).layoutWeight(1)
        GradientBadge({
          text: item.typeName,
          color: EclipseUtil.lunarBadgeColor(item.type),
          icon: EclipseUtil.lunarEmoji(item.type)
        })
      }.width('100%')
      if (item.type !== 'B') {
        Row({ space: 8 }) {
          Text('食分').fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.TEXT_SECONDARY)
          Progress({ value: EclipseUtil.magnitudeProgress(item.magnitude) * 100, total: 100, type: ProgressType.Linear })
            .layoutWeight(1).height(6)
            .color(EclipseUtil.lunarBadgeColor(item.type))
            .backgroundColor(AppColors.DIVIDER)
          Text(item.magnitude.toFixed(3))
            .fontSize(AppDimens.FONT_CAPTION).fontWeight(FontWeight.Bold).fontColor(AppColors.ON_SURFACE)
        }.width('100%')
      } else {
        Text('半影食 · 月亮仅穿过地球半影')
          .fontSize(11).fontColor(AppColors.TEXT_SECONDARY)
      }
      Text('点击 → 过程动画 · 数据')
        .fontSize(11).fontColor(AppColors.PRIMARY).fontWeight(FontWeight.Medium)
    }
    .width('100%').padding(AppDimens.PADDING_MD)
    .backgroundColor(AppColors.SURFACE).borderRadius(AppDimens.BORDER_RADIUS_MD)
    .shadow({ radius: 4, color: AppColors.CARD_SHADOW_COLOR, offsetY: 1 })
    .onClick(() => {
      this.lunarSelected = item
      this.solarSelected = null
      this.sheetMode = 'lunar'
      this.showSheet = true
    })
  }

  // ─── Solar Sheet ────────────────────────────────────────
  @State solarTab: number = 0
  @State solarPath: SolarEclipsePath | null = null
  @State solarLocal: LocalSolarEclipse | null = null
  @State solarDitu0: number[] = []
  @State solarDitu1: number[] = []
  @State solarDitu2: number[] = []
  @State solarUseBig: boolean = true
  @State solarCity: ObserverCity = Cities.DEFAULT
  @State solarJd: number = 0
  // dataVersion: 数据加载完成后递增, 驱动 Canvas 通过 @Watch 重绘
  @State solarDataVersion: number = 0

  @Builder solarSheet() {
    if (this.solarSelected) {
      Column() {
        Row({ space: 8 }) {
          Text(EclipseUtil.solarEmoji(this.solarSelected!.type)).fontSize(26)
          Column({ space: 2 }) {
            Text(this.solarSelected!.typeName + '日食').fontSize(16)
              .fontWeight(FontWeight.Bold).fontColor(AppColors.PRIMARY)
            Text(EclipseUtil.jdTdToLocal(this.solarSelected!.jd, 8, true) + ' 北京时')
              .fontSize(11).fontColor(AppColors.TEXT_SECONDARY).fontFamily('monospace')
            Text('Saros ' + this.solarSelected!.saros + '  ·  #' + this.solarSelected!.sarosMember
              + '  ·  ' + this.solarSelected!.season)
              .fontSize(10).fontColor(AppColors.LUNAR_TEXT)
          }.alignItems(HorizontalAlign.Start).layoutWeight(1)
          GradientBadge({
            text: this.solarSelected!.typeName,
            color: EclipseUtil.solarBadgeColor(this.solarSelected!.type),
            icon: EclipseUtil.solarEmoji(this.solarSelected!.type)
          })
        }.width('100%').padding({ left: AppDimens.PADDING_MD, right: AppDimens.PADDING_MD })

        Row({ space: 8 }) {
          Text('📅 导出 ICS')
            .fontSize(11).fontColor(AppColors.PRIMARY).fontWeight(FontWeight.Medium)
            .backgroundColor(AppColors.BACKGROUND).borderRadius(8)
            .padding({ left: 10, right: 10, top: 6, bottom: 6 })
            .onClick(() => this.exportSolarICS())
        }.width('100%').padding({ left: AppDimens.PADDING_MD, right: AppDimens.PADDING_MD, top: 6 })

        this.tabBar(['🗺  食带地图', '🌞  本地观测', '📊  数据'], this.solarTab,
          (i: number) => { this.solarTab = i })

        Scroll() {
          Column({ space: 10 }) {
            if (this.solarTab === 0) this.solarTabMap()
            else if (this.solarTab === 1) this.solarTabLocal()
            else this.solarTabData()
            Blank().height(20)
          }.width('100%').padding(AppDimens.PADDING_MD)
        }.layoutWeight(1)
      }.width('100%').height('100%')
      .onAppear(() => this.loadSolarDetail())
    }
  }

  @Builder solarTabMap() {
    Row({ space: 8 }) {
      Text('地图').fontSize(11).fontColor(AppColors.TEXT_SECONDARY)
      Toggle({ type: ToggleType.Switch, isOn: this.solarUseBig })
        .onChange((v: boolean) => { this.solarUseBig = v })
      Text(this.solarUseBig ? '大图 (ditu1+2)' : '小图 (ditu0)')
        .fontSize(11).fontColor(AppColors.TEXT_SECONDARY)
    }
    if (this.solarPath) {
      SolarPathMapCanvas({
        worldMapRad: this.solarUseBig && this.solarDitu1.length > 0 ? this.solarDitu1 : this.solarDitu0,
        bordersRad: this.solarUseBig && this.solarDitu2.length > 0 ? this.solarDitu2 : [],
        path: this.solarPath,
        currentJd: this.solarJd,
        dataVersion: this.solarDataVersion
      })
      EclipseTimelineBar({
        jdBegin: this.solarPath.partialStart.jd,
        jdEnd: this.solarPath.partialEnd.jd,
        currentJd: $solarJd,
        marks: this.solarMarks()
      })
      PlaybackController({
        jdBegin: this.solarPath.partialStart.jd,
        jdEnd: this.solarPath.partialEnd.jd,
        currentJd: $solarJd
      })
    } else {
      Text('正在加载食带数据…')
        .fontSize(AppDimens.FONT_BODY).fontColor(AppColors.TEXT_SECONDARY)
    }
  }

  @Builder solarTabLocal() {
    LocationStrip({ selected: $solarCity })
    if (this.solarLocal && this.solarLocal.frames.length > 0) {
      SolarLocalDiscCanvas({
        frames: this.solarLocal.frames,
        currentJd: this.solarJd,
        dataVersion: this.solarDataVersion,
        tzHours: this.solarCity.tzHours
      })
      EclipseTimelineBar({
        jdBegin: this.solarLocal.frames[0].jd,
        jdEnd: this.solarLocal.frames[this.solarLocal.frames.length - 1].jd,
        currentJd: $solarJd,
        marks: this.localSolarMarks()
      })
      PlaybackController({
        jdBegin: this.solarLocal.frames[0].jd,
        jdEnd: this.solarLocal.frames[this.solarLocal.frames.length - 1].jd,
        currentJd: $solarJd
      })
      this.paramGrid([
        ['当地最大食分', this.solarLocal.maxMagnitude.toFixed(3)],
        ['初亏 C₁',     EclipseUtil.jdTdToLocal(this.solarLocal.tC1, this.solarCity.tzHours)],
        ['食既 C₂',     EclipseUtil.jdTdToLocal(this.solarLocal.tC2, this.solarCity.tzHours)],
        ['食甚',        EclipseUtil.jdTdToLocal(this.solarLocal.tMax, this.solarCity.tzHours)],
        ['生光 C₃',     EclipseUtil.jdTdToLocal(this.solarLocal.tC3, this.solarCity.tzHours)],
        ['复圆 C₄',     EclipseUtil.jdTdToLocal(this.solarLocal.tC4, this.solarCity.tzHours)]
      ])
    } else {
      Text('该地点不可见').fontSize(AppDimens.FONT_BODY).fontColor(AppColors.TEXT_SECONDARY)
    }
  }

  @Builder solarTabData() {
    this.paramGrid([
      ['日期', this.fmtDate(this.solarSelected!.year, this.solarSelected!.month, this.solarSelected!.day)],
      ['时刻 (TD)', this.pad2(this.solarSelected!.hour) + ':' + this.pad2(this.solarSelected!.minute)],
      ['北京时', EclipseUtil.jdTdToLocal(this.solarSelected!.jd, 8)],
      ['食分', this.solarSelected!.magnitude.toFixed(4)],
      ['γ', (this.solarSelected!.gamma >= 0 ? '+' : '') + this.solarSelected!.gamma.toFixed(4)],
      ['最长时长', EclipseUtil.formatDuration(this.solarSelected!.duration)],
      ['食带宽', this.solarSelected!.pathWidth > 0 ? this.solarSelected!.pathWidth.toFixed(1) + ' km' : '—'],
      ['类型代码', this.solarSelected!.type],
      ['中心点', (this.solarSelected!.centerLon < 99 && this.solarSelected!.centerLat < 99)
        ? EclipseUtil.lonLabel(this.solarSelected!.centerLon) + ', ' + EclipseUtil.latLabel(this.solarSelected!.centerLat)
        : '无中心'],
      ['Saros', this.solarSelected!.saros + ' (#' + this.solarSelected!.sarosMember + ')'],
      ['食季', this.solarSelected!.season]
    ])
    if (this.solarPath) {
      Divider().color(AppColors.DIVIDER)
      Text('全球路径关键点 (TD)').fontSize(AppDimens.FONT_SUBTITLE)
        .fontWeight(FontWeight.Bold).fontColor(AppColors.PRIMARY)
      this.eventPoint('偏食始 P₁', this.solarPath.partialStart)
      this.eventPoint('中心始 C₁', this.solarPath.centralStart)
      this.eventPoint('食甚 MAX',  this.solarPath.greatestEclipse)
      this.eventPoint('中心终 C₂', this.solarPath.centralEnd)
      this.eventPoint('偏食终 P₄', this.solarPath.partialEnd)
    }
  }

  // 延迟到 sheet 弹出动画完成后再做重计算, 避免在动画帧里阻塞主线程.
  // 每步独立 try/catch: 即使某次 NAPI 抛异常也不会污染整个 sheet.
  loadSolarDetail(): void {
    if (!this.solarSelected) return
    const it = this.solarSelected
    setTimeout(() => {
      try {
        if (this.solarDitu0.length === 0) this.solarDitu0 = NativeBridge.getWorldMapDitu0()
        if (this.solarDitu1.length === 0) this.solarDitu1 = NativeBridge.getWorldMapData(1)
        if (this.solarDitu2.length === 0) this.solarDitu2 = NativeBridge.getWorldMapData(2)
      } catch (e) {
        this.safeToast('地图数据加载失败: ' + (e as Error).message)
      }
      try {
        this.solarPath = NativeBridge.getSolarEclipsePath(it.year, it.month, it.day)
        this.solarJd = this.solarPath ? this.solarPath.maxEclipseJd : it.jd
      } catch (e) {
        this.solarPath = null
        this.safeToast('食带计算失败: ' + (e as Error).message)
      }
      try {
        this.solarLocal = NativeBridge.getLocalSolarEclipse(
          it.year, it.month, it.day,
          this.solarCity.longitude, this.solarCity.latitude, 60)
      } catch (e) {
        this.solarLocal = null
      }
      // 数据全部就位, 通知 Canvas 重绘.
      this.solarDataVersion = this.solarDataVersion + 1
    }, 16)
  }

  solarMarks(): TimelineMark[] {
    if (!this.solarPath) return []
    const p = this.solarPath
    const r: TimelineMark[] = []
    if (p.partialStart && p.partialStart.jd > 0) r.push(new TimelineMark('P1', p.partialStart.jd, AppColors.PATH_PENUMBRA))
    if (p.centralStart && p.centralStart.jd > 0 && p.centralStart.longitude < 99) r.push(new TimelineMark('C1', p.centralStart.jd, AppColors.PATH_UMBRA))
    if (p.greatestEclipse && p.greatestEclipse.jd > 0) r.push(new TimelineMark('MAX', p.greatestEclipse.jd, AppColors.PATH_CENTER))
    if (p.centralEnd && p.centralEnd.jd > 0 && p.centralEnd.longitude < 99) r.push(new TimelineMark('C2', p.centralEnd.jd, AppColors.PATH_UMBRA))
    if (p.partialEnd && p.partialEnd.jd > 0) r.push(new TimelineMark('P4', p.partialEnd.jd, AppColors.PATH_PENUMBRA))
    return r
  }
  localSolarMarks(): TimelineMark[] {
    if (!this.solarLocal) return []
    const l = this.solarLocal
    const r: TimelineMark[] = []
    if (l.tC1 > 0) r.push(new TimelineMark('初亏', l.tC1, AppColors.PATH_PENUMBRA))
    if (l.tC2 > 0) r.push(new TimelineMark('食既', l.tC2, AppColors.PATH_UMBRA))
    if (l.tMax > 0) r.push(new TimelineMark('食甚', l.tMax, AppColors.PATH_CENTER))
    if (l.tC3 > 0) r.push(new TimelineMark('生光', l.tC3, AppColors.PATH_UMBRA))
    if (l.tC4 > 0) r.push(new TimelineMark('复圆', l.tC4, AppColors.PATH_PENUMBRA))
    return r
  }

  exportSolarICS(): void {
    if (!this.solarSelected) return
    const it = this.solarSelected
    try {
      const ctx = this.getUIContext().getHostContext() as common.UIAbilityContext
      const cacheDir = ctx.cacheDir
      const begin = this.solarPath ? this.solarPath.partialStart.jd : it.jd - 1.5 / 24
      const end = this.solarPath ? this.solarPath.partialEnd.jd : it.jd + 1.5 / 24
      const ics = EclipseShareUtil.buildICS(
        'solar-' + it.year + it.month + it.day + '@sxwnl',
        it.typeName + ' (Saros ' + it.saros + ')',
        '类型: ' + it.typeName + ' (' + it.type + ')\n食分: ' + it.magnitude.toFixed(3)
          + '\nSaros: ' + it.saros + ' (#' + it.sarosMember + ')\n' + it.season,
        begin, end)
      const path = EclipseShareUtil.writeICS(cacheDir, ics,
        'solar_' + it.year + '_' + it.month + '_' + it.day + '.ics')
      this.safeToast('ICS 已保存: ' + path)
    } catch (e) {
      this.safeToast('ICS 导出失败')
    }
  }

  /** 包装 PromptAction.showToast, 永不外抛. */
  private safeToast(msg: string): void {
    try {
      this.getUIContext().getPromptAction().showToast({ message: msg, duration: 2500 })
    } catch (_) {
      // 静默吞掉, Toast 本身失败无需再 Toast
    }
  }

  // ─── Lunar Sheet ──────────────────────────────────────
  @State lunarTab: number = 0
  @State lunarDetail: LunarEclipseDetail | null = null
  @State lunarJd: number = 0
  @State lunarDataVersion: number = 0

  @Builder lunarSheet() {
    if (this.lunarSelected) {
      Column() {
        Row({ space: 8 }) {
          Text(EclipseUtil.lunarEmoji(this.lunarSelected!.type)).fontSize(26)
          Column({ space: 2 }) {
            Text(this.lunarSelected!.typeName + '月食').fontSize(16)
              .fontWeight(FontWeight.Bold).fontColor(AppColors.PRIMARY)
            Text(EclipseUtil.jdTdToLocal(this.lunarSelected!.jd, 8, true) + ' 北京时')
              .fontSize(11).fontColor(AppColors.TEXT_SECONDARY).fontFamily('monospace')
            Text('Saros ' + this.lunarSelected!.saros + '  ·  #' + this.lunarSelected!.sarosMember
              + '  ·  ' + this.lunarSelected!.season)
              .fontSize(10).fontColor(AppColors.LUNAR_TEXT)
          }.alignItems(HorizontalAlign.Start).layoutWeight(1)
          GradientBadge({
            text: this.lunarSelected!.typeName,
            color: EclipseUtil.lunarBadgeColor(this.lunarSelected!.type),
            icon: EclipseUtil.lunarEmoji(this.lunarSelected!.type)
          })
        }.width('100%').padding({ left: AppDimens.PADDING_MD, right: AppDimens.PADDING_MD })

        Row({ space: 8 }) {
          Text('📅 导出 ICS')
            .fontSize(11).fontColor(AppColors.PRIMARY).fontWeight(FontWeight.Medium)
            .backgroundColor(AppColors.BACKGROUND).borderRadius(8)
            .padding({ left: 10, right: 10, top: 6, bottom: 6 })
            .onClick(() => this.exportLunarICS())
        }.width('100%').padding({ left: AppDimens.PADDING_MD, right: AppDimens.PADDING_MD, top: 6 })

        this.tabBar(['🌙  过程动画', '📊  数据'], this.lunarTab,
          (i: number) => { this.lunarTab = i })

        Scroll() {
          Column({ space: 10 }) {
            if (this.lunarTab === 0) this.lunarTabAnim()
            else this.lunarTabData()
            Blank().height(20)
          }.width('100%').padding(AppDimens.PADDING_MD)
        }.layoutWeight(1)
      }.width('100%').height('100%')
      .onAppear(() => this.loadLunarDetail())
    }
  }

  @Builder lunarTabAnim() {
    if (this.lunarDetail && this.lunarDetail.frames.length > 0) {
      LunarDiscCanvas({
        frames: this.lunarDetail.frames,
        currentJd: this.lunarJd,
        dataVersion: this.lunarDataVersion
      })
      EclipseTimelineBar({
        jdBegin: this.lunarDetail.frames[0].jd,
        jdEnd: this.lunarDetail.frames[this.lunarDetail.frames.length - 1].jd,
        currentJd: $lunarJd,
        marks: this.lunarMarks()
      })
      PlaybackController({
        jdBegin: this.lunarDetail.frames[0].jd,
        jdEnd: this.lunarDetail.frames[this.lunarDetail.frames.length - 1].jd,
        currentJd: $lunarJd
      })
    }
  }

  @Builder lunarTabData() {
    this.paramGrid([
      ['日期', this.fmtDate(this.lunarSelected!.year, this.lunarSelected!.month, this.lunarSelected!.day)],
      ['时刻 (TD)', this.pad2(this.lunarSelected!.hour) + ':' + this.pad2(this.lunarSelected!.minute)],
      ['北京时', EclipseUtil.jdTdToLocal(this.lunarSelected!.jd, 8)],
      ['类型', this.lunarSelected!.typeName],
      ['食分', this.lunarSelected!.type === 'B' ? '—' : this.lunarSelected!.magnitude.toFixed(4)],
      ['Saros', this.lunarSelected!.saros + ' (#' + this.lunarSelected!.sarosMember + ')'],
      ['食季', this.lunarSelected!.season]
    ])
    if (this.lunarDetail) {
      Divider().color(AppColors.DIVIDER)
      Text('关键时刻 (TD)').fontSize(AppDimens.FONT_SUBTITLE)
        .fontWeight(FontWeight.Bold).fontColor(AppColors.PRIMARY)
      this.paramGrid([
        ['半影始 P₁', EclipseUtil.jdToTime(this.lunarDetail.tP1)],
        ['初亏 U₁', EclipseUtil.jdToTime(this.lunarDetail.tU1)],
        ['食既 U₂', EclipseUtil.jdToTime(this.lunarDetail.tU2)],
        ['食甚',    EclipseUtil.jdToTime(this.lunarDetail.tMax)],
        ['生光 U₃', EclipseUtil.jdToTime(this.lunarDetail.tU3)],
        ['复圆 U₄', EclipseUtil.jdToTime(this.lunarDetail.tU4)],
        ['半影终 P₄', EclipseUtil.jdToTime(this.lunarDetail.tP4)]
      ])
    }
  }

  loadLunarDetail(): void {
    if (!this.lunarSelected) return
    const it = this.lunarSelected
    setTimeout(() => {
      try {
        this.lunarDetail = NativeBridge.getLunarEclipseDetail(it.year, it.month, it.day, 120)
      } catch (e) {
        this.lunarDetail = null
        this.safeToast('月食过程计算失败: ' + (e as Error).message)
      }
      this.lunarJd = this.lunarDetail ? this.lunarDetail.tMax : it.jd
      this.lunarDataVersion = this.lunarDataVersion + 1
    }, 16)
  }

  lunarMarks(): TimelineMark[] {
    if (!this.lunarDetail) return []
    const d = this.lunarDetail
    const r: TimelineMark[] = []
    if (d.tP1 > 0) r.push(new TimelineMark('P1', d.tP1, AppColors.PATH_PENUMBRA))
    if (d.tU1 > 0) r.push(new TimelineMark('U1', d.tU1, AppColors.PATH_UMBRA))
    if (d.tU2 > 0) r.push(new TimelineMark('U2', d.tU2, AppColors.MOON_BLOOD))
    if (d.tMax > 0) r.push(new TimelineMark('MAX', d.tMax, AppColors.PATH_CENTER))
    if (d.tU3 > 0) r.push(new TimelineMark('U3', d.tU3, AppColors.MOON_BLOOD))
    if (d.tU4 > 0) r.push(new TimelineMark('U4', d.tU4, AppColors.PATH_UMBRA))
    if (d.tP4 > 0) r.push(new TimelineMark('P4', d.tP4, AppColors.PATH_PENUMBRA))
    return r
  }

  exportLunarICS(): void {
    if (!this.lunarSelected) return
    const it = this.lunarSelected
    try {
      const ctx = this.getUIContext().getHostContext() as common.UIAbilityContext
      const cacheDir = ctx.cacheDir
      const begin = this.lunarDetail && this.lunarDetail.tP1 > 0 ? this.lunarDetail.tP1 : it.jd - 1.5 / 24
      const end = this.lunarDetail && this.lunarDetail.tP4 > 0 ? this.lunarDetail.tP4 : it.jd + 1.5 / 24
      const ics = EclipseShareUtil.buildICS(
        'lunar-' + it.year + it.month + it.day + '@sxwnl',
        it.typeName + ' (Saros ' + it.saros + ')',
        '类型: ' + it.typeName + ' (' + it.type + ')\n食分: ' + it.magnitude.toFixed(3)
          + '\nSaros: ' + it.saros + ' (#' + it.sarosMember + ')\n' + it.season,
        begin, end)
      const path = EclipseShareUtil.writeICS(cacheDir, ics,
        'lunar_' + it.year + '_' + it.month + '_' + it.day + '.ics')
      this.safeToast('ICS 已保存: ' + path)
    } catch (e) {
      this.safeToast('ICS 导出失败')
    }
  }

  // ─── 共用 ─────────────────────────────────────────────

  @Builder tabBar(tabs: string[], selected: number, onSelect: (i: number) => void) {
    Row() {
      ForEach(tabs, (t: string, idx: number) => {
        Column({ space: 4 }) {
          Text(t).fontSize(AppDimens.FONT_BODY)
            .fontWeight(selected === idx ? FontWeight.Bold : FontWeight.Normal)
            .fontColor(selected === idx ? AppColors.PRIMARY : AppColors.TEXT_SECONDARY)
          Row().width('100%').height(2)
            .backgroundColor(selected === idx ? AppColors.PRIMARY : '#00000000')
        }
        .layoutWeight(1).padding({ top: 8, bottom: 8 })
        .onClick(() => onSelect(idx))
      })
    }
    .width('100%').backgroundColor(AppColors.SURFACE)
  }

  @Builder paramGrid(pairs: string[][]) {
    Column({ space: 4 }) {
      ForEach(this.chunk(pairs, 2), (row: string[][]) => {
        Row({ space: 10 }) {
          ForEach(row, (kv: string[]) => {
            Row() {
              Text(kv[0]).fontSize(11).fontColor(AppColors.TEXT_SECONDARY).width(80)
              Text(kv[1]).fontSize(12).fontColor(AppColors.ON_SURFACE)
                .fontWeight(FontWeight.Medium).fontFamily('monospace')
            }.layoutWeight(1)
            .padding({ left: 8, right: 8, top: 5, bottom: 5 })
            .backgroundColor(AppColors.BACKGROUND).borderRadius(6)
          })
          if (row.length === 1) Blank().layoutWeight(1)
        }.width('100%')
      })
    }.width('100%')
  }

  @Builder eventPoint(label: string, p: EclipseGeoPoint) {
    Row() {
      Text(label).fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.TEXT_SECONDARY).width(110)
      Column() {
        Text(p.longitude < 99
          ? p.longitude.toFixed(2) + '°,  ' + p.latitude.toFixed(2) + '°' : '—')
          .fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.ON_SURFACE)
          .fontFamily('monospace')
        Text(EclipseUtil.jdToDateTime(p.jd))
          .fontSize(11).fontColor(AppColors.LUNAR_TEXT).fontFamily('monospace')
      }.alignItems(HorizontalAlign.Start).layoutWeight(1)
    }
    .width('100%').padding({ left: 10, right: 10, top: 6, bottom: 6 })
    .backgroundColor(AppColors.BACKGROUND).borderRadius(6)
  }

  // ─── helpers ───────────────────────────────────────────
  fmtDateTime(item: SolarEclipseItem | LunarEclipseItem): string {
    return this.pad4(item.year) + '-' + this.pad2(item.month) + '-' + this.pad2(item.day)
      + ' ' + this.pad2(item.hour) + ':' + this.pad2(item.minute)
  }
  fmtDate(y: number, mo: number, d: number): string {
    return this.pad4(y) + '-' + this.pad2(mo) + '-' + this.pad2(d)
  }
  pad2(v: number): string {
    const s = Math.abs(v).toString()
    return (s.length >= 2 ? '' : '0') + s
  }
  pad4(v: number): string {
    const s = Math.abs(v).toString()
    const sign = v < 0 ? '-' : ''
    return sign + (s.length >= 4 ? s : '0'.repeat(4 - s.length) + s)
  }
  chunk<T>(arr: T[], n: number): T[][] {
    const r: T[][] = []
    for (let i = 0; i < arr.length; i += n) r.push(arr.slice(i, i + n))
    return r
  }
}

