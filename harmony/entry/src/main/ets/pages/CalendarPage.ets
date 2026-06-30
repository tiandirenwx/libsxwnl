import { AppColors, AppDimens, WEEK_NAMES } from '../common/Constants'
import { YearUtil } from '../common/YearUtil'
import { NativeBridge, DayInfo, DayRTS, DayAlmanac, AlmanacTopic,
         GeoProvince, GeoCity, TimezoneGroup } from '../services/NativeBridge'
import { AlmanacInlineRow, AlmanacSheet, AlmanacTopicsSheet } from '../components/AlmanacComponents'
import { Location, getDefaultLocation, formatLonLat, formatTimezone } from '../common/Locations'

// 月相事件 / 节气事件 (用于月历底部摘要)
interface MonthEvent {
  day: number      // 公历日
  time: string     // 时刻 HH:MM:SS
  name: string     // 名称: 朔月/望月/上弦/下弦 或 节气名
}

// ────────────────────────────────────────────────────────────
//  月历页面: 月份网格 + hover/click 浮层详情
//  支持: 年/月输入跳转, 上/下年, 上/下月, 今天
//  详情字段对应 sxwnl/lunar.js showMessD() 中的展示
// ────────────────────────────────────────────────────────────
@Component
export struct CalendarPage {
  @State currentYear: number = new Date().getFullYear()
  @State currentMonth: number = new Date().getMonth() + 1
  @State yearInput: string = YearUtil.astroYearToStr(new Date().getFullYear())
  @State monthInput: string = `${new Date().getMonth() + 1}`
  @State monthDays: DayInfo[] = []
  // popupKey: 当前展开浮层的单元格 key("Y-M-D"), 同时只有一个浮层展开
  @State popupKey: string = ''
  @State yearGZ: string = ''
  @State shengXiao: string = ''
  @State huangdiYear: number = 0
  @State todayYear: number = new Date().getFullYear()
  @State todayMonth: number = new Date().getMonth() + 1
  @State todayDay: number = new Date().getDate()

  // ── 底栏数据 ──
  // 当前用于计算 RTS 的目标日 (默认今天, 用户点击其它日则跟随)
  @State rtsDay: number = new Date().getDate()
  @State rts: DayRTS | null = null
  @State moonEvents: MonthEvent[] = []   // 当月朔/望/上弦/下弦
  @State jieQiEvents: MonthEvent[] = []  // 当月节气

  // ── 顶栏地点选择 (对应 JS 原版 Sel1/Sel2 国内 + Sel_zhou/Sel_dq 国际) ──
  //
  //   国内省市切换: 改变 location (经纬度+时区) → 日月升降跟着重算
  //   国际洲国切换: 仅用于顶栏"外地实时时间"显示, 不影响日历主时区
  //                 (与 JS 版本一致, 日月升降固定按本地系统时区)
  @State location: Location = getDefaultLocation()
  @State provinces: GeoProvince[] = []
  @State cities: GeoCity[] = []           // 当前省下城市
  @State provIdx: number = 0
  @State cityIdx: number = 0
  @State provOpts: SelectOption[] = []    // 省 SelectOption 缓存 (避免 build 内 map)
  @State cityOpts: SelectOption[] = []

  @State continents: string[] = []
  @State countries: TimezoneGroup[] = []  // 当前洲下国家
  @State contIdx: number = 0
  @State countryIdx: number = 0
  @State contOpts: SelectOption[] = []
  @State countryOpts: SelectOption[] = []
  private allTzGroups: TimezoneGroup[] = []

  // 双时钟 (每秒刷新)
  @State localClock: string = ''          // 本地系统时区当前时间
  @State intlClock: string = ''           // 选中国际时区当前时间
  private clockTimer: number = -1

  // ── 老黄历 (按需懒加载, 浮层中显示摘要, 抽屉中显示全量) ──
  // popupAlmanacKey 与 popupKey 对齐, 用于缓存"当前展开浮层"的老黄历摘要,
  // 避免每次重绘都跨 JNI 调用.
  @State popupAlmanacKey: string = ''
  @State popupAlmanac: DayAlmanac | null = null

  // 抽屉态: 当 sheetDay 非空时 bindSheet 显示, 关闭后清空
  //   sheetMode: 'almanac' 显示当日老黄历; 'topics' 显示静态知识(董公总论/口诀);
  //              'location' 城市选择
  //   ArkTS 规则下同一组件只能挂一个 bindSheet, 所以用 mode 复用同一 sheet 容器.
  @State sheetDay: DayInfo | null = null
  @State sheetAlmanac: DayAlmanac | null = null
  @State sheetVisible: boolean = false
  @State sheetMode: string = 'almanac'   // 'almanac' | 'topics'

  // 静态知识 (董公总论/口诀/方位) — 全局只取一次, 常驻内存
  @State topics: AlmanacTopic[] = []

  aboutToAppear() {
    this.initTopBars()
    this.loadMonthData()
    this.tickClock()
    // ArkTS 全局 setInterval 返回 number, 用于稍后清理
    this.clockTimer = setInterval((): void => { this.tickClock() }, 1000)
  }

  aboutToDisappear() {
    if (this.clockTimer >= 0) {
      clearInterval(this.clockTimer)
      this.clockTimer = -1
    }
  }

  // ── 顶栏数据初始化 ──────────────────────────────────────
  initTopBars() {
    this.provinces = NativeBridge.geoListProvinces()
    const popts: SelectOption[] = []
    for (let i = 0; i < this.provinces.length; i++) popts.push({ value: this.provinces[i].province })
    this.provOpts = popts

    // 默认匹配 location 当前值 (getDefaultLocation 一般为北京·天安门)
    let pIdx = 0
    for (let i = 0; i < this.provinces.length; i++) {
      if (this.provinces[i].province === this.location.province) { pIdx = i; break }
    }
    this.provIdx = pIdx
    this.loadCitiesOfProv(pIdx, true)
    let cIdx = 0
    for (let i = 0; i < this.cities.length; i++) {
      if (this.cities[i].district === this.location.district) { cIdx = i; break }
    }
    this.cityIdx = cIdx
    if (this.cities.length > 0) this.location = this.cities[cIdx]

    this.allTzGroups = NativeBridge.geoListTimezones()
    const conts: string[] = []
    const seen: Set<string> = new Set<string>()
    for (let i = 0; i < this.allTzGroups.length; i++) {
      const c = this.allTzGroups[i].continent
      if (!seen.has(c)) {
        seen.add(c)
        conts.push(c)
      }
    }
    this.continents = conts
    const copts: SelectOption[] = []
    for (let i = 0; i < conts.length; i++) copts.push({ value: conts[i] })
    this.contOpts = copts

    // 默认选"亚洲·中国"以匹配本地体验
    let defaultCont = 0
    for (let i = 0; i < conts.length; i++) {
      if (conts[i] === '亚洲') { defaultCont = i; break }
    }
    this.contIdx = defaultCont
    this.filterCountries(this.continents[this.contIdx])
    let defaultCountry = 0
    for (let i = 0; i < this.countries.length; i++) {
      if (this.countries[i].country === '中国') { defaultCountry = i; break }
    }
    this.countryIdx = defaultCountry
  }

  // 加载某省份下的城市. initialLoad=true 时不重置 cityIdx/location, 由调用方自定;
  // initialLoad=false 时默认选第一个城市并触发 RTS 重算.
  private loadCitiesOfProv(provIdx: number, initialLoad: boolean) {
    if (provIdx < 0 || provIdx >= this.provinces.length) return
    this.cities = NativeBridge.geoListCities(this.provinces[provIdx].province)
    const opts: SelectOption[] = []
    for (let i = 0; i < this.cities.length; i++) opts.push({ value: this.cities[i].district })
    this.cityOpts = opts
    if (initialLoad) return
    this.cityIdx = 0
    if (this.cities.length > 0) {
      this.location = this.cities[0]
      this.loadRTS()
    }
  }

  private filterCountries(continent: string) {
    const out: TimezoneGroup[] = []
    for (let i = 0; i < this.allTzGroups.length; i++) {
      if (this.allTzGroups[i].continent === continent) out.push(this.allTzGroups[i])
    }
    this.countries = out
    const opts: SelectOption[] = []
    for (let i = 0; i < out.length; i++) opts.push({ value: out[i].country })
    this.countryOpts = opts
  }

  // ── 顶栏切换回调 ────────────────────────────────────────
  onProvinceChange(idx: number) {
    if (idx < 0 || idx >= this.provinces.length) return
    this.provIdx = idx
    this.loadCitiesOfProv(idx, false)
  }

  onCityChange(idx: number) {
    if (idx < 0 || idx >= this.cities.length) return
    this.cityIdx = idx
    this.location = this.cities[idx]
    this.loadRTS()
  }

  onContinentChange(idx: number) {
    if (idx < 0 || idx >= this.continents.length) return
    this.contIdx = idx
    this.filterCountries(this.continents[idx])
    this.countryIdx = 0
    this.tickClock()
  }

  onCountryChange(idx: number) {
    if (idx < 0 || idx >= this.countries.length) return
    this.countryIdx = idx
    this.tickClock()
  }

  // 每秒刷新双时钟. 外地时间 = UTC + 选中时区
  tickClock() {
    const now = new Date()
    this.localClock = this.fmtClock(now)

    if (this.countries.length === 0 || this.countryIdx >= this.countries.length) {
      this.intlClock = ''
      return
    }
    const tz = this.countries[this.countryIdx].timezone
    // UTC 毫秒 = 本地毫秒 + 本地 offset(分钟) * 60000
    const utcMs = now.getTime() + now.getTimezoneOffset() * 60000
    const intl = new Date(utcMs + tz * 3600000)
    this.intlClock = this.fmtClock(intl)
  }

  private fmtClock(d: Date): string {
    const pad = (n: number): string => n < 10 ? `0${n}` : `${n}`
    return `${d.getMonth() + 1}/${d.getDate()} ` +
           `${pad(d.getHours())}:${pad(d.getMinutes())}:${pad(d.getSeconds())}`
  }

  loadMonthData() {
    try {
      this.monthDays = NativeBridge.getMonthData(this.currentYear, this.currentMonth)
      if (this.monthDays.length > 0) {
        const first = this.monthDays[0]
        this.yearGZ = first.yearGZ
        this.shengXiao = first.shengXiao
        this.huangdiYear = first.huangdiYear
      }
      this.yearInput = YearUtil.astroYearToStr(this.currentYear)
      this.monthInput = `${this.currentMonth}`
      this.collectMonthEvents()
      // 切月时 rtsDay 自动夹紧
      if (this.rtsDay > this.monthDays.length) this.rtsDay = 1
      this.loadRTS()
    } catch (e) {
      console.error('Failed to load month data: ' + JSON.stringify(e))
    }
  }

  // 计算/刷新当日日月升降
  loadRTS() {
    try {
      this.rts = NativeBridge.calcDayRTS(
        this.currentYear, this.currentMonth, this.rtsDay,
        this.location.longitude, this.location.latitude, this.location.timezone)
    } catch (e) {
      console.error('calcDayRTS failed: ' + JSON.stringify(e))
      this.rts = null
    }
  }

  // 从当月 DayInfo 派生: 朔/望/上弦/下弦, 以及 24 节气
  collectMonthEvents() {
    const moon: MonthEvent[] = []
    const jq: MonthEvent[] = []
    for (let i = 0; i < this.monthDays.length; i++) {
      const d = this.monthDays[i]
      const yueName = d.yueXiangName
      if (yueName !== '') {
        const yt = this.extractTime(d.yueXiangTime)
        // JS 月历命名: 朔→朔月, 望→望月
        let label = yueName
        if (yueName === '朔') label = '朔月'
        else if (yueName === '望') label = '望月'
        moon.push({ day: d.solarDay, time: yt, name: label })
      }
      if (d.jieQiName !== '') {
        jq.push({ day: d.solarDay, time: this.extractTime(d.jieQiTime), name: d.jieQiName })
      }
    }
    this.moonEvents = moon
    this.jieQiEvents = jq
  }

  navigate(deltaYear: number, deltaMonth: number) {
    let y = this.currentYear + deltaYear
    let m = this.currentMonth + deltaMonth
    while (m <= 0)  { m += 12; y-- }
    while (m > 12) { m -= 12; y++ }
    this.currentYear = y
    this.currentMonth = m
    this.popupKey = ''
    this.rtsDay = 1
    this.loadMonthData()
  }

  goToday() {
    const now = new Date()
    this.currentYear = now.getFullYear()
    this.currentMonth = now.getMonth() + 1
    this.todayYear = this.currentYear
    this.todayMonth = this.currentMonth
    this.todayDay = now.getDate()
    this.rtsDay = this.todayDay
    this.popupKey = ''
    this.loadMonthData()
  }

  // 应用输入框中的年/月跳转 (年支持: 2026 / -211 / B212 / *212)
  applyYearMonthInput() {
    const y = YearUtil.yearStrToAstro(this.yearInput)
    const m = parseInt(this.monthInput, 10)
    if (YearUtil.isAstroYearValid(y) && !isNaN(m) && m >= 1 && m <= 12) {
      this.currentYear = y
      this.currentMonth = m
      this.popupKey = ''
      this.rtsDay = 1
      this.loadMonthData()
    } else {
      // 输入非法, 恢复显示当前值
      this.yearInput = YearUtil.astroYearToStr(this.currentYear)
      this.monthInput = `${this.currentMonth}`
    }
  }

  dayKey(day: DayInfo): string {
    return `${day.solarYear}-${day.solarMonth}-${day.solarDay}`
  }

  isToday(day: DayInfo): boolean {
    return day.solarYear === this.todayYear &&
      day.solarMonth === this.todayMonth &&
      day.solarDay === this.todayDay
  }

  isPopupOpen(day: DayInfo): boolean {
    return this.popupKey === this.dayKey(day)
  }

  getFirstDayOffset(): number {
    if (this.monthDays.length === 0) return 0
    return this.monthDays[0].weekDay
  }

  // 摘要文本: 节日 > 节气 > 朔望初一 > 农历日名
  getCellSubText(day: DayInfo): string {
    if (day.holiday && day.holiday.length > 0) return this.shortName(day.holiday)
    if (day.major && day.major.length > 0)     return this.shortName(day.major)
    if (day.jieQiName && day.jieQiName.length > 0) return day.jieQiName
    if (day.lunarDayName === '初一') return day.lunarMonthName
    return day.lunarDayName
  }

  shortName(s: string): string {
    let first = s.split(' ')[0].split(',')[0].split('，')[0]
    if (first.length > 4) first = first.substring(0, 4)
    return first
  }

  extractTime(s: string): string {
    const m = s.match(/(\d{2}:\d{2}:\d{2})/)
    return m !== null ? m[1] : s.trim()
  }

  // 儒略日格式化(与 JS 版 index.htm 一致): "完整JD(J2000起算相对值)"
  // JS: 'JD '+(ob.d0+J2000)+'('+ob.d0+')'
  formatJD(jd: number): string {
    if (jd === undefined || jd === null || isNaN(jd)) return '-'
    const fullJD = Math.round(jd)
    const d0 = fullJD - 2451545
    return `${fullJD}(${d0})`
  }

  build() {
    Column() {
      this.headerSection()
      this.navSection()
      this.topLocationStrip()
      this.yearInfoBar()
      this.weekHeader()
      this.calendarGrid()
      this.bottomInfoBar()
    }
    .width('100%')
    .height('100%')
    .backgroundColor(AppColors.BACKGROUND)
    // ── 底部抽屉: 老黄历 / 董公择日经典知识 共用 ──
    //   通过 sheetMode 切换内容; 单个 bindSheet 避免多 sheet 互相覆盖.
    .bindSheet($$this.sheetVisible, this.unifiedSheetBuilder(), {
      height: SheetSize.LARGE,
      dragBar: true,
      title: { title: this.sheetTitle() },
      onDisappear: (): void => { this.closeAlmanacSheet() }
    })
  }

  sheetTitle(): string {
    return this.sheetMode === 'topics' ? '董公择日 · 经典知识' : '老黄历'
  }

  @Builder
  unifiedSheetBuilder() {
    if (this.sheetMode === 'topics') {
      AlmanacTopicsSheet({ topics: this.topics })
    } else {
      AlmanacSheet({
        day: this.sheetDay,
        almanac: this.sheetAlmanac,
        onOpenTopics: (): void => { this.openTopicsSheet() }
      })
    }
  }

  // ── 顶部标题 ────────────────────────────────────────────
  @Builder
  headerSection() {
    Row() {
      Text(`${YearUtil.astroYearToStr(this.currentYear)}年${this.currentMonth}月`)
        .fontSize(AppDimens.FONT_TITLE)
        .fontWeight(FontWeight.Bold)
        .fontColor(AppColors.ON_PRIMARY)
        .layoutWeight(1)
        .textAlign(TextAlign.Start)

      Row() {
        Text('农历')
          .fontSize(AppDimens.FONT_CAPTION)
          .fontColor(AppColors.ON_PRIMARY + 'B0')
        Text(`${this.yearGZ}年`)
          .fontSize(AppDimens.FONT_CAPTION)
          .fontColor(AppColors.ACCENT)
          .fontWeight(FontWeight.Medium)
          .margin({ left: 4 })
      }
    }
    .width('100%')
    .padding({
      left: AppDimens.PADDING_LG,
      right: AppDimens.PADDING_LG,
      top: AppDimens.PADDING_XL,
      bottom: AppDimens.PADDING_SM
    })
    .alignItems(VerticalAlign.Center)
    .linearGradient({
      direction: GradientDirection.Right,
      colors: [[AppColors.GRADIENT_START, 0.0], [AppColors.GRADIENT_END, 1.0]]
    })
  }

  // ── 导航条: 上下年/月 + 年月输入 + 今天 ────────────────
  @Builder
  navSection() {
    Row() {
      // 上一年
      this.navIconButton('«', () => { this.navigate(-1, 0) })
      // 上一月
      this.navIconButton('‹', () => { this.navigate(0, -1) })

      // 年输入: 接受 2026 / -211 / B212 / *212 等格式
      TextInput({ text: this.yearInput, placeholder: 'YYYY/B212' })
        .type(InputType.Normal)
        .maxLength(6)
        .textAlign(TextAlign.Center)
        .fontSize(AppDimens.FONT_BODY)
        .fontColor(AppColors.ON_PRIMARY)
        .placeholderColor(AppColors.ON_PRIMARY + '60')
        .backgroundColor(AppColors.PRIMARY_LIGHT + '60')
        .borderRadius(AppDimens.BORDER_RADIUS_SM)
        .width(76)
        .height(32)
        .margin({ left: 6 })
        .onChange((v: string) => { this.yearInput = v })
        .onSubmit(() => { this.applyYearMonthInput() })

      Text('年')
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(AppColors.ON_PRIMARY)
        .margin({ left: 2, right: 2 })

      // 月输入
      TextInput({ text: this.monthInput })
        .type(InputType.Number)
        .maxLength(2)
        .textAlign(TextAlign.Center)
        .fontSize(AppDimens.FONT_BODY)
        .fontColor(AppColors.ON_PRIMARY)
        .backgroundColor(AppColors.PRIMARY_LIGHT + '60')
        .borderRadius(AppDimens.BORDER_RADIUS_SM)
        .width(40)
        .height(32)
        .onChange((v: string) => { this.monthInput = v })
        .onSubmit(() => { this.applyYearMonthInput() })

      Text('月')
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(AppColors.ON_PRIMARY)
        .margin({ left: 2, right: 6 })

      // 下一月
      this.navIconButton('›', () => { this.navigate(0, 1) })
      // 下一年
      this.navIconButton('»', () => { this.navigate(1, 0) })

      Blank()

      Button('今天')
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(AppColors.PRIMARY)
        .backgroundColor(AppColors.ACCENT)
        .borderRadius(AppDimens.BORDER_RADIUS_LG)
        .height(30)
        .padding({ left: 12, right: 12 })
        .onClick(() => { this.goToday() })
    }
    .width('100%')
    .padding({
      left: AppDimens.PADDING_LG,
      right: AppDimens.PADDING_LG,
      top: AppDimens.PADDING_SM,
      bottom: AppDimens.PADDING_MD
    })
    .alignItems(VerticalAlign.Center)
    .linearGradient({
      direction: GradientDirection.Right,
      colors: [[AppColors.GRADIENT_START, 0.0], [AppColors.GRADIENT_END, 1.0]]
    })
  }

  @Builder
  navIconButton(icon: string, onTap: () => void) {
    Button({ type: ButtonType.Circle }) {
      Text(icon)
        .fontSize(18)
        .fontWeight(FontWeight.Bold)
        .fontColor(AppColors.ON_PRIMARY)
    }
    .width(30)
    .height(30)
    .margin({ right: 4 })
    .backgroundColor(AppColors.PRIMARY_LIGHT + '80')
    .onClick(onTap)
  }

  @Builder
  yearInfoBar() {
    Row() {
      Text('☰')
        .fontSize(14)
        .fontColor(AppColors.ACCENT)
        .margin({ right: 6 })
      Text(`${this.yearGZ}年`)
        .fontSize(AppDimens.FONT_BODY)
        .fontColor(AppColors.PRIMARY)
        .fontWeight(FontWeight.Medium)
      Text('  |  ')
        .fontSize(AppDimens.FONT_BODY)
        .fontColor(AppColors.DIVIDER)
      Text(`生肖${this.shengXiao}`)
        .fontSize(AppDimens.FONT_BODY)
        .fontColor(AppColors.PRIMARY)
      Blank()
      if (this.huangdiYear > 0) {
        Text(`黄帝${this.huangdiYear}年`)
          .fontSize(AppDimens.FONT_CAPTION)
          .fontColor(AppColors.TEXT_SECONDARY)
      }
    }
    .width('100%')
    .padding({
      left: AppDimens.PADDING_LG,
      right: AppDimens.PADDING_LG,
      top: AppDimens.PADDING_SM,
      bottom: AppDimens.PADDING_SM
    })
    .backgroundColor(AppColors.SURFACE)
    .border({ width: { bottom: 0.5 }, color: AppColors.DIVIDER })
  }

  @Builder
  weekHeader() {
    Row() {
      ForEach(WEEK_NAMES, (name: string, index: number) => {
        Text(name)
          .fontSize(AppDimens.FONT_CAPTION)
          .fontColor(index === 0 || index === 6 ?
            AppColors.WEEKEND_COLOR : AppColors.TEXT_SECONDARY)
          .fontWeight(FontWeight.Medium)
          .width('14.28%')
          .textAlign(TextAlign.Center)
          .padding({ top: AppDimens.PADDING_SM, bottom: AppDimens.PADDING_SM })
      }, (name: string, index: number) => `wn-${index}`)
    }
    .width('100%')
    .backgroundColor(AppColors.SURFACE)
    .padding({ left: AppDimens.PADDING_SM, right: AppDimens.PADDING_SM })
  }

  @Builder
  calendarGrid() {
    Scroll() {
      Flex({ wrap: FlexWrap.Wrap }) {
        ForEach(Array.from({ length: this.getFirstDayOffset() }, (_: undefined, i: number) => i),
          (i: number) => {
            Column()
              .width('14.28%')
              .height(64)
          }, (i: number) => `blank-${i}`)

        ForEach(this.monthDays, (day: DayInfo) => {
          this.dayCell(day)
        }, (day: DayInfo): string => this.dayKey(day))
      }
      .width('100%')
      .padding({ left: AppDimens.PADDING_SM, right: AppDimens.PADDING_SM })
    }
    .layoutWeight(1)
    .scrollBar(BarState.Off)
    .backgroundColor(AppColors.SURFACE)
  }

  @Builder
  dayCell(day: DayInfo) {
    Column() {
      Row() {
        Text(`${day.solarDay}`)
          .fontSize(AppDimens.FONT_SUBTITLE)
          .fontWeight(this.isToday(day) ? FontWeight.Bold : FontWeight.Normal)
          .fontColor(this.getDaySolarColor(day))
        if (day.yueXiangName === '朔' || day.yueXiangName === '望') {
          Text('●')
            .fontSize(8)
            .fontColor(day.yueXiangName === '朔' ? '#505050' : '#F0B000')
            .margin({ left: 2 })
        } else if (day.yueXiangName === '上弦' || day.yueXiangName === '下弦') {
          Text(day.yueXiangName === '上弦' ? '◑' : '◐')
            .fontSize(8)
            .fontColor('#F0B000')
            .margin({ left: 2 })
        }
      }
      .alignItems(VerticalAlign.Center)

      Text(this.getCellSubText(day))
        .fontSize(AppDimens.FONT_SMALL)
        .fontColor(this.getDayLunarColor(day))
        .fontWeight(day.jieQiName !== '' || day.holiday !== '' ?
          FontWeight.Medium : FontWeight.Normal)
        .maxLines(1)
        .textOverflow({ overflow: TextOverflow.Clip })
    }
    .width('14.28%')
    .height(64)
    .justifyContent(FlexAlign.Center)
    .alignItems(HorizontalAlign.Center)
    .borderRadius(AppDimens.BORDER_RADIUS_SM)
    .backgroundColor(this.getDayCellBg(day))
    .border(this.isPopupOpen(day) && !this.isToday(day) ?
      { width: 1.5, color: AppColors.ACCENT } :
      { width: 0, color: Color.Transparent })
    .onClick(() => {
      const key = this.dayKey(day)
      this.popupKey = (this.popupKey === key) ? '' : key
      // 点击切换底部日月升降目标日
      this.rtsDay = day.solarDay
      this.loadRTS()
      if (this.popupKey !== '') this.ensurePopupAlmanac(day)
    })
    .onHover((isHover: boolean) => {
      if (isHover) {
        this.popupKey = this.dayKey(day)
        this.rtsDay = day.solarDay
        this.loadRTS()
        this.ensurePopupAlmanac(day)
      }
    })
    .bindPopup(this.isPopupOpen(day), {
      builder: () => { this.dayPopupContent(day) },
      placement: Placement.Top,
      autoCancel: true,
      enableArrow: true,
      popupColor: '#14163A',                       // 深午夜蓝, 与 Android/iOS 一致
      mask: { color: '#000000A0' },                // 深色遮罩, 凸显浮层
      radius: AppDimens.BORDER_RADIUS_MD,
      onStateChange: (e: PopupStateChangeParam) => {
        if (!e.isVisible && this.popupKey === this.dayKey(day)) {
          this.popupKey = ''
          this.resetRtsToToday()
        }
      }
    })
  }

  // 浮层关闭: 当前月若为今月 → 日月升降回今天, 否则 → 1 号
  resetRtsToToday() {
    const resetDay: number = (this.currentYear === this.todayYear
                          && this.currentMonth === this.todayMonth) ? this.todayDay : 1
    if (resetDay !== this.rtsDay) {
      this.rtsDay = resetDay
      this.loadRTS()
    }
  }

  // 浮层老黄历摘要: 仅在 popup 切换时调用一次, 避免每帧跨 JNI.
  //   key 用日期, 防止 popupKey 抖动时重复请求.
  ensurePopupAlmanac(day: DayInfo) {
    const key = this.dayKey(day)
    if (this.popupAlmanacKey === key && this.popupAlmanac !== null) return
    try {
      this.popupAlmanac = NativeBridge.getAlmanac(day.solarYear, day.solarMonth, day.solarDay)
      this.popupAlmanacKey = key
    } catch (e) {
      console.error('getAlmanac failed: ' + JSON.stringify(e))
      this.popupAlmanac = null
    }
  }

  // 打开底部抽屉显示完整老黄历
  openAlmanacSheet(day: DayInfo) {
    try {
      const a = NativeBridge.getAlmanac(day.solarYear, day.solarMonth, day.solarDay)
      if (a === null) return
      this.sheetDay = day
      this.sheetAlmanac = a
      this.sheetMode = 'almanac'
      this.sheetVisible = true
    } catch (e) {
      console.error('openAlmanacSheet failed: ' + JSON.stringify(e))
    }
  }

  closeAlmanacSheet() {
    this.sheetVisible = false
    this.sheetMode = 'almanac'
    this.sheetDay = null
    this.sheetAlmanac = null
  }

  // 切换到"董公择日静态知识"; 第一次打开懒加载, 后续命中缓存.
  //   sheet 已显示, 仅替换 mode -> 内容刷新, 不重新做出场动画
  openTopicsSheet() {
    try {
      if (this.topics.length === 0) {
        this.topics = NativeBridge.getAlmanacTopics()
      }
      this.sheetMode = 'topics'
    } catch (e) {
      console.error('openTopicsSheet failed: ' + JSON.stringify(e))
    }
  }

  getDaySolarColor(day: DayInfo): string {
    if (this.isToday(day)) return AppColors.ON_PRIMARY
    if (day.isOffDay || day.holiday !== '') return AppColors.WEEKEND_COLOR
    if (day.weekDay === 0 || day.weekDay === 6) return AppColors.WEEKEND_COLOR
    return AppColors.ON_SURFACE
  }

  getDayLunarColor(day: DayInfo): string {
    if (day.holiday !== '') return AppColors.WEEKEND_COLOR
    if (day.jieQiName !== '') return AppColors.JIEQI_COLOR
    if (this.isToday(day)) return AppColors.ON_PRIMARY + 'CC'
    return AppColors.LUNAR_TEXT
  }

  getDayCellBg(day: DayInfo): string {
    if (this.isToday(day)) return AppColors.PRIMARY
    if (this.isPopupOpen(day)) return AppColors.ACCENT + '20'
    return 'transparent'
  }

  // ─────────────────────────────────────────────────────────────
  //  浮层内容: 深色紧凑卡片 (与 Android / iOS 同款)
  //   - 容器宽 240~300, padding 12, 行间距 1px, 字号紧缩一档
  //   - 分隔线 margin 各 3, 拉近视觉密度
  //   - 标题居中 + "今天"徽章右贴; 三纳音加粗
  // ─────────────────────────────────────────────────────────────
  @Builder
  dayPopupContent(day: DayInfo) {
    Column() {
      // 1) 标题行: 公历 + 今天徽章 (徽章右贴)
      Row() {
        Text(`${YearUtil.astroYearToStr(day.solarYear)}年${day.solarMonth}月${day.solarDay}日`)
          .fontSize(AppDimens.FONT_SUBTITLE)
          .fontWeight(FontWeight.Bold)
          .fontColor('#FFFFFF')
          .layoutWeight(1)
          .textAlign(TextAlign.Center)
        if (this.isToday(day)) {
          Text('今天')
            .fontSize(AppDimens.FONT_SMALL)
            .fontColor(AppColors.PRIMARY)
            .fontWeight(FontWeight.Medium)
            .backgroundColor(AppColors.ACCENT)
            .borderRadius(4)
            .padding({ left: 5, right: 5, top: 1, bottom: 1 })
        }
      }
      .width('100%')
      .alignItems(VerticalAlign.Center)
      .margin({ bottom: 6 })

      // 2) 黄帝纪年 / 星期 / 星座
      this.popupInfoLine(`黄帝${day.huangdiYear}年`, day.weekName, day.constellationName, '#FFFFFF', false)
      // 3) 干支年 / 农历月 / 农历日
      this.popupInfoLine(
        `${day.yearGZ}年`,
        `${day.isLeapMonth ? '闰' : ''}${day.lunarMonthName}`,
        `${day.lunarDayName}日`,
        '#FFFFFF', false
      )
      // 4) 干支年月日 (三柱)
      this.popupInfoLine(`${day.yearGZ}年`, `${day.monthGZ}月`, `${day.dayGZ}日`, '#FFFFFF', false)
      // 5) 三柱纳音 (高亮+加粗)
      this.popupInfoLine(day.yearNaYin, day.monthNaYin, day.dayNaYin, '#FFD54F', true)

      Divider().color('#5A5F8C').margin({ top: 4, bottom: 3 })

      // 6) 生肖 · 建除 (居中)
      this.popupCenterRow('生肖', day.shengXiao, '建除', day.jian12Name)
      // 7) 回历 + JD (居中)
      this.popupCenterText(`回历[${day.moslemYear}年${day.moslemMonth}月${day.moslemDay}日]`,
        '#DCE0F2', AppDimens.FONT_SMALL)
      this.popupCenterText(`JD ${this.formatJD(day.julianDay)}`,
        '#DCE0F2', AppDimens.FONT_SMALL)

      // 8) 节气/月相/节日
      if (day.jieQiName !== '' || day.yueXiangName !== '' ||
          day.holiday !== '' || day.major !== '' ||
          day.minor !== '' || day.misc !== '') {
        Divider().color('#5A5F8C').margin({ top: 4, bottom: 3 })
        if (day.jieQiName !== '') {
          this.popupEventLine('🌿', `节气 ${day.jieQiName}`,
            day.jieQiTime !== '' ? this.extractTime(day.jieQiTime) : '', '#80E0A7')
        }
        if (day.yueXiangName !== '') {
          this.popupEventLine('🌙', `月相 ${day.yueXiangName}`,
            day.yueXiangTime !== '' ? this.extractTime(day.yueXiangTime) : '', '#FFD54F')
        }
        if (day.holiday !== '')  this.popupEventLine('🎉', day.holiday, '', '#FF9D9D')
        if (day.major !== '')    this.popupEventLine('🎊', day.major, '', '#FFD54F')
        if (day.minor !== '')    this.popupEventLine('·',  day.minor, '', '#DCE0F2')
        if (day.misc !== '')     this.popupEventLine('🍂', day.misc, '', '#80E0A7')
      }

      // 9) 老黄历摘要 + 详情入口
      Divider().color('#5A5F8C').margin({ top: 4, bottom: 3 })
      AlmanacInlineRow({ almanac: this.popupAlmanac })
      Row() {
        Text('📜 老黄历详情 ›')
          .fontSize(AppDimens.FONT_CAPTION)
          .fontColor(AppColors.ACCENT)
          .fontWeight(FontWeight.Medium)
      }
      .width('100%')
      .justifyContent(FlexAlign.Center)
      .padding({ top: 6, bottom: 2 })
      .onClick(() => { this.openAlmanacSheet(day) })
    }
    .padding({ left: 12, right: 12, top: 12, bottom: 12 })
    .constraintSize({ minWidth: 200, maxWidth: 240 })
    .backgroundColor('#14163A')
    .borderRadius(AppDimens.BORDER_RADIUS_MD)
    .border({ width: 0.5, color: '#3B3F73', radius: AppDimens.BORDER_RADIUS_MD })
    .alignItems(HorizontalAlign.Center)
  }

  // ── 居中紧凑三列行 (黄帝/星期/星座, 干支/农历, 三柱纳音) ──
  //  说明: 不再三等分整宽, 改为水平居中 + 小间隔, 视觉更紧凑
  @Builder
  popupInfoLine(a: string, b: string, c: string, color: string, bold: boolean) {
    Row() {
      this.popupCell(a, color, bold)
      this.popupCell(b, color, bold)
      this.popupCell(c, color, bold)
    }
    .width('100%')
    .justifyContent(FlexAlign.Center)
    .margin({ top: 1, bottom: 1 })
  }

  @Builder
  popupCell(text: string, color: string, bold: boolean) {
    Text(text && text.length > 0 ? text : '—')
      .fontSize(AppDimens.FONT_CAPTION)
      .fontWeight(bold ? FontWeight.Medium : FontWeight.Normal)
      .fontColor(text && text.length > 0 ? color : '#DCE0F2')
      .maxLines(1)
      .textAlign(TextAlign.Center)
      .textOverflow({ overflow: TextOverflow.Ellipsis })
      .margin({ left: 6, right: 6 })
  }

  @Builder
  popupCenterRow(l1: string, v1: string, l2: string, v2: string) {
    Row() {
      Text(l1).fontSize(AppDimens.FONT_SMALL).fontColor('#DCE0F2').margin({ right: 4 })
      Text(v1).fontSize(AppDimens.FONT_CAPTION).fontColor('#FFFFFF')
        .fontWeight(FontWeight.Medium)
      Text(' · ').fontSize(AppDimens.FONT_CAPTION).fontColor('#DCE0F2')
        .margin({ left: 6, right: 6 })
      Text(l2).fontSize(AppDimens.FONT_SMALL).fontColor('#DCE0F2').margin({ right: 4 })
      Text(v2).fontSize(AppDimens.FONT_CAPTION).fontColor('#FFFFFF')
        .fontWeight(FontWeight.Medium)
    }
    .width('100%')
    .justifyContent(FlexAlign.Center)
    .margin({ top: 1, bottom: 1 })
  }

  @Builder
  popupCenterText(text: string, color: string, size: number) {
    Row() {
      Text(text)
        .fontSize(size)
        .fontColor(color)
        .maxLines(1)
        .textAlign(TextAlign.Center)
    }
    .width('100%')
    .justifyContent(FlexAlign.Center)
    .margin({ top: 1, bottom: 1 })
  }

  @Builder
  popupEventLine(icon: string, text: string, time: string, color: string) {
    Row() {
      Text(icon).fontSize(11).fontColor(color).margin({ right: 4 })
      Text(text)
        .fontSize(AppDimens.FONT_CAPTION).fontColor(color)
        .fontWeight(FontWeight.Medium)
        .maxLines(2)
        .textOverflow({ overflow: TextOverflow.Ellipsis })
      if (time !== '') {
        Text(` ${time}`).fontSize(AppDimens.FONT_SMALL).fontColor('#DCE0F2').margin({ left: 4 })
      }
    }
    .width('100%')
    .justifyContent(FlexAlign.Center)
    .alignItems(VerticalAlign.Center)
    .margin({ top: 1, bottom: 1 })
  }

  // ─────────────────────────────────────────────────────────────
  //  顶部地点条 — 两行内联 Select, 对齐 JS 原版 indexmp.htm 顶部
  //
  //    第 1 行: [洲 ▼] [国家 ▼]  · UTC±x  当地时间    (国际, 仅展示)
  //    第 2 行: [省 ▼] [市/县 ▼] · 经X°E 纬Y°N         (国内, 改 location)
  //
  //  · 国内行决定日月升降算用的经纬度
  //  · 国际行只影响顶栏外地时钟; 日历主时区始终用系统本地时区, 与 JS 一致
  // ─────────────────────────────────────────────────────────────
  @Builder
  topLocationStrip() {
    Column() {
      // ── 国际行 (展示外地时区+实时时间) ────────────
      Row() {
        if (this.continents.length > 0) {
          this.continentSelect()
          this.countrySelect()
        }
        Blank()
        Text(this.countries.length > 0 ? formatTimezone(this.countries[this.countryIdx].timezone) : '')
          .fontSize(AppDimens.FONT_SMALL)
          .fontColor(AppColors.PRIMARY)
          .margin({ left: 6, right: 6 })
        Text(this.intlClock)
          .fontSize(AppDimens.FONT_SMALL)
          .fontColor(AppColors.ACCENT)
          .fontWeight(FontWeight.Medium)
      }
      .width('100%')
      .padding({ left: 8, right: AppDimens.PADDING_MD, top: 4, bottom: 4 })

      // ── 国内行 (改 location, 触发日月升降重算) ────
      Row() {
        if (this.provinces.length > 0) {
          this.provinceSelect()
          this.citySelect()
        }
        Blank()
        Text(formatLonLat(this.location.longitude, this.location.latitude))
          .fontSize(AppDimens.FONT_SMALL)
          .fontColor(AppColors.TEXT_SECONDARY)
      }
      .width('100%')
      .padding({ left: 8, right: AppDimens.PADDING_MD, top: 4, bottom: 4 })
    }
    .width('100%')
    .backgroundColor(AppColors.SURFACE)
    .border({ width: { bottom: 0.5 }, color: AppColors.DIVIDER })
  }

  @Builder
  continentSelect() {
    Select(this.contOpts)
      .selected(this.contIdx)
      .value(this.continents[this.contIdx])
      .font({ size: AppDimens.FONT_SMALL })
      .fontColor(AppColors.ON_SURFACE)
      .selectedOptionFont({ size: AppDimens.FONT_SMALL, weight: FontWeight.Medium })
      .optionFont({ size: AppDimens.FONT_SMALL })
      .onSelect((idx: number) => { this.onContinentChange(idx) })
  }

  @Builder
  countrySelect() {
    Select(this.countryOpts)
      .selected(this.countryIdx)
      .value(this.countries.length > 0 ? this.countries[this.countryIdx].country : '')
      .font({ size: AppDimens.FONT_SMALL })
      .fontColor(AppColors.ON_SURFACE)
      .selectedOptionFont({ size: AppDimens.FONT_SMALL, weight: FontWeight.Medium })
      .optionFont({ size: AppDimens.FONT_SMALL })
      .onSelect((idx: number) => { this.onCountryChange(idx) })
  }

  @Builder
  provinceSelect() {
    Select(this.provOpts)
      .selected(this.provIdx)
      .value(this.provinces[this.provIdx].province)
      .font({ size: AppDimens.FONT_SMALL })
      .fontColor(AppColors.ON_SURFACE)
      .selectedOptionFont({ size: AppDimens.FONT_SMALL, weight: FontWeight.Medium })
      .optionFont({ size: AppDimens.FONT_SMALL })
      .onSelect((idx: number) => { this.onProvinceChange(idx) })
  }

  @Builder
  citySelect() {
    Select(this.cityOpts)
      .selected(this.cityIdx)
      .value(this.cities.length > 0 ? this.cities[this.cityIdx].district : '')
      .font({ size: AppDimens.FONT_SMALL })
      .fontColor(AppColors.PRIMARY)
      .selectedOptionFont({ size: AppDimens.FONT_SMALL, weight: FontWeight.Medium })
      .optionFont({ size: AppDimens.FONT_SMALL })
      .onSelect((idx: number) => { this.onCityChange(idx) })
  }

  // ─────────────────────────────────────────────────────────────
  //  底部信息条 — 紧凑布局, 不再占用过多高度, 保证日历主体完整可见
  //
  //  · 单行日月时刻 (1 行 6 列: 日出/中天/日落/月出/月中/月落)
  //  · 1 行 (天亮/白昼/天黑)
  //  · 月相 + 节气 (各 1 行流式)
  // ─────────────────────────────────────────────────────────────
  @Builder
  bottomInfoBar() {
    Column() {
      // 日月升降 紧凑两行
      this.sunMoonCompact()
      // 月相 / 节气 — 一行流式排列
      if (this.moonEvents.length > 0) {
        this.eventLine(`月相`, this.moonEvents, AppColors.PRIMARY_LIGHT)
      }
      if (this.jieQiEvents.length > 0) {
        this.eventLine(`节气`, this.jieQiEvents, AppColors.JIEQI_COLOR)
      }
    }
    .width('100%')
    .padding({
      left: AppDimens.PADDING_MD,
      right: AppDimens.PADDING_MD,
      top: 6,
      bottom: 6,
    })
    .backgroundColor(AppColors.SURFACE)
    .border({ width: { top: 0.5 }, color: AppColors.DIVIDER })
  }

  // 日月升降紧凑卡 — 标题 + 一行 6 列 + 一行 3 列
  @Builder
  sunMoonCompact() {
    // 标题
    Row() {
      Text(`${this.currentMonth}月${this.rtsDay}日 · 日月升降`)
        .fontSize(AppDimens.FONT_SMALL)
        .fontColor(AppColors.PRIMARY)
        .fontWeight(FontWeight.Medium)
      Blank()
      Text(this.localClock)
        .fontSize(10)
        .fontColor(AppColors.TEXT_SECONDARY)
    }
    .width('100%')
    .margin({ bottom: 2 })

    // 一行 6 列 (太阳出/中/落 + 月亮出/中/落)
    Row() {
      this.tinyCell('☀↑', this.rts !== null ? this.rts!.sunRise     : '--:--:--')
      this.tinyCell('☀',  this.rts !== null ? this.rts!.sunMeridian : '--:--:--')
      this.tinyCell('☀↓', this.rts !== null ? this.rts!.sunSet      : '--:--:--')
      this.tinyCell('☾↑', this.rts !== null ? this.rts!.moonRise     : '--:--:--')
      this.tinyCell('☾',  this.rts !== null ? this.rts!.moonMeridian : '--:--:--')
      this.tinyCell('☾↓', this.rts !== null ? this.rts!.moonSet      : '--:--:--')
    }
    .width('100%')
    .margin({ top: 2, bottom: 2 })

    // 一行 3 列 (晨/昏/白昼)
    Row() {
      this.tinyCell('晨', this.rts !== null ? this.rts!.civilDawn   : '--:--:--')
      this.tinyCell('日', this.rts !== null ? this.rts!.dayLength   : '--:--:--')
      this.tinyCell('昏', this.rts !== null ? this.rts!.civilDusk   : '--:--:--')
      this.tinyCell('白', this.rts !== null ? this.rts!.lightLength : '--:--:--')
    }
    .width('100%')
    .margin({ top: 2 })
  }

  // 6/3 列等宽紧凑单元 (label + value)
  @Builder
  tinyCell(label: string, value: string) {
    Row() {
      Text(label)
        .fontSize(10)
        .fontColor(AppColors.TEXT_SECONDARY)
        .margin({ right: 3 })
      Text(value)
        .fontSize(10)
        .fontColor(value === '--:--:--' ? AppColors.TEXT_SECONDARY : AppColors.ON_SURFACE)
        .fontWeight(FontWeight.Medium)
    }
    .layoutWeight(1)
    .justifyContent(FlexAlign.Center)
  }

  // 月相/节气一行流式
  @Builder
  eventLine(title: string, events: MonthEvent[], dotColor: string) {
    Flex({ wrap: FlexWrap.Wrap, alignItems: ItemAlign.Center }) {
      Text(`${title}: `)
        .fontSize(10)
        .fontColor(AppColors.PRIMARY)
        .fontWeight(FontWeight.Medium)
      ForEach(events, (e: MonthEvent, idx: number) => {
        Row() {
          Text('●')
            .fontSize(8)
            .fontColor(dotColor)
            .margin({ right: 2 })
          Text(`${e.name}${this.currentMonth}/${e.day}`)
            .fontSize(10)
            .fontColor(AppColors.ON_SURFACE)
        }
        .margin({ right: 8 })
      }, (e: MonthEvent, idx: number): string => `${title}-${e.name}-${e.day}`)
    }
    .width('100%')
    .margin({ top: 4 })
  }


}
