import { AppColors, AppDimens, WEEK_NAMES } from '../common/Constants'
import { YearUtil } from '../common/YearUtil'
import { NativeBridge, DayInfo, DayRTS } from '../services/NativeBridge'

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

  // 站点参数 (默认北京)
  private siteLon: number = 116.3833
  private siteLat: number = 39.9
  private siteTz: number = 8

  aboutToAppear() {
    this.loadMonthData()
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
        this.siteLon, this.siteLat, this.siteTz)
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
      this.yearInfoBar()
      this.weekHeader()
      this.calendarGrid()
      this.bottomInfoBar()
    }
    .width('100%')
    .height('100%')
    .backgroundColor(AppColors.BACKGROUND)
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
    })
    .onHover((isHover: boolean) => {
      if (isHover) {
        this.popupKey = this.dayKey(day)
        this.rtsDay = day.solarDay
        this.loadRTS()
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
  //  月历底部信息栏: 当日日月升降 (跟随选中日) + 当月朔望/节气 (按月)
  //  对应 sxwnl/lunar.js 的 RTS1() 文字模块
  // ─────────────────────────────────────────────────────────────
  @Builder
  bottomInfoBar() {
    Column() {
      // 标题: 选中日 (北京)
      Row() {
        Text(`${YearUtil.astroYearToStr(this.currentYear)}年${this.currentMonth}月${this.rtsDay}日 · 日月升降`)
          .fontSize(AppDimens.FONT_CAPTION)
          .fontColor(AppColors.PRIMARY)
          .fontWeight(FontWeight.Medium)
        Blank()
        Text(`北京 (${this.siteLon.toFixed(2)}°E ${this.siteLat.toFixed(2)}°N)`)
          .fontSize(AppDimens.FONT_SMALL)
          .fontColor(AppColors.TEXT_SECONDARY)
      }
      .width('100%')
      .margin({ bottom: 6 })

      // 日出/日落/中天
      this.rtsTripleRow('日出',     this.rts !== null ? this.rts!.sunRise     : '--:--:--',
                        '日落',     this.rts !== null ? this.rts!.sunSet      : '--:--:--',
                        '中天',     this.rts !== null ? this.rts!.sunMeridian : '--:--:--')
      this.rtsTripleRow('月出',     this.rts !== null ? this.rts!.moonRise     : '--:--:--',
                        '月落',     this.rts !== null ? this.rts!.moonSet      : '--:--:--',
                        '月中',     this.rts !== null ? this.rts!.moonMeridian : '--:--:--')
      this.rtsDoubleRow('晨起天亮', this.rts !== null ? this.rts!.civilDawn   : '--:--:--',
                        '晚上天黑', this.rts !== null ? this.rts!.civilDusk   : '--:--:--')
      this.rtsDoubleRow('日照时间', this.rts !== null ? this.rts!.dayLength   : '--:--:--',
                        '白天时间', this.rts !== null ? this.rts!.lightLength : '--:--:--')

      // 月相 / 节气 — 按月固定不变
      if (this.moonEvents.length > 0 || this.jieQiEvents.length > 0) {
        Divider().color(AppColors.DIVIDER).margin({ top: 8, bottom: 6 })

        Text(`${this.currentMonth}月月相与节气`)
          .fontSize(AppDimens.FONT_SMALL)
          .fontColor(AppColors.TEXT_SECONDARY)
          .margin({ bottom: 4 })

        // 月相 (按 2 列布局)
        this.eventGrid(this.moonEvents, AppColors.PRIMARY_LIGHT)

        if (this.jieQiEvents.length > 0) {
          // 节气分隔与列表
          this.eventGrid(this.jieQiEvents, AppColors.JIEQI_COLOR)
        }
      }
    }
    .width('100%')
    .padding(AppDimens.PADDING_MD)
    .backgroundColor(AppColors.SURFACE)
    .border({ width: { top: 0.5 }, color: AppColors.DIVIDER })
  }

  // 三列等宽 (用于"日出/日落/中天"这种紧凑信息)
  @Builder
  rtsTripleRow(l1: string, v1: string, l2: string, v2: string, l3: string, v3: string) {
    Row() {
      this.rtsLabeled(l1, v1)
      this.rtsLabeled(l2, v2)
      this.rtsLabeled(l3, v3)
    }
    .width('100%')
    .margin({ top: 2, bottom: 2 })
  }

  // 两列等宽 (用于"晨起天亮/晚上天黑"这种)
  @Builder
  rtsDoubleRow(l1: string, v1: string, l2: string, v2: string) {
    Row() {
      this.rtsLabeled(l1, v1)
      this.rtsLabeled(l2, v2)
    }
    .width('100%')
    .margin({ top: 2, bottom: 2 })
  }

  @Builder
  rtsLabeled(label: string, value: string) {
    Row() {
      Text(label)
        .fontSize(AppDimens.FONT_SMALL)
        .fontColor(AppColors.TEXT_SECONDARY)
        .margin({ right: 4 })
      Text(value)
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(value === '--:--:--' ? AppColors.TEXT_SECONDARY : AppColors.ON_SURFACE)
        .fontWeight(FontWeight.Medium)
    }
    .layoutWeight(1)
    .alignItems(VerticalAlign.Center)
  }

  // 月相/节气栅格 (两列, JS 月历下方"3日 09:17:44 朔月 / 18日 16:47:28 望月"风格)
  @Builder
  eventGrid(events: MonthEvent[], color: string) {
    Flex({ wrap: FlexWrap.Wrap }) {
      ForEach(events, (e: MonthEvent, idx: number) => {
        Row() {
          Text(this.padNum(e.day, 2) + '日')
            .fontSize(AppDimens.FONT_SMALL)
            .fontColor(AppColors.TEXT_SECONDARY)
            .width(28)
          Text(e.time)
            .fontSize(AppDimens.FONT_SMALL)
            .fontColor(AppColors.ON_SURFACE)
            .margin({ left: 2, right: 4 })
          Text(e.name)
            .fontSize(AppDimens.FONT_SMALL)
            .fontColor(color)
            .fontWeight(FontWeight.Medium)
        }
        .width('50%')
        .padding({ top: 2, bottom: 2 })
        .alignItems(VerticalAlign.Center)
      }, (e: MonthEvent, idx: number): string => `${e.name}-${e.day}`)
    }
    .width('100%')
  }

  padNum(n: number, width: number): string {
    const s = '' + n
    return s.length >= width ? s : ' '.repeat(width - s.length) + s
  }
}
