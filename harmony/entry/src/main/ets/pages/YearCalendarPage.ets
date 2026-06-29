import { AppColors, AppDimens } from '../common/Constants'
import { YearUtil } from '../common/YearUtil'
import { NativeBridge, YearCalMonth, YearCalJieQi } from '../services/NativeBridge'

// ────────────────────────────────────────────────────────────
//  年历页面: 按农历月聚合的一年数据
//  对应 sxwnl/lunar.js 中的 nianLi2HTML(年版"时刻/干支"模式):
//    每行展示一个农历月: [闰?][月名][大/小][朔日干支][朔日公历日期]
//                       该月内出现的节气: [日序][干支][节气名][日期]
// ────────────────────────────────────────────────────────────
@Component
export struct YearCalendarPage {
  @State currentYear: number = new Date().getFullYear()
  @State yearInput: string = YearUtil.astroYearToStr(new Date().getFullYear())
  @State months: YearCalMonth[] = []
  @State loading: boolean = false

  aboutToAppear() {
    this.reload()
  }

  reload() {
    this.loading = true
    try {
      this.months = NativeBridge.getYearCalendar(this.currentYear)
    } catch (e) {
      console.error('YearCalendar load failed: ' + JSON.stringify(e))
      this.months = []
    }
    this.loading = false
    this.yearInput = YearUtil.astroYearToStr(this.currentYear)
  }

  // 年份位移: delta 可为 -10/-1/+1/+10 等
  shiftYear(delta: number) {
    this.currentYear += delta
    this.reload()
  }

  goCurYear() {
    this.currentYear = new Date().getFullYear()
    this.reload()
  }

  // 应用输入框中的年份跳转 (支持: 2026 / -211 / B212 / *212)
  applyYearInput() {
    const y = YearUtil.yearStrToAstro(this.yearInput)
    if (YearUtil.isAstroYearValid(y)) {
      this.currentYear = y
      this.reload()
    } else {
      this.yearInput = YearUtil.astroYearToStr(this.currentYear)
    }
  }

  build() {
    Column() {
      this.headerSection()
      this.navSection()
      this.summaryBar()

      if (this.loading) {
        Column() {
          LoadingProgress().width(40).height(40).color(AppColors.PRIMARY)
          Text('正在计算 ' + YearUtil.astroYearToStr(this.currentYear) + ' 年历...')
            .fontSize(AppDimens.FONT_CAPTION)
            .fontColor(AppColors.TEXT_SECONDARY)
            .margin({ top: 8 })
        }
        .width('100%').layoutWeight(1)
        .justifyContent(FlexAlign.Center)
      } else {
        Scroll() {
          Column() {
            ForEach(this.months, (m: YearCalMonth, idx: number) => {
              this.monthRow(m, idx)
            })

            Text('数据来源: 寿星万年历内核(C++) · 节气精确到秒')
              .fontSize(AppDimens.FONT_SMALL)
              .fontColor(AppColors.TEXT_SECONDARY)
              .margin({ top: AppDimens.PADDING_LG, bottom: AppDimens.PADDING_LG })
              .textAlign(TextAlign.Center)
              .width('100%')
          }
          .width('100%')
          .padding({
            left: AppDimens.PADDING_MD,
            right: AppDimens.PADDING_MD,
            top: AppDimens.PADDING_MD
          })
        }
        .layoutWeight(1)
        .scrollBar(BarState.Auto)
      }
    }
    .width('100%')
    .height('100%')
    .backgroundColor(AppColors.BACKGROUND)
  }

  // ── 顶部标题 ────────────────────────────────────────────
  @Builder
  headerSection() {
    Row() {
      Column() {
        Text(`${YearUtil.astroYearToStr(this.currentYear)} 年`)
          .fontSize(AppDimens.FONT_TITLE)
          .fontWeight(FontWeight.Bold)
          .fontColor(AppColors.ON_PRIMARY)
        Text('农历年视图 (按朔月排列)')
          .fontSize(AppDimens.FONT_SMALL)
          .fontColor(AppColors.ON_PRIMARY + 'B0')
      }
      .alignItems(HorizontalAlign.Start)
      .layoutWeight(1)
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

  // ── 导航条: -10/-1/+1/+10 年 + 年输入 + 本年 ─────────
  @Builder
  navSection() {
    Row() {
      this.navIconButton('«', () => { this.shiftYear(-10) })
      this.navIconButton('‹', () => { this.shiftYear(-1) })

      TextInput({ text: this.yearInput, placeholder: 'YYYY/B212' })
        .type(InputType.Normal)
        .maxLength(6)
        .textAlign(TextAlign.Center)
        .fontSize(AppDimens.FONT_BODY)
        .fontColor(AppColors.ON_PRIMARY)
        .placeholderColor(AppColors.ON_PRIMARY + '60')
        .backgroundColor(AppColors.PRIMARY_LIGHT + '60')
        .borderRadius(AppDimens.BORDER_RADIUS_SM)
        .width(80)
        .height(32)
        .margin({ left: 6 })
        .onChange((v: string) => { this.yearInput = v })
        .onSubmit(() => { this.applyYearInput() })

      Text('年')
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(AppColors.ON_PRIMARY)
        .margin({ left: 2, right: 6 })

      this.navIconButton('›', () => { this.shiftYear(1) })
      this.navIconButton('»', () => { this.shiftYear(10) })

      Blank()

      Button('本年')
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(AppColors.PRIMARY)
        .backgroundColor(AppColors.ACCENT)
        .borderRadius(AppDimens.BORDER_RADIUS_LG)
        .height(30)
        .padding({ left: 12, right: 12 })
        .onClick(() => { this.goCurYear() })
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
  summaryBar() {
    Row() {
      Text(`本年共 ${this.months.length} 个农历月`)
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(AppColors.TEXT_SECONDARY)
      if (this.hasLeapMonth()) {
        Text(`  · 含 ${this.leapMonthName()}`)
          .fontSize(AppDimens.FONT_CAPTION)
          .fontColor(AppColors.JIEQI_COLOR)
          .fontWeight(FontWeight.Medium)
      }
      Blank()
      Text(`24 节气 / ${this.totalJieQi()} 项`)
        .fontSize(AppDimens.FONT_SMALL)
        .fontColor(AppColors.TEXT_SECONDARY)
    }
    .width('100%')
    .padding({
      left: AppDimens.PADDING_LG, right: AppDimens.PADDING_LG,
      top: AppDimens.PADDING_SM,  bottom: AppDimens.PADDING_SM
    })
    .backgroundColor(AppColors.SURFACE)
    .border({ width: { bottom: 0.5 }, color: AppColors.DIVIDER })
  }

  hasLeapMonth(): boolean {
    return this.months.some((m: YearCalMonth): boolean => m.isLeap)
  }
  leapMonthName(): string {
    const leap = this.months.find((m: YearCalMonth): boolean => m.isLeap)
    return leap !== undefined ? leap.monthName : ''
  }
  totalJieQi(): number {
    let n = 0
    for (let i = 0; i < this.months.length; i++) n += this.months[i].jieQi.length
    return n
  }

  // 单月一行
  @Builder
  monthRow(m: YearCalMonth, idx: number) {
    Column() {
      // 标题行: 月名 + 大小 + 朔日干支 + 朔日公历日期
      Row() {
        Text(this.formatMonthLabel(m))
          .fontSize(AppDimens.FONT_SUBTITLE)
          .fontWeight(FontWeight.Bold)
          .fontColor(m.isLeap ? AppColors.JIEQI_COLOR : AppColors.PRIMARY)
          .width(70)

        Column() {
          Text(`朔 ${m.shuoGZ}`)
            .fontSize(AppDimens.FONT_BODY)
            .fontColor(AppColors.ON_SURFACE)
            .fontWeight(FontWeight.Medium)
          Text(this.formatShuoDate(m))
            .fontSize(AppDimens.FONT_SMALL)
            .fontColor(AppColors.TEXT_SECONDARY)
        }
        .alignItems(HorizontalAlign.Start)
        .margin({ left: AppDimens.PADDING_MD })

        Blank()

        // 月大小标签
        Text(m.dayCount > 29 ? '大' : '小')
          .fontSize(AppDimens.FONT_CAPTION)
          .fontColor(AppColors.ON_PRIMARY)
          .backgroundColor(m.dayCount > 29 ? AppColors.PRIMARY : AppColors.TEXT_SECONDARY)
          .padding({ left: 6, right: 6, top: 2, bottom: 2 })
          .borderRadius(4)
          .margin({ right: 4 })

        Text(`${m.dayCount} 天`)
          .fontSize(AppDimens.FONT_SMALL)
          .fontColor(AppColors.TEXT_SECONDARY)
      }
      .width('100%')
      .alignItems(VerticalAlign.Center)
      .padding({ top: AppDimens.PADDING_MD, bottom: AppDimens.PADDING_SM })

      // 节气列表
      if (m.jieQi && m.jieQi.length > 0) {
        Column() {
          ForEach(m.jieQi, (jq: YearCalJieQi) => {
            this.jieQiRow(jq, m)
          })
        }
        .width('100%')
        .padding({
          left: AppDimens.PADDING_MD,
          right: AppDimens.PADDING_MD,
          top: AppDimens.PADDING_SM,
          bottom: AppDimens.PADDING_SM
        })
        .backgroundColor(AppColors.BACKGROUND)
        .borderRadius(AppDimens.BORDER_RADIUS_SM)
        .margin({ top: AppDimens.PADDING_XS })
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
    .borderRadius(AppDimens.BORDER_RADIUS_MD)
    .margin({ bottom: AppDimens.PADDING_SM })
  }

  @Builder
  jieQiRow(jq: YearCalJieQi, m: YearCalMonth) {
    Row() {
      Text(jq.dayName)
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(AppColors.TEXT_SECONDARY)
        .width(36)

      Text(jq.gz)
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(AppColors.ON_SURFACE)
        .fontWeight(FontWeight.Medium)
        .width(40)

      Text(jq.name)
        .fontSize(AppDimens.FONT_BODY)
        .fontColor(AppColors.JIEQI_COLOR)
        .fontWeight(FontWeight.Medium)
        .width(48)

      Text(this.formatJieQiDate(jq, m))
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(AppColors.TEXT_SECONDARY)
        .layoutWeight(1)
    }
    .width('100%')
    .padding({ top: 2, bottom: 2 })
  }

  formatMonthLabel(m: YearCalMonth): string {
    return m.monthName
  }

  formatShuoDate(m: YearCalMonth): string {
    return `${YearUtil.astroYearToStr(m.solarYear)}-${this.pad2(m.solarMonth)}-${this.pad2(m.solarDay)}`
  }

  formatJieQiDate(jq: YearCalJieQi, m: YearCalMonth): string {
    const y = YearUtil.astroYearToStr(m.solarYear)
    const date = `${this.pad2(jq.solarMonth)}-${this.pad2(jq.solarDay)}`
    return jq.time !== '' ? `${date} ${jq.time}` : `${y}-${date}`
  }

  pad2(n: number): string {
    return n < 10 ? '0' + n : '' + n
  }
}
