import { AppColors, AppDimens } from '../common/Constants'
import { NativeBridge, LunarMonth, ReverseItem } from '../services/NativeBridge'
import { BaziCalc, DateSelection, MODE_SOLAR, MODE_LUNAR, MODE_REVERSE } from '../model/BaziInput'

/**
 * 出生时间选择弹窗(底部抽屉)。
 * 三个 Tab: 公历 / 农历 / 四柱反推。
 * 选择完成后通过 onConfirm 回填一条 DateSelection 记录(反推结果归一为公历)。
 */
@CustomDialog
export struct DateTimePicker {
  controller: CustomDialogController
  init: DateSelection = {
    inputMode: MODE_SOLAR, year: 1990, month: 1, day: 1,
    hour: 12, minute: 0, isLeap: false, isSpec: false
  }
  onConfirm: (sel: DateSelection) => void = () => {}

  @State inputMode: number = MODE_SOLAR
  @State year: number = 1990
  @State yearText: string = '1990'
  @State month: number = 1
  @State day: number = 1
  @State hour: number = 12
  @State minute: number = 0
  @State isLeap: boolean = false
  @State isSpec: boolean = false
  @State lunarMonths: LunarMonth[] = []
  // 当前年月"真实存在"的日号(公历跳过改历缺口/正确月长; 农历为 1..N)
  @State validDays: number[] = []

  // ── 四柱反推 ──
  @State revIdx: number[] = [0, 2, 0, 0]
  @State revHourUnknown: boolean = false
  @State revStartText: string = '1900'
  @State revEndText: string = '2100'
  @State revResults: ReverseItem[] = []
  @State revSearched: boolean = false

  private readonly tabs: string[] = ['公历', '农历', '四柱反推']
  private readonly pillarNames: string[] = ['年柱', '月柱', '日柱', '时柱']

  aboutToAppear() {
    this.inputMode = this.init.inputMode === MODE_LUNAR ? MODE_LUNAR : MODE_SOLAR
    this.year = this.init.year
    this.yearText = `${this.init.year}`
    this.month = this.init.month
    this.day = this.init.day
    this.hour = this.init.hour
    this.minute = this.init.minute
    this.isLeap = this.init.isLeap
    this.isSpec = this.init.isSpec
    const arr: number[] = this.revIdx.slice()
    arr[1] = BaziCalc.remapByZhi(BaziCalc.validMonths(arr[0]), arr[1])
    arr[3] = BaziCalc.remapByZhi(BaziCalc.validHours(arr[2]), arr[3])
    this.revIdx = arr
    this.refreshLunarMonths()
    this.refreshValidDays()
  }

  refreshLunarMonths() {
    this.lunarMonths = NativeBridge.getLunarMonths(this.year)
  }

  // 重算当前年月有效日列表, 并把已选日吸附到有效日(取不超过它的最大有效日)
  refreshValidDays() {
    let days: number[] = []
    if (this.inputMode === MODE_LUNAR) {
      const n: number = BaziCalc.lunarDaysInMonth(this.year, this.month, this.isLeap, this.isSpec)
      for (let i = 1; i <= n; i++) days.push(i)
    } else {
      days = NativeBridge.getSolarMonthValidDays(this.year, this.month)
      if (days.length === 0) {
        const n: number = BaziCalc.solarDaysInMonth(this.year, this.month)
        for (let i = 1; i <= n; i++) days.push(i)
      }
    }
    this.validDays = days
    if (days.indexOf(this.day) < 0) {
      let pick: number = days.length > 0 ? days[0] : 1
      for (let i = 0; i < days.length; i++) {
        if (days[i] <= this.day) pick = days[i]
      }
      this.day = pick
    }
  }

  // 日轮选项文本 — 农历模式显示"初一/廿九"式名称, 公历显示"1日"
  dayNames(): string[] {
    const r: string[] = []
    for (let i = 0; i < this.validDays.length; i++) {
      const d: number = this.validDays[i]
      r.push(this.inputMode === MODE_LUNAR ? BaziCalc.lunarDayName(d) : `${d}日`)
    }
    return r
  }

  selectedDayIndex(): number {
    const i: number = this.validDays.indexOf(this.day)
    return i < 0 ? 0 : i
  }

  // ── 选项数据 ──
  numRange(from: number, to: number, suffix: string): string[] {
    const r: string[] = []
    for (let i = from; i <= to; i++) r.push(`${i}${suffix}`)
    return r
  }

  lunarMonthNames(): string[] {
    const r: string[] = []
    for (let i = 0; i < this.lunarMonths.length; i++) r.push(this.lunarMonths[i].name)
    return r
  }

  selectedLunarIndex(): number {
    for (let i = 0; i < this.lunarMonths.length; i++) {
      const m: LunarMonth = this.lunarMonths[i]
      if (m.month === this.month && m.isLeap === this.isLeap) return i
    }
    return 0
  }

  onYearInput(v: string) {
    this.yearText = v
    const y: number = BaziCalc.parseYear(v)
    if (!isNaN(y)) {
      this.year = y
      if (this.inputMode === MODE_LUNAR) this.refreshLunarMonths()
      this.refreshValidDays()
    }
  }

  confirmDate() {
    const sel: DateSelection = {
      inputMode: this.inputMode, year: this.year, month: this.month,
      day: this.day, hour: this.hour, minute: this.minute,
      isLeap: this.isLeap, isSpec: this.isSpec
    }
    this.onConfirm(sel)
    this.controller.close()
  }

  build() {
    Column() {
      // 顶部操作条: 取消 / Tab / 确认
      Row() {
        Text('取消').fontSize(AppDimens.FONT_BODY).fontColor(AppColors.TEXT_SECONDARY)
          .onClick(() => this.controller.close())
        Blank()
        Row() {
          ForEach(this.tabs, (t: string, i: number) => {
            Text(t)
              .fontSize(AppDimens.FONT_CAPTION)
              .fontColor(this.inputMode === i ? AppColors.ON_PRIMARY : AppColors.TEXT_SECONDARY)
              .fontWeight(this.inputMode === i ? FontWeight.Bold : FontWeight.Normal)
              .backgroundColor(this.inputMode === i ? AppColors.ACCENT : Color.Transparent)
              .borderRadius(AppDimens.BORDER_RADIUS_SM)
              .padding({ left: 12, right: 12, top: 6, bottom: 6 })
              .onClick(() => {
                this.inputMode = i
                if (i === MODE_LUNAR) this.refreshLunarMonths()
                if (i === MODE_SOLAR || i === MODE_LUNAR) this.refreshValidDays()
              })
          })
        }
        Blank()
        Text('确认').fontSize(AppDimens.FONT_BODY).fontWeight(FontWeight.Bold).fontColor(AppColors.ACCENT)
          .visibility(this.inputMode === MODE_REVERSE ? Visibility.Hidden : Visibility.Visible)
          .onClick(() => this.confirmDate())
      }
      .width('100%').padding({ left: AppDimens.PADDING_LG, right: AppDimens.PADDING_LG, top: AppDimens.PADDING_MD, bottom: AppDimens.PADDING_MD })
      .alignItems(VerticalAlign.Center)
      .border({ width: { bottom: 0.5 }, color: AppColors.DIVIDER })

      if (this.inputMode === MODE_REVERSE) {
        this.reverseBody()
      } else {
        this.dateWheels()
      }
    }
    .width('100%')
    .backgroundColor(AppColors.SURFACE)
    .borderRadius({ topLeft: AppDimens.BORDER_RADIUS_LG, topRight: AppDimens.BORDER_RADIUS_LG })
  }

  // ── 公历/农历: 年(输入) + 月日时分(滚轮) ──
  @Builder
  dateWheels() {
    Column() {
      Row() {
        Text('年份').fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.TEXT_SECONDARY)
        TextInput({ placeholder: '如 1990 / B212 / -211', text: this.yearText })
          .type(InputType.Normal).fontSize(AppDimens.FONT_BODY).layoutWeight(1).height(38)
          .margin({ left: 8 })
          .borderRadius(AppDimens.BORDER_RADIUS_SM).backgroundColor(AppColors.BACKGROUND)
          .onChange((v: string) => this.onYearInput(v))
        Text(BaziCalc.formatYear(this.year)).fontSize(AppDimens.FONT_CAPTION)
          .fontColor(AppColors.PRIMARY).margin({ left: 8 })
      }
      .width('100%').alignItems(VerticalAlign.Center)
      .padding({ left: AppDimens.PADDING_LG, right: AppDimens.PADDING_LG, top: AppDimens.PADDING_MD })

      Row() {
        if (this.inputMode === MODE_LUNAR) {
          TextPicker({ range: this.lunarMonthNames(), selected: this.selectedLunarIndex() })
            .layoutWeight(1)
            .onChange((value: string | string[], index: number | number[]) => {
              const i: number = index as number
              const m: LunarMonth = this.lunarMonths[i]
              this.month = m.month
              this.isLeap = m.isLeap
              this.isSpec = m.isSpec
              this.refreshValidDays()
            })
        } else {
          TextPicker({ range: this.numRange(1, 12, '月'), selected: this.month - 1 })
            .layoutWeight(1)
            .onChange((value: string | string[], index: number | number[]) => {
              this.month = (index as number) + 1
              this.isLeap = false
              this.isSpec = false
              this.refreshValidDays()
            })
        }
        TextPicker({ range: this.dayNames(), selected: this.selectedDayIndex() })
          .layoutWeight(1)
          .onChange((value: string | string[], index: number | number[]) => {
            const i: number = index as number
            if (i >= 0 && i < this.validDays.length) this.day = this.validDays[i]
          })
        TextPicker({ range: this.numRange(0, 23, '时'), selected: this.hour })
          .layoutWeight(1)
          .onChange((value: string | string[], index: number | number[]) => { this.hour = index as number })
        TextPicker({ range: this.numRange(0, 59, '分'), selected: this.minute })
          .layoutWeight(1)
          .onChange((value: string | string[], index: number | number[]) => { this.minute = index as number })
      }
      .width('100%').height(180)
      .padding({ left: AppDimens.PADDING_SM, right: AppDimens.PADDING_SM })
    }
    .width('100%').padding({ bottom: AppDimens.PADDING_LG })
  }

  // ── 四柱反推 ──
  @Builder
  reverseBody() {
    Column() {
      Row() {
        ForEach(this.pillarNames, (name: string, i: number) => {
          Column() {
            Text(name).fontSize(AppDimens.FONT_SMALL).fontColor(AppColors.TEXT_SECONDARY)
            Text(this.revHourUnknown && i === 3 ? '未知' : BaziCalc.jiaZiName(this.revIdx[i]))
              .fontSize(20).fontWeight(FontWeight.Bold)
              .fontColor(this.revHourUnknown && i === 3 ? AppColors.TEXT_SECONDARY : AppColors.PRIMARY)
              .margin({ top: 4 })
            Text('▾').fontSize(AppDimens.FONT_SMALL).fontColor(AppColors.TEXT_SECONDARY)
          }
          .layoutWeight(1).height(70).justifyContent(FlexAlign.Center)
          .backgroundColor(AppColors.BACKGROUND).borderRadius(AppDimens.BORDER_RADIUS_SM)
          .margin({ right: i < 3 ? 6 : 0 })
          .onClick(() => this.pickPillar(i))
        })
      }
      .width('100%').margin({ bottom: AppDimens.PADDING_MD })

      Row() {
        Toggle({ type: ToggleType.Switch, isOn: this.revHourUnknown })
          .onChange((v: boolean) => { this.revHourUnknown = v }).selectedColor(AppColors.ACCENT).width(44).height(24)
        Text('时辰未知').fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.TEXT_SECONDARY).margin({ left: 4 })
      }
      .width('100%').margin({ bottom: AppDimens.PADDING_MD })

      Row() {
        Text('搜索范围').fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.TEXT_SECONDARY)
        TextInput({ text: this.revStartText })
          .fontSize(AppDimens.FONT_CAPTION).width(80).height(36).margin({ left: 6 })
          .borderRadius(AppDimens.BORDER_RADIUS_SM).backgroundColor(AppColors.BACKGROUND)
          .onChange((v: string) => { this.revStartText = v })
        Text('—').fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.TEXT_SECONDARY).margin({ left: 6, right: 6 })
        TextInput({ text: this.revEndText })
          .fontSize(AppDimens.FONT_CAPTION).width(80).height(36)
          .borderRadius(AppDimens.BORDER_RADIUS_SM).backgroundColor(AppColors.BACKGROUND)
          .onChange((v: string) => { this.revEndText = v })
        Blank()
        Button('反推').height(36).fontSize(AppDimens.FONT_CAPTION)
          .fontColor(AppColors.PRIMARY_DARK).backgroundColor(AppColors.ACCENT)
          .borderRadius(AppDimens.BORDER_RADIUS_SM)
          .onClick(() => this.doReverse())
      }
      .width('100%').alignItems(VerticalAlign.Center).margin({ bottom: AppDimens.PADDING_SM })

      if (this.revSearched) {
        Text(`匹配结果（${this.revResults.length} 条${this.revResults.length >= 60 ? '，已截断' : ''}），点击选用`)
          .fontSize(AppDimens.FONT_SMALL).fontColor(AppColors.TEXT_SECONDARY)
          .width('100%').margin({ bottom: AppDimens.PADDING_XS })
        List() {
          ForEach(this.revResults, (item: ReverseItem) => {
            ListItem() {
              Row() {
                Text(`${BaziCalc.formatYear(item.year)}${item.month}月${item.day}日`
                  + (item.hour >= 0 ? ` ${item.hour}时` : ''))
                  .fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.ON_SURFACE).layoutWeight(1)
                Text('选用 ›').fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.ACCENT)
              }
              .width('100%').padding({ top: 8, bottom: 8 })
              .border({ width: { bottom: 0.5 }, color: AppColors.DIVIDER })
              .onClick(() => this.applyReverse(item))
            }
          })
        }
        .width('100%').height(160).scrollBar(BarState.Off)
      }
    }
    .width('100%')
    .padding({ left: AppDimens.PADDING_LG, right: AppDimens.PADDING_LG, top: AppDimens.PADDING_MD, bottom: AppDimens.PADDING_LG })
  }

  pickPillar(i: number) {
    if (i === 3 && this.revHourUnknown) return
    let list: number[] = []
    if (i === 1) list = BaziCalc.validMonths(this.revIdx[0])
    else if (i === 3) list = BaziCalc.validHours(this.revIdx[2])
    else { for (let n = 0; n < 60; n++) list.push(n) }
    const names: string[] = []
    for (let k = 0; k < list.length; k++) names.push(BaziCalc.jiaZiName(list[k]))
    let sel: number = list.indexOf(((this.revIdx[i] % 60) + 60) % 60)
    if (sel < 0) sel = 0
    this.getUIContext().showTextPickerDialog({
      range: names,
      selected: sel,
      onAccept: (v: TextPickerResult) => {
        const arr: number[] = this.revIdx.slice()
        arr[i] = list[v.index as number]
        if (i === 0) arr[1] = BaziCalc.remapByZhi(BaziCalc.validMonths(arr[0]), arr[1])
        if (i === 2) arr[3] = BaziCalc.remapByZhi(BaziCalc.validHours(arr[2]), arr[3])
        this.revIdx = arr
      }
    })
  }

  doReverse() {
    const sy: number = BaziCalc.parseYear(this.revStartText)
    const ey: number = BaziCalc.parseYear(this.revEndText)
    if (isNaN(sy) || isNaN(ey) || sy > ey) {
      this.revResults = []
      this.revSearched = true
      return
    }
    const g: number[] = []
    const z: number[] = []
    for (let i = 0; i < 4; i++) {
      const idx: number = ((this.revIdx[i] % 60) + 60) % 60
      g.push(idx % 10)
      z.push(idx % 12)
    }
    const hg: number = this.revHourUnknown ? -1 : g[3]
    const hz: number = this.revHourUnknown ? -1 : z[3]
    try {
      this.revResults = NativeBridge.baziReverse(g[0], z[0], g[1], z[1], g[2], z[2], hg, hz, sy, ey)
    } catch (e) {
      this.revResults = []
    }
    this.revSearched = true
  }

  applyReverse(item: ReverseItem) {
    const sel: DateSelection = {
      inputMode: MODE_SOLAR, year: item.year, month: item.month, day: item.day,
      hour: item.hour >= 0 ? item.hour : 12, minute: 0, isLeap: false, isSpec: false
    }
    this.onConfirm(sel)
    this.controller.close()
  }
}
