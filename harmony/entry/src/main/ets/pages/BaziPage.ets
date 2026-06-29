import { AppColors, AppDimens } from '../common/Constants'
import { NativeBridge, BaziResult, LunarMonth } from '../services/NativeBridge'
import {
  BaziForm, BaziCalc, DateSelection, LiFaOption, LIFA_OPTIONS, LIFA_DING_QI, MODE_SOLAR
} from '../model/BaziInput'
import { Cities, Region, ObserverCity } from '../common/Cities'
import { DateTimePicker } from '../components/DateTimePicker'
import { CityPickerDialog } from '../components/CityPickerDialog'
import { BaziResultView, BaziResultArg } from '../components/BaziResultView'

@Component
export struct BaziPage {
  @State name: string = ''
  @State gender: number = 0
  @State lifa: number = LIFA_DING_QI
  @State astEnabled: boolean = false
  // 出生地状态: 默认选择「北京市 / 天安门」
  // 手动开关关闭=省份+区县级联下拉(自动填经纬度); 开启=直接手动输入经纬度
  @State manualGeo: boolean = false
  @State regionIdx: number = 0
  @State cityIdx: number = 0
  @State longitude: number = Cities.DEFAULT.longitude
  @State latitude: number = Cities.DEFAULT.latitude
  @State sel: DateSelection = {
    inputMode: MODE_SOLAR, year: 1990, month: 1, day: 1,
    hour: 12, minute: 0, isLeap: false, isSpec: false
  }
  @State lunarMonths: LunarMonth[] = []
  @State picked: boolean = false

  private pageStack: NavPathStack = new NavPathStack()

  private dateDialog: CustomDialogController = new CustomDialogController({
    builder: DateTimePicker({
      init: this.sel,
      onConfirm: (sel: DateSelection): void => this.onPickConfirm(sel)
    }),
    alignment: DialogAlignment.Bottom,
    customStyle: true,
    autoCancel: true
  })

  private cityDialog: CustomDialogController = new CustomDialogController({
    builder: CityPickerDialog({
      initRegionIdx: this.regionIdx,
      initCityIdx: this.cityIdx,
      onConfirm: (rIdx: number, cIdx: number, city: ObserverCity): void =>
        this.onCityDialogConfirm(rIdx, cIdx, city)
    }),
    alignment: DialogAlignment.Bottom,
    customStyle: true,
    autoCancel: true
  })

  aboutToAppear() {
    try {
      this.getUIContext().getFont().registerFont({ familyName: 'WenYue', familySrc: $rawfile('WenYue.otf') })
    } catch (e) {
      console.error('registerFont failed: ' + JSON.stringify(e))
    }
    this.lunarMonths = NativeBridge.getLunarMonths(this.sel.year)
  }

  onPickConfirm(sel: DateSelection) {
    this.sel = sel
    this.lunarMonths = NativeBridge.getLunarMonths(sel.year)
    this.picked = true
  }

  recordText(): string {
    if (!this.picked) return '请选择出生时间'
    return BaziCalc.formatRecord(this.sel, this.lunarMonths)
  }

  paiPan() {
    // 八字时间一律按北京时间处理, 经纬度仅用于真太阳时方程修正
    const form: BaziForm = {
      name: this.name, gender: this.gender, lifa: this.lifa,
      astEnabled: this.astEnabled,
      longitude: this.longitude, latitude: this.latitude,
      sel: this.sel
    }
    try {
      const result: BaziResult = NativeBridge.calcBazi(BaziCalc.toParams(form))
      let lifaLabel: string = '定气'
      for (let i = 0; i < LIFA_OPTIONS.length; i++) {
        if (LIFA_OPTIONS[i].value === this.lifa) lifaLabel = LIFA_OPTIONS[i].label
      }
      const arg: BaziResultArg = {
        result: result, birthYear: this.sel.year, astEnabled: this.astEnabled,
        longitude: this.longitude, latitude: this.latitude, lifaLabel: lifaLabel
      }
      this.pageStack.pushPathByName('result', arg)
    } catch (e) {
      console.error('Bazi calculation failed: ' + JSON.stringify(e))
      try {
        this.getUIContext().getPromptAction().showToast({ message: '排盘失败' })
      } catch (err) {
        console.error('showToast failed: ' + JSON.stringify(err))
      }
    }
  }

  build() {
    Navigation(this.pageStack) {
      Scroll() {
        Column() {
          this.headerBanner()
          this.inputCard()
        }
        .width('100%').padding({ bottom: 32 })
      }
      .width('100%').height('100%')
      .backgroundColor(AppColors.BACKGROUND)
      .scrollBar(BarState.Off)
    }
    .hideTitleBar(true)
    .mode(NavigationMode.Stack)
    .navDestination(this.resultPage)
  }

  // 结果页(新开页面)
  @Builder
  resultPage(name: string, param: object) {
    NavDestination() {
      BaziResultView({ arg: param as BaziResultArg })
    }
    .title((param as BaziResultArg).result.gender !== undefined
      && (param as BaziResultArg).result.gender.indexOf('女') >= 0 ? '排盘 · 坤造' : '排盘 · 乾造')
    .backgroundColor(AppColors.BACKGROUND)
  }

  @Builder
  headerBanner() {
    Column() {
      Text('☯ 八字排盘')
        .fontSize(AppDimens.FONT_TITLE).fontWeight(FontWeight.Bold).fontColor(AppColors.ON_PRIMARY)
    }
    .width('100%')
    .padding({ top: AppDimens.PADDING_XL, bottom: AppDimens.PADDING_LG })
    .alignItems(HorizontalAlign.Center)
    .linearGradient({
      direction: GradientDirection.RightBottom,
      colors: [[AppColors.GRADIENT_START, 0.0], [AppColors.GRADIENT_END, 1.0]]
    })
  }

  @Builder
  inputCard() {
    Column() {
      // 姓名 + 性别
      Row() {
        TextInput({ placeholder: '姓名(选填)', text: this.name })
          .fontSize(AppDimens.FONT_BODY).layoutWeight(1).height(44)
          .borderRadius(AppDimens.BORDER_RADIUS_SM).backgroundColor(AppColors.BACKGROUND)
          .onChange((v: string) => { this.name = v })
        Row() {
          this.genderChip('男', 0)
          this.genderChip('女', 1)
        }
        .margin({ left: 10 })
        .border({ width: 1, color: AppColors.DIVIDER, radius: AppDimens.BORDER_RADIUS_SM })
        .clip(true)
      }
      .width('100%').margin({ bottom: AppDimens.PADDING_MD })

      // 出生时间(一条记录, 点击弹出选择器)
      Row() {
        Text(this.recordText())
          .fontSize(AppDimens.FONT_SUBTITLE)
          .fontColor(this.picked ? AppColors.ON_SURFACE : AppColors.TEXT_SECONDARY)
          .layoutWeight(1)
        Text('▾').fontSize(AppDimens.FONT_BODY).fontColor(AppColors.TEXT_SECONDARY)
      }
      .width('100%').height(48).padding({ left: 14, right: 14 })
      .backgroundColor(AppColors.BACKGROUND).borderRadius(AppDimens.BORDER_RADIUS_SM)
      .alignItems(VerticalAlign.Center)
      .margin({ bottom: AppDimens.PADDING_MD })
      .onClick(() => this.dateDialog.open())

      // 出生地: 城市选择(默认) 或 手动输入经纬度(开关切换)
      this.birthPlaceSection()

      // 历法选择: 定气 / 平气定夏至 / 平气定冬至
      Row() {
        Text('历法').fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.TEXT_SECONDARY).width(48)
        Row() {
          ForEach(LIFA_OPTIONS, (opt: LiFaOption, i: number) => {
            Text(opt.label)
              .fontSize(AppDimens.FONT_CAPTION)
              .fontColor(this.lifa === opt.value ? AppColors.ON_PRIMARY : AppColors.TEXT_SECONDARY)
              .fontWeight(this.lifa === opt.value ? FontWeight.Bold : FontWeight.Normal)
              .backgroundColor(this.lifa === opt.value ? AppColors.PRIMARY : AppColors.BACKGROUND)
              .textAlign(TextAlign.Center).layoutWeight(1).height(34)
              .borderRadius(AppDimens.BORDER_RADIUS_SM)
              .margin({ right: i < LIFA_OPTIONS.length - 1 ? 6 : 0 })
              .onClick(() => { this.lifa = opt.value })
          })
        }
        .layoutWeight(1)
      }
      .width('100%').alignItems(VerticalAlign.Center).margin({ bottom: AppDimens.PADDING_MD })

      // 真太阳时
      Row() {
        Checkbox().select(this.astEnabled).selectedColor(AppColors.PRIMARY).width(20).height(20)
          .onChange((v: boolean) => { this.astEnabled = v })
        Text('使用真太阳时').fontSize(AppDimens.FONT_BODY).fontColor(AppColors.ON_SURFACE).margin({ left: 6 })
      }
      .width('100%').margin({ bottom: AppDimens.PADDING_LG })

      Button('排  盘')
        .width('100%').height(48).fontSize(AppDimens.FONT_SUBTITLE).fontWeight(FontWeight.Bold)
        .fontColor(AppColors.PRIMARY_DARK).backgroundColor(AppColors.ACCENT)
        .borderRadius(AppDimens.BORDER_RADIUS_LG)
        .onClick(() => this.paiPan())
    }
    .width('100%').padding(AppDimens.PADDING_LG)
    .margin({ left: AppDimens.PADDING_MD, right: AppDimens.PADDING_MD, top: AppDimens.PADDING_MD })
    .backgroundColor(AppColors.SURFACE).borderRadius(AppDimens.BORDER_RADIUS_LG)
    .shadow({ radius: 8, color: AppColors.CARD_SHADOW_COLOR, offsetY: 2 })
  }

  // ── 出生地选择 ──────────────────────────────────────────
  //   关闭手动开关: 整行点击 → 弹出底部两栏(省/区县)联动选择器 CityPickerDialog
  //   开启手动开关: 直接输入经纬度
  //   时间一律按北京时间处理, 经纬度仅用于真太阳时方程修正
  @Builder
  birthPlaceSection() {
    Column() {
      Row() {
        Text('出生地').fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.TEXT_SECONDARY).width(48)
        Text(this.birthPlaceSummary())
          .fontSize(AppDimens.FONT_SMALL).fontColor(AppColors.TEXT_SECONDARY).layoutWeight(1)
          .maxLines(1).textOverflow({ overflow: TextOverflow.Ellipsis })
        Text('手动')
          .fontSize(AppDimens.FONT_SMALL).fontColor(AppColors.TEXT_SECONDARY).margin({ right: 4 })
        Toggle({ type: ToggleType.Switch, isOn: this.manualGeo })
          .selectedColor(AppColors.PRIMARY).width(36).height(20)
          .onChange((v: boolean): void => { this.manualGeo = v })
      }
      .width('100%').alignItems(VerticalAlign.Center).margin({ bottom: 6 })

      if (this.manualGeo) {
        Row() {
          Text('经度').fontSize(AppDimens.FONT_SMALL).fontColor(AppColors.TEXT_SECONDARY)
          TextInput({ text: `${this.longitude}` })
            .type(InputType.NUMBER_DECIMAL).fontSize(AppDimens.FONT_CAPTION)
            .layoutWeight(1).height(38).margin({ left: 4, right: 8 })
            .borderRadius(AppDimens.BORDER_RADIUS_SM).backgroundColor(AppColors.BACKGROUND)
            .onChange((v: string): void => { const n: number = Number(v); if (!isNaN(n)) this.longitude = n })
          Text('纬度').fontSize(AppDimens.FONT_SMALL).fontColor(AppColors.TEXT_SECONDARY)
          TextInput({ text: `${this.latitude}` })
            .type(InputType.NUMBER_DECIMAL).fontSize(AppDimens.FONT_CAPTION)
            .layoutWeight(1).height(38).margin({ left: 4 })
            .borderRadius(AppDimens.BORDER_RADIUS_SM).backgroundColor(AppColors.BACKGROUND)
            .onChange((v: string): void => { const n: number = Number(v); if (!isNaN(n)) this.latitude = n })
        }
        .width('100%').alignItems(VerticalAlign.Center)
      } else {
        // 触发器: 点击整行打开两栏弹窗
        Row() {
          Text(this.regionCityLabel())
            .fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.ON_SURFACE)
            .layoutWeight(1).maxLines(1).textOverflow({ overflow: TextOverflow.Ellipsis })
          Text('▾').fontSize(AppDimens.FONT_SMALL).fontColor(AppColors.TEXT_SECONDARY)
        }
        .width('100%').height(38).padding({ left: 14, right: 14 })
        .backgroundColor(AppColors.BACKGROUND).borderRadius(AppDimens.BORDER_RADIUS_SM)
        .alignItems(VerticalAlign.Center)
        .onClick((): void => this.openCityDialog())
      }
    }
    .width('100%').margin({ bottom: AppDimens.PADDING_MD })
  }

  // 概要文字 "北京市 · 天安门 · 116.38°E, 39.90°N" 或手动模式 "116.38°E, 39.90°N"
  birthPlaceSummary(): string {
    const ew: string = this.longitude >= 0 ? 'E' : 'W'
    const ns: string = this.latitude  >= 0 ? 'N' : 'S'
    const lon: string = Math.abs(this.longitude).toFixed(2)
    const lat: string = Math.abs(this.latitude).toFixed(2)
    if (this.manualGeo) {
      return `${lon}°${ew}, ${lat}°${ns}`
    }
    const r: Region = Cities.REGIONS[this.regionIdx]
    const c: ObserverCity = r.cities[this.cityIdx]
    return `${r.name} · ${c.name} · ${lon}°${ew}, ${lat}°${ns}`
  }

  // 触发行的紧凑标签
  regionCityLabel(): string {
    const r: Region = Cities.REGIONS[this.regionIdx]
    const c: ObserverCity = r.cities[this.cityIdx]
    return `${r.name}  ▸  ${c.name}`
  }

  openCityDialog(): void {
    // 重建 dialog 控制器, 把当前选中传进去作为初始值
    this.cityDialog = new CustomDialogController({
      builder: CityPickerDialog({
        initRegionIdx: this.regionIdx,
        initCityIdx: this.cityIdx,
        onConfirm: (rIdx: number, cIdx: number, city: ObserverCity): void =>
          this.onCityDialogConfirm(rIdx, cIdx, city)
      }),
      alignment: DialogAlignment.Bottom,
      customStyle: true,
      autoCancel: true
    })
    this.cityDialog.open()
  }

  onCityDialogConfirm(rIdx: number, cIdx: number, city: ObserverCity): void {
    this.regionIdx = rIdx
    this.cityIdx = cIdx
    this.longitude = city.longitude
    this.latitude  = city.latitude
  }

  @Builder
  genderChip(label: string, value: number) {
    Text(label)
      .fontSize(AppDimens.FONT_BODY).fontWeight(this.gender === value ? FontWeight.Bold : FontWeight.Normal)
      .fontColor(this.gender === value ? AppColors.ON_PRIMARY : AppColors.TEXT_SECONDARY)
      .backgroundColor(this.gender === value ? AppColors.PRIMARY : AppColors.SURFACE)
      .textAlign(TextAlign.Center).width(48).height(44)
      .onClick(() => { this.gender = value })
  }
}
