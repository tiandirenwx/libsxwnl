import { AppColors, AppDimens } from '../common/Constants'
import { Cities, Region, ObserverCity } from '../common/Cities'

// 出生地选择弹窗 — 左省右市两栏 List 联动
//
// 用法:
//   private cityDialog = new CustomDialogController({
//     builder: CityPickerDialog({
//       initRegionIdx: this.regionIdx,
//       initCityIdx: this.cityIdx,
//       onConfirm: (rIdx, cIdx, city) => { ... }
//     }),
//     alignment: DialogAlignment.Bottom,
//     customStyle: true,
//     autoCancel: true
//   })
@CustomDialog
export struct CityPickerDialog {
  controller: CustomDialogController
  initRegionIdx: number = 0
  initCityIdx: number = 0
  onConfirm: (regionIdx: number, cityIdx: number, city: ObserverCity) => void =
    (regionIdx: number, cityIdx: number, city: ObserverCity): void => {}

  @State tmpRegionIdx: number = 0     // 弹窗内左栏选中(临时)
  @State tmpCityIdx: number = 0       // 弹窗内右栏选中(临时, 用户切换省时跟随左栏首条)

  aboutToAppear() {
    this.tmpRegionIdx = this.initRegionIdx
    this.tmpCityIdx = this.initCityIdx
  }

  build() {
    Column() {
      // 标题栏
      Row() {
        Text('取消')
          .fontSize(AppDimens.FONT_BODY).fontColor(AppColors.TEXT_SECONDARY)
          .padding({ left: 16, right: 16, top: 12, bottom: 12 })
          .onClick((): void => { this.controller.close() })
        Text('选择出生地')
          .fontSize(AppDimens.FONT_SUBTITLE).fontWeight(FontWeight.Bold)
          .fontColor(AppColors.ON_SURFACE).layoutWeight(1).textAlign(TextAlign.Center)
        Text('确定')
          .fontSize(AppDimens.FONT_BODY).fontColor(AppColors.PRIMARY)
          .padding({ left: 16, right: 16, top: 12, bottom: 12 })
          .onClick((): void => {
            const region: Region = Cities.REGIONS[this.tmpRegionIdx]
            const city: ObserverCity = region.cities[this.tmpCityIdx]
            this.onConfirm(this.tmpRegionIdx, this.tmpCityIdx, city)
            this.controller.close()
          })
      }
      .width('100%')
      .border({ width: { bottom: 0.5 }, color: AppColors.DIVIDER })

      // 当前选中预览
      Row() {
        Text(this.currentSummary())
          .fontSize(AppDimens.FONT_SMALL).fontColor(AppColors.TEXT_SECONDARY)
          .layoutWeight(1).maxLines(1).textOverflow({ overflow: TextOverflow.Ellipsis })
      }
      .width('100%').padding({ left: 16, right: 16, top: 8, bottom: 8 })

      // 左省 + 右市 两栏 List
      Row() {
        // 左栏: 省/直辖市
        List({ space: 0 }) {
          ForEach(Cities.REGIONS, (r: Region, idx: number): void => {
            ListItem() {
              Row() {
                Text(r.name)
                  .fontSize(AppDimens.FONT_CAPTION)
                  .fontColor(this.tmpRegionIdx === idx ? AppColors.PRIMARY : AppColors.ON_SURFACE)
                  .fontWeight(this.tmpRegionIdx === idx ? FontWeight.Bold : FontWeight.Normal)
                  .layoutWeight(1)
              }
              .width('100%').height(44).padding({ left: 14, right: 8 })
              .backgroundColor(this.tmpRegionIdx === idx ? AppColors.BACKGROUND : AppColors.SURFACE)
              .alignItems(VerticalAlign.Center)
              .onClick((): void => {
                if (this.tmpRegionIdx !== idx) {
                  this.tmpRegionIdx = idx
                  this.tmpCityIdx = 0
                }
              })
            }
          }, (r: Region, idx: number): string => `r_${idx}_${r.name}`)
        }
        .layoutWeight(1).height('100%')
        .backgroundColor(AppColors.SURFACE)
        .scrollBar(BarState.Auto)

        // 右栏: 区/县 (随左栏切换)
        List({ space: 0 }) {
          ForEach(Cities.REGIONS[this.tmpRegionIdx].cities,
            (c: ObserverCity, idx: number): void => {
              ListItem() {
                Row() {
                  Text(c.name)
                    .fontSize(AppDimens.FONT_CAPTION)
                    .fontColor(this.tmpCityIdx === idx ? AppColors.PRIMARY : AppColors.ON_SURFACE)
                    .fontWeight(this.tmpCityIdx === idx ? FontWeight.Bold : FontWeight.Normal)
                    .layoutWeight(1)
                }
                .width('100%').height(44).padding({ left: 14, right: 8 })
                .backgroundColor(this.tmpCityIdx === idx ? AppColors.BACKGROUND : AppColors.SURFACE)
                .alignItems(VerticalAlign.Center)
                .onClick((): void => {
                  this.tmpCityIdx = idx
                })
              }
            },
            (c: ObserverCity, idx: number): string => `c_${this.tmpRegionIdx}_${idx}_${c.name}`)
        }
        .layoutWeight(2).height('100%')
        .backgroundColor(AppColors.BACKGROUND)
        .scrollBar(BarState.Auto)
      }
      .width('100%').height(360)
    }
    .width('100%')
    .backgroundColor(AppColors.SURFACE)
    .borderRadius({ topLeft: 16, topRight: 16 })
  }

  currentSummary(): string {
    const r: Region = Cities.REGIONS[this.tmpRegionIdx]
    const c: ObserverCity = r.cities[this.tmpCityIdx]
    const ew: string = c.longitude >= 0 ? 'E' : 'W'
    const ns: string = c.latitude >= 0 ? 'N' : 'S'
    const lon: string = Math.abs(c.longitude).toFixed(2)
    const lat: string = Math.abs(c.latitude).toFixed(2)
    return `${r.name} · ${c.name} · ${lon}°${ew}, ${lat}°${ns}`
  }
}
