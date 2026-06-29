import { AppColors, AppDimens } from '../common/Constants'
import { CalendarPage } from './CalendarPage'
import { YearCalendarPage } from './YearCalendarPage'
import { BaziPage } from './BaziPage'
import { EclipsePage } from './EclipsePage'

@Entry
@Component
struct Index {
  @State currentIndex: number = 0

  @Builder
  tabBuilder(title: string, icon: string, index: number) {
    Column() {
      Text(icon)
        .fontSize(22)
        .fontColor(this.currentIndex === index ? AppColors.ACCENT : AppColors.TEXT_SECONDARY)
        .margin({ bottom: 2 })
      Text(title)
        .fontSize(AppDimens.FONT_SMALL)
        .fontColor(this.currentIndex === index ? AppColors.ACCENT : AppColors.TEXT_SECONDARY)
        .fontWeight(this.currentIndex === index ? FontWeight.Medium : FontWeight.Normal)
    }
    .width('100%')
    .height(AppDimens.TAB_BAR_HEIGHT)
    .justifyContent(FlexAlign.Center)
    .alignItems(HorizontalAlign.Center)
  }

  build() {
    Column() {
      Tabs({ barPosition: BarPosition.End, index: this.currentIndex }) {
        TabContent() {
          CalendarPage()
        }
        .tabBar(this.tabBuilder('月历', '📅', 0))

        TabContent() {
          YearCalendarPage()
        }
        .tabBar(this.tabBuilder('年历', '📆', 1))

        TabContent() {
          BaziPage()
        }
        .tabBar(this.tabBuilder('八字', '☯', 2))

        TabContent() {
          EclipsePage()
        }
        .tabBar(this.tabBuilder('日月食', '🌑', 3))
      }
      .barMode(BarMode.Fixed)
      .barWidth('100%')
      .barHeight(AppDimens.TAB_BAR_HEIGHT)
      .animationDuration(200)
      .onChange((index: number) => {
        this.currentIndex = index
      })
    }
    .width('100%')
    .height('100%')
    .backgroundColor(AppColors.BACKGROUND)
  }
}
