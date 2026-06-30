import { AppColors, AppDimens } from '../common/Constants'
import {
  AlmanacQuote, AlmanacTopic, DayAlmanac, DayInfo, EventAdvice, LuckyHour, ShenSha
} from '../services/NativeBridge'

const ZHI_NAMES: string[] = ['子','丑','寅','卯','辰','巳','午','未','申','酉','戌','亥']
// 各地支对应的现代时段, 标注以 23:00 起算; 子为 23-01.
const ZHI_HOURS: string[] = [
  '23:00-01:00', '01:00-03:00', '03:00-05:00', '05:00-07:00',
  '07:00-09:00', '09:00-11:00', '11:00-13:00', '13:00-15:00',
  '15:00-17:00', '17:00-19:00', '19:00-21:00', '21:00-23:00',
]

// ════════════════════════════════════════════════════════════════
//  AlmanacComponents — 老黄历展示组件
//
//  · AlmanacSheet:    底部抽屉式完整老黄历
//  · AlmanacInlineRow: popup 内可嵌入的一行摘要 (一眼速览)
//  · AlmanacBadge:    小色块, 用于"吉/凶"等关键标记
//
//  数据来源: NativeBridge.getAlmanac (-> capi -> sxwnl::almanac)
//  设计哲学: 大色块强信号 (吉绿 / 凶红 / 黄道金), 字号紧凑
// ════════════════════════════════════════════════════════════════

class AlmanacPalette {
  static readonly LUCKY: string  = '#3DA76A'   // 吉 - 绿
  static readonly UNLUCKY: string = '#C75450'  // 凶 - 红
  static readonly HUANG_DAO: string = '#D4A015' // 黄道 - 金
  static readonly HEI_DAO: string = '#666666'   // 黑道 - 灰
  static readonly LABEL: string = '#9097B5'
  static readonly VALUE: string = '#FFFFFF'
  static readonly SUB:   string = '#DCE0F2'
  static readonly DIVIDER: string = '#3B3F73'
  static readonly BG_DEEP: string = '#14163A'
  static readonly BG_SECTION: string = '#1B1D44'
}

// ── 多行摘要 (供 CalendarPage 的 popup 直接调用) ─────────────
//   第 1 行: 二十八宿 + 十二神 + 冲煞
//   第 2 行: 宜 X / 忌 Y (神煞投票后聚合)
@Component
export struct AlmanacInlineRow {
  almanac: DayAlmanac | null = null

  build() {
    if (this.almanac !== null) {
      Column({ space: 2 }) {
        Row({ space: 6 }) {
          Text('📜').fontSize(11)
          Text(`${this.almanac.xiu}${this.almanac.xiuZheng}${this.almanac.xiuAnimal}`)
            .fontSize(AppDimens.FONT_CAPTION)
            .fontColor(this.almanac.xiuLuck === '吉'
              ? AlmanacPalette.LUCKY : AlmanacPalette.UNLUCKY)
            .fontWeight(FontWeight.Medium)
          Text('·').fontSize(AppDimens.FONT_SMALL).fontColor(AlmanacPalette.SUB)
          Text(`${this.almanac.twelveGod}${this.almanac.huangHei}`)
            .fontSize(AppDimens.FONT_CAPTION)
            .fontColor(this.almanac.isHuangDao
              ? AlmanacPalette.HUANG_DAO : AlmanacPalette.HEI_DAO)
          Text('·').fontSize(AppDimens.FONT_SMALL).fontColor(AlmanacPalette.SUB)
          Text(`冲${this.almanac.chongShengXiao}煞${this.almanac.sha}`)
            .fontSize(AppDimens.FONT_CAPTION)
            .fontColor(AlmanacPalette.SUB)
        }
        .width('100%')
        .justifyContent(FlexAlign.Center)
        .alignItems(VerticalAlign.Center)

        if ((this.almanac.yi && this.almanac.yi.length > 0)
          || (this.almanac.ji && this.almanac.ji.length > 0)) {
          Row({ space: 6 }) {
            if (this.almanac.yi.length > 0) {
              Text('宜').fontSize(AppDimens.FONT_CAPTION)
                .fontColor('#FFFFFF').fontWeight(FontWeight.Medium)
                .backgroundColor(AlmanacPalette.LUCKY)
                .padding({ left: 5, right: 5, top: 1, bottom: 1 })
                .borderRadius(4)
              Text(this.almanac.yi.slice(0, 3).join(' '))
                .fontSize(AppDimens.FONT_CAPTION)
                .fontColor(AlmanacPalette.SUB)
            }
            if (this.almanac.ji.length > 0) {
              Text('忌').fontSize(AppDimens.FONT_CAPTION)
                .fontColor('#FFFFFF').fontWeight(FontWeight.Medium)
                .backgroundColor(AlmanacPalette.UNLUCKY)
                .padding({ left: 5, right: 5, top: 1, bottom: 1 })
                .borderRadius(4)
                .margin({ left: this.almanac.yi.length > 0 ? 6 : 0 })
              Text(this.almanac.ji.slice(0, 3).join(' '))
                .fontSize(AppDimens.FONT_CAPTION)
                .fontColor(AlmanacPalette.SUB)
            }
          }
          .width('100%')
          .justifyContent(FlexAlign.Center)
          .alignItems(VerticalAlign.Center)
        }
      }
      .width('100%')
      .margin({ top: 1, bottom: 1 })
    }
  }
}

// ── 底部抽屉: 完整老黄历卡片 ────────────────────────────────
//   状态 (day/almanac) 全部由父级 (CalendarPage) 注入,
//   显隐由父级 bindSheet 控制 - 本组件保持纯展示, 不持有任何业务逻辑.
@Component
export struct AlmanacSheet {
  day: DayInfo | null = null
  almanac: DayAlmanac | null = null
  onOpenTopics?: () => void

  build() {
    Scroll() {
      Column() {
        if (this.day !== null && this.almanac !== null) {
          this.headerSection()
          this.yiJiSection()
          this.shenShaSection()
          this.luckyHoursSection()
          this.eventsSection()
          this.xiuSection()
          this.directionSection()
          this.pengZuSection()
          this.quotesSection()
          this.notesSection()
          this.topicsEntrySection()
        } else {
          Text('数据加载中...')
            .fontSize(AppDimens.FONT_BODY)
            .fontColor(AlmanacPalette.SUB)
            .padding(40)
        }
      }
      .width('100%')
      .padding({ left: AppDimens.PADDING_LG, right: AppDimens.PADDING_LG,
                 top: AppDimens.PADDING_MD, bottom: AppDimens.PADDING_XL })
    }
    .width('100%')
    .height('100%')
    .scrollBar(BarState.Auto)
    .backgroundColor(AlmanacPalette.BG_DEEP)
  }

  // ── 顶部: 日期 + 干支 + 12 神(黄道/黑道大徽章) ───────────
  @Builder
  headerSection() {
    Row() {
      Column() {
        Text(`${this.day!.solarYear}年${this.day!.solarMonth}月${this.day!.solarDay}日`)
          .fontSize(AppDimens.FONT_SUBTITLE)
          .fontWeight(FontWeight.Bold)
          .fontColor(AlmanacPalette.VALUE)
        Text(`${this.day!.lunarMonthName}${this.day!.lunarDayName} · ${this.day!.yearGZ}年 ${this.day!.monthGZ}月 ${this.day!.dayGZ}日`)
          .fontSize(AppDimens.FONT_CAPTION)
          .fontColor(AlmanacPalette.SUB)
          .margin({ top: 4 })
      }
      .alignItems(HorizontalAlign.Start)
      .layoutWeight(1)

      Column() {
        Text(this.almanac!.twelveGod)
          .fontSize(AppDimens.FONT_SUBTITLE)
          .fontWeight(FontWeight.Bold)
          .fontColor('#FFFFFF')
        Text(this.almanac!.huangHei)
          .fontSize(AppDimens.FONT_SMALL)
          .fontColor('#FFFFFF')
          .margin({ top: 2 })
      }
      .padding({ left: 12, right: 12, top: 6, bottom: 6 })
      .borderRadius(AppDimens.BORDER_RADIUS_MD)
      .backgroundColor(this.almanac!.isHuangDao
        ? AlmanacPalette.HUANG_DAO : AlmanacPalette.HEI_DAO)
      .alignItems(HorizontalAlign.Center)
    }
    .width('100%')
    .alignItems(VerticalAlign.Center)
    .margin({ bottom: AppDimens.PADDING_LG })
  }

  // ── 二十八宿 · 冲煞 ───────────────────────────────────────
  @Builder
  xiuSection() {
    Column() {
      this.sectionTitle('星宿 · 冲煞')
      Row({ space: 12 }) {
        Column() {
          Text(`${this.almanac!.xiu}${this.almanac!.xiuZheng}${this.almanac!.xiuAnimal}`)
            .fontSize(AppDimens.FONT_TITLE)
            .fontWeight(FontWeight.Bold)
            .fontColor(this.almanac!.xiuLuck === '吉'
              ? AlmanacPalette.LUCKY : AlmanacPalette.UNLUCKY)
          Text(`${this.almanac!.xiuGong}宫 · ${this.almanac!.xiuLuck}`)
            .fontSize(AppDimens.FONT_SMALL)
            .fontColor(AlmanacPalette.SUB)
            .margin({ top: 4 })
        }
        .layoutWeight(1)
        .alignItems(HorizontalAlign.Start)

        Divider().vertical(true).height(40).color(AlmanacPalette.DIVIDER)

        Column() {
          Text(`冲${this.almanac!.chongShengXiao}(${this.almanac!.chongGanZhi})`)
            .fontSize(AppDimens.FONT_BODY)
            .fontWeight(FontWeight.Medium)
            .fontColor(AlmanacPalette.VALUE)
          Text(`煞${this.almanac!.sha}`)
            .fontSize(AppDimens.FONT_SMALL)
            .fontColor(AppColors.WEEKEND_COLOR + 'CC')
            .margin({ top: 4 })
        }
        .layoutWeight(1)
        .alignItems(HorizontalAlign.Start)
      }
      .width('100%')
      .alignItems(VerticalAlign.Center)
    }
    .width('100%')
    .padding(AppDimens.PADDING_MD)
    .backgroundColor(AlmanacPalette.BG_SECTION)
    .borderRadius(AppDimens.BORDER_RADIUS_MD)
    .margin({ bottom: AppDimens.PADDING_MD })
  }

  // ── 吉神方位 (五项, 2x3 网格) ──────────────────────────────
  @Builder
  directionSection() {
    Column() {
      this.sectionTitle('吉神方位')
      Flex({ wrap: FlexWrap.Wrap, justifyContent: FlexAlign.SpaceBetween }) {
        this.directionItem('喜神', this.almanac!.xiShenFang)
        this.directionItem('财神', this.almanac!.caiShenFang)
        this.directionItem('福神', this.almanac!.fuShenFang)
        this.directionItem('阳贵', this.almanac!.yangGuiFang)
        this.directionItem('阴贵', this.almanac!.yinGuiFang)
      }
      .width('100%')
    }
    .width('100%')
    .padding(AppDimens.PADDING_MD)
    .backgroundColor(AlmanacPalette.BG_SECTION)
    .borderRadius(AppDimens.BORDER_RADIUS_MD)
    .margin({ bottom: AppDimens.PADDING_MD })
  }

  @Builder
  directionItem(label: string, value: string) {
    Column() {
      Text(label)
        .fontSize(AppDimens.FONT_SMALL)
        .fontColor(AlmanacPalette.LABEL)
      Text(value && value.length > 0 ? value : '—')
        .fontSize(AppDimens.FONT_BODY)
        .fontWeight(FontWeight.Medium)
        .fontColor(AppColors.ACCENT)
        .margin({ top: 4 })
    }
    .width('30%')
    .alignItems(HorizontalAlign.Center)
    .padding({ top: 6, bottom: 6 })
  }

  // ── 彭祖百忌 ──────────────────────────────────────────────
  @Builder
  pengZuSection() {
    Column() {
      this.sectionTitle('彭祖百忌')
      this.pengZuLine(this.almanac!.pengZuGan)
      this.pengZuLine(this.almanac!.pengZuZhi)
    }
    .width('100%')
    .padding(AppDimens.PADDING_MD)
    .backgroundColor(AlmanacPalette.BG_SECTION)
    .borderRadius(AppDimens.BORDER_RADIUS_MD)
    .margin({ bottom: AppDimens.PADDING_MD })
  }

  @Builder
  pengZuLine(text: string) {
    if (text.length > 0) {
      Row() {
        Text('·').fontSize(AppDimens.FONT_BODY).fontColor(AppColors.ACCENT)
          .margin({ right: 6 })
        Text(text)
          .fontSize(AppDimens.FONT_BODY)
          .fontColor(AlmanacPalette.VALUE)
          .layoutWeight(1)
      }
      .width('100%')
      .alignItems(VerticalAlign.Top)
      .margin({ top: 2, bottom: 2 })
    }
  }

  // ── 择日典籍语录 (董公择日要诀等) ───────────────────────────
  @Builder
  quotesSection() {
    if (this.almanac!.quotes && this.almanac!.quotes.length > 0) {
      Column() {
        this.sectionTitle('择日典籍')
        ForEach(this.almanac!.quotes, (q: AlmanacQuote, idx: number) => {
          this.quoteCard(q, idx)
        }, (q: AlmanacQuote, idx: number) => `${q.source}_${idx}`)
      }
      .width('100%')
      .padding(AppDimens.PADDING_MD)
      .backgroundColor(AlmanacPalette.BG_SECTION)
      .borderRadius(AppDimens.BORDER_RADIUS_MD)
      .margin({ bottom: AppDimens.PADDING_MD })
    }
  }

  @Builder
  quoteCard(q: AlmanacQuote, idx: number) {
    Column() {
      Row({ space: 8 }) {
        Text(q.title)
          .fontSize(AppDimens.FONT_BODY)
          .fontWeight(FontWeight.Bold)
          .fontColor(AlmanacPalette.VALUE)
          .layoutWeight(1)
        if (q.luck && q.luck.length > 0) {
          Text(q.luck)
            .fontSize(AppDimens.FONT_SMALL)
            .fontColor('#FFFFFF')
            .fontWeight(FontWeight.Medium)
            .backgroundColor(this.luckColor(q.luck))
            .padding({ left: 8, right: 8, top: 2, bottom: 2 })
            .borderRadius(8)
        }
      }
      .width('100%')

      Text(q.source)
        .fontSize(AppDimens.FONT_SMALL)
        .fontColor(AlmanacPalette.LABEL)
        .margin({ top: 4 })

      Text(q.body)
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(AlmanacPalette.SUB)
        .lineHeight(20)
        .margin({ top: 8 })
        .width('100%')
        .textAlign(TextAlign.Start)
    }
    .width('100%')
    .padding(AppDimens.PADDING_SM)
    .backgroundColor(AlmanacPalette.BG_DEEP)
    .borderRadius(AppDimens.BORDER_RADIUS_SM)
    .margin({ top: idx === 0 ? 0 : 8 })
    .alignItems(HorizontalAlign.Start)
  }

  // ── 宜忌 (最显眼, 大色块) ───────────────────────────────
  @Builder
  yiJiSection() {
    if ((this.almanac!.yi && this.almanac!.yi.length > 0)
      || (this.almanac!.ji && this.almanac!.ji.length > 0)) {
      Row({ space: 8 }) {
        // 宜
        Column() {
          Text('宜').fontSize(AppDimens.FONT_TITLE)
            .fontColor('#FFFFFF').fontWeight(FontWeight.Bold)
            .backgroundColor(AlmanacPalette.LUCKY)
            .padding({ left: 12, right: 12, top: 4, bottom: 4 })
            .borderRadius(AppDimens.BORDER_RADIUS_SM)
            .margin({ bottom: 8 })
          Text(this.almanac!.yi.length > 0 ? this.almanac!.yi.join('  ') : '—')
            .fontSize(AppDimens.FONT_BODY)
            .fontColor(AlmanacPalette.VALUE)
            .textAlign(TextAlign.Center)
            .width('100%')
            .lineHeight(22)
        }
        .layoutWeight(1)
        .padding(AppDimens.PADDING_MD)
        .backgroundColor(AlmanacPalette.BG_SECTION)
        .borderRadius(AppDimens.BORDER_RADIUS_MD)
        .alignItems(HorizontalAlign.Center)

        Column() {
          Text('忌').fontSize(AppDimens.FONT_TITLE)
            .fontColor('#FFFFFF').fontWeight(FontWeight.Bold)
            .backgroundColor(AlmanacPalette.UNLUCKY)
            .padding({ left: 12, right: 12, top: 4, bottom: 4 })
            .borderRadius(AppDimens.BORDER_RADIUS_SM)
            .margin({ bottom: 8 })
          Text(this.almanac!.ji.length > 0 ? this.almanac!.ji.join('  ') : '—')
            .fontSize(AppDimens.FONT_BODY)
            .fontColor(AlmanacPalette.VALUE)
            .textAlign(TextAlign.Center)
            .width('100%')
            .lineHeight(22)
        }
        .layoutWeight(1)
        .padding(AppDimens.PADDING_MD)
        .backgroundColor(AlmanacPalette.BG_SECTION)
        .borderRadius(AppDimens.BORDER_RADIUS_MD)
        .alignItems(HorizontalAlign.Center)
      }
      .width('100%')
      .alignItems(VerticalAlign.Top)
      .margin({ bottom: AppDimens.PADDING_MD })
    }
  }

  // ── 神煞 (吉/凶徽章列表) ───────────────────────────────
  @Builder
  shenShaSection() {
    if (this.almanac!.shenSha && this.almanac!.shenSha.length > 0) {
      Column() {
        this.sectionTitle(`神煞 (${this.almanac!.shenSha.length})`)
        Flex({ wrap: FlexWrap.Wrap }) {
          ForEach(this.almanac!.shenSha, (s: ShenSha) => {
            Text(s.name)
              .fontSize(AppDimens.FONT_CAPTION)
              .fontColor('#FFFFFF')
              .fontWeight(s.weight >= 3 ? FontWeight.Bold : FontWeight.Medium)
              .backgroundColor(s.isLucky ? AlmanacPalette.LUCKY : AlmanacPalette.UNLUCKY)
              .padding({ left: 8, right: 8, top: 3, bottom: 3 })
              .borderRadius(10)
              .margin({ right: 6, bottom: 6 })
              .opacity(s.weight === 1 ? 0.78 : 1.0)
          }, (s: ShenSha) => s.name)
        }
        .width('100%')
      }
      .width('100%')
      .padding(AppDimens.PADDING_MD)
      .backgroundColor(AlmanacPalette.BG_SECTION)
      .borderRadius(AppDimens.BORDER_RADIUS_MD)
      .margin({ bottom: AppDimens.PADDING_MD })
      .alignItems(HorizontalAlign.Start)
    }
  }

  // ── 吉曜时法 ────────────────────────────────────────────
  @Builder
  luckyHoursSection() {
    if (this.almanac!.luckyHours && this.almanac!.luckyHours.length > 0) {
      Column() {
        this.sectionTitle('吉曜时法')
        Flex({ wrap: FlexWrap.Wrap, justifyContent: FlexAlign.Start }) {
          ForEach(this.almanac!.luckyHours, (h: LuckyHour) => {
            this.luckyHourCard(h)
          }, (h: LuckyHour) => `${h.name}_${h.zhi}`)
        }
        .width('100%')
      }
      .width('100%')
      .padding(AppDimens.PADDING_MD)
      .backgroundColor(AlmanacPalette.BG_SECTION)
      .borderRadius(AppDimens.BORDER_RADIUS_MD)
      .margin({ bottom: AppDimens.PADDING_MD })
      .alignItems(HorizontalAlign.Start)
    }
  }

  @Builder
  luckyHourCard(h: LuckyHour) {
    Column() {
      Text(h.name)
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(AlmanacPalette.LABEL)
      Text(`${ZHI_NAMES[h.zhi]}时`)
        .fontSize(AppDimens.FONT_BODY)
        .fontWeight(FontWeight.Bold)
        .fontColor(AlmanacPalette.HUANG_DAO)
        .margin({ top: 2 })
      Text(ZHI_HOURS[h.zhi])
        .fontSize(10)
        .fontColor(AlmanacPalette.SUB)
        .margin({ top: 1 })
    }
    .width('31%')
    .alignItems(HorizontalAlign.Center)
    .padding({ top: 6, bottom: 6 })
    .margin({ right: '2%', bottom: 6 })
    .backgroundColor(AlmanacPalette.BG_DEEP)
    .borderRadius(AppDimens.BORDER_RADIUS_SM)
  }

  // ── 用事择吉 ────────────────────────────────────────────
  @Builder
  eventsSection() {
    if (this.almanac!.events && this.almanac!.events.length > 0) {
      Column() {
        this.sectionTitle('用事择吉')
        ForEach(this.almanac!.events, (e: EventAdvice) => {
          this.eventRow(e)
        }, (e: EventAdvice) => e.event)
      }
      .width('100%')
      .padding(AppDimens.PADDING_MD)
      .backgroundColor(AlmanacPalette.BG_SECTION)
      .borderRadius(AppDimens.BORDER_RADIUS_MD)
      .margin({ bottom: AppDimens.PADDING_MD })
      .alignItems(HorizontalAlign.Start)
    }
  }

  @Builder
  eventRow(e: EventAdvice) {
    Row({ space: 8 }) {
      Text(e.suitable ? '宜' : '忌')
        .fontSize(AppDimens.FONT_SMALL)
        .fontColor('#FFFFFF')
        .fontWeight(FontWeight.Bold)
        .backgroundColor(e.suitable ? AlmanacPalette.LUCKY : AlmanacPalette.UNLUCKY)
        .padding({ left: 6, right: 6, top: 2, bottom: 2 })
        .borderRadius(4)
      Text(e.event)
        .fontSize(AppDimens.FONT_BODY)
        .fontColor(AlmanacPalette.VALUE)
        .fontWeight(FontWeight.Medium)
        .width(64)
      Text(e.reason)
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(AlmanacPalette.SUB)
        .layoutWeight(1)
    }
    .width('100%')
    .alignItems(VerticalAlign.Center)
    .margin({ top: 4, bottom: 4 })
  }

  // ── 特别提示 (三煞方位/节气) ────────────────────────────
  @Builder
  notesSection() {
    if (this.almanac!.notes && this.almanac!.notes.length > 0) {
      Column() {
        this.sectionTitle('特别提示')
        ForEach(this.almanac!.notes, (n: string) => {
          Row() {
            Text('·').fontSize(AppDimens.FONT_BODY).fontColor(AppColors.ACCENT)
              .margin({ right: 6 })
            Text(n).fontSize(AppDimens.FONT_CAPTION).fontColor(AlmanacPalette.SUB)
              .layoutWeight(1)
          }
          .width('100%').alignItems(VerticalAlign.Top)
          .margin({ top: 2, bottom: 2 })
        }, (n: string, i: number) => `${i}_${n}`)
      }
      .width('100%')
      .padding(AppDimens.PADDING_MD)
      .backgroundColor(AlmanacPalette.BG_SECTION)
      .borderRadius(AppDimens.BORDER_RADIUS_MD)
      .margin({ bottom: AppDimens.PADDING_MD })
      .alignItems(HorizontalAlign.Start)
    }
  }

  // ── 百科入口 (董公总论/口诀/方位 静态知识) ─────────────
  @Builder
  topicsEntrySection() {
    if (this.onOpenTopics) {
      Row() {
        Text('📖 董公择日 · 总论 · 口诀 · 方位')
          .fontSize(AppDimens.FONT_BODY)
          .fontColor(AppColors.ACCENT)
          .fontWeight(FontWeight.Medium)
        Blank()
        Text('›').fontSize(20).fontColor(AppColors.ACCENT)
      }
      .width('100%')
      .padding(AppDimens.PADDING_MD)
      .backgroundColor(AlmanacPalette.BG_SECTION)
      .borderRadius(AppDimens.BORDER_RADIUS_MD)
      .onClick((): void => {
        const cb = this.onOpenTopics
        if (cb) cb()
      })
      .margin({ bottom: AppDimens.PADDING_MD })
    }
  }

  // 吉/凶/混/平 对应的徽章色
  luckColor(luck: string): string {
    if (luck === '吉') return AlmanacPalette.LUCKY
    if (luck === '凶') return AlmanacPalette.UNLUCKY
    if (luck === '混') return AlmanacPalette.HUANG_DAO
    return AlmanacPalette.HEI_DAO   // 平 / 空
  }

  // ── 章节标题 ──────────────────────────────────────────────
  @Builder
  sectionTitle(t: string) {
    Row() {
      Text(t)
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(AlmanacPalette.LABEL)
        .fontWeight(FontWeight.Medium)
      Blank()
    }
    .width('100%')
    .margin({ bottom: 6 })
  }
}

// ════════════════════════════════════════════════════════════════
//  AlmanacTopicsSheet — 静态知识抽屉
//  数据为常量 (董公总论/口诀/方位), 与日期无关. 由 NativeBridge.getAlmanacTopics 一次取回.
// ════════════════════════════════════════════════════════════════

@Component
export struct AlmanacTopicsSheet {
  topics: AlmanacTopic[] = []

  build() {
    Scroll() {
      Column({ space: AppDimens.PADDING_MD }) {
        Text('董公择日 · 经典知识')
          .fontSize(AppDimens.FONT_SUBTITLE)
          .fontWeight(FontWeight.Bold)
          .fontColor(AlmanacPalette.VALUE)
          .margin({ top: AppDimens.PADDING_MD, bottom: 4 })

        if (this.topics.length === 0) {
          Text('数据加载中...')
            .fontSize(AppDimens.FONT_BODY)
            .fontColor(AlmanacPalette.SUB)
            .padding(40)
        }

        ForEach(this.topics, (t: AlmanacTopic) => {
          this.topicCard(t)
        }, (t: AlmanacTopic) => `${t.category}_${t.title}`)
      }
      .width('100%')
      .padding({ left: AppDimens.PADDING_LG, right: AppDimens.PADDING_LG,
                 bottom: AppDimens.PADDING_XL })
    }
    .width('100%')
    .height('100%')
    .scrollBar(BarState.Auto)
    .backgroundColor(AlmanacPalette.BG_DEEP)
  }

  @Builder
  topicCard(t: AlmanacTopic) {
    Column() {
      Row({ space: 6 }) {
        Text(t.category)
          .fontSize(AppDimens.FONT_SMALL)
          .fontColor('#FFFFFF')
          .backgroundColor(AppColors.ACCENT)
          .padding({ left: 6, right: 6, top: 1, bottom: 1 })
          .borderRadius(4)
        Text(t.title)
          .fontSize(AppDimens.FONT_BODY)
          .fontWeight(FontWeight.Bold)
          .fontColor(AlmanacPalette.VALUE)
          .layoutWeight(1)
      }
      .width('100%').alignItems(VerticalAlign.Center)

      Text(t.body)
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(AlmanacPalette.SUB)
        .lineHeight(22)
        .margin({ top: 8 })
        .width('100%')
        .textAlign(TextAlign.Start)
    }
    .width('100%')
    .padding(AppDimens.PADDING_MD)
    .backgroundColor(AlmanacPalette.BG_SECTION)
    .borderRadius(AppDimens.BORDER_RADIUS_MD)
    .alignItems(HorizontalAlign.Start)
  }
}
