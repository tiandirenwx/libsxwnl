import { AppColors, AppDimens } from '../common/Constants'
import { BaziResult, BaziColumn, CangGanItem, LiuNianItem } from '../services/NativeBridge'
import { GAN, ZHI, BaziCalc } from '../model/BaziInput'
import { image } from '@kit.ImageKit'
import { fileIo as fs, fileUri } from '@kit.CoreFileKit'
import { common } from '@kit.AbilityKit'
import { photoAccessHelper } from '@kit.MediaLibraryKit'
import { systemShare } from '@kit.ShareKit'
import { uniformTypeDescriptor as utd } from '@kit.ArkData'

/** 结果页入参 (class 形式, 满足 ArkTS 严格类型规则) */
export class BaziResultArg {
  result: BaziResult = emptyBaziResult()
  birthYear: number = 1990
  astEnabled: boolean = false
  longitude: number = 116.4
  latitude: number = 39.9
  lifaLabel: string = '定气'
}

/** 构造一个全空 BaziResult, 用于 @Component 默认值占位. */
function emptyBaziResult(): BaziResult {
  return {
    userName: '', gender: '', solarBirth: '', lunarBirth: '', dateOfBirth: '',
    shengXiao: '', age: '', lifa: '', dingQiType: '', ast: '', jieQi: '',
    qiYun: '', jiaoYun: '',
    columns: [], daYunColumns: [], wuXingCount: [], wuXingStatus: [],
    siLing: '', liuNian: []
  } as BaziResult
}

/**
 * 八字排盘结果视图: 顶部切换 简洁版 / 专业版。
 * 简洁版严格参考"宣纸命书"排版(参考图四)。
 */
@Component
export struct BaziResultView {
  arg: BaziResultArg = new BaziResultArg()
  @State proMode: boolean = false
  /** 专业版"大运详情"选中的下标 (默认= 当前大运在 daYunColumns 里的位置) */
  @State selectedDaYunIdx: number = -1

  private readonly colTitles: string[] = ['年柱', '月柱', '日柱', '时柱', '大运', '流年']
  private readonly elemColors: string[] = ['#2E7D32', '#D32F2F', '#8D6E40', '#C9A227', '#1565C0']
  private readonly LABEL_W: number = 32
  // 简洁版配色(宣纸命书: 墨色为主, 朱红高亮)
  private readonly inkColor: string = '#3A2A1C'
  private readonly inkSoft: string = '#6B5640'
  private readonly redColor: string = '#C0392B'
  private readonly redBg: string = '#F6DAD2'
  private readonly goldLine: string = '#C9A86A'

  build() {
    Scroll() {
      Column() {
        this.versionToggle()
        if (this.proMode) {
          // 专业版命盘整体套一层带 id 的容器供截图(保存/分享)
          Column() {
            this.infoBar()
            this.paiPanTable()
            this.summaryBar()
            this.liuNianCard()
          }
          .width('100%').id('baziProCard').backgroundColor(AppColors.BACKGROUND)
          this.saveShareRow('baziProCard')
        } else {
          this.simpleCard()
        }
      }
      .width('100%').padding({ bottom: 12 })
    }
    .width('100%').height('100%')
    .backgroundColor(AppColors.BACKGROUND)
    .scrollBar(BarState.Off)
    .align(Alignment.Top)
  }

  // ───────────────────── 版本切换 ─────────────────────
  @Builder
  versionToggle() {
    Row() {
      Text('简洁版')
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(!this.proMode ? AppColors.ON_PRIMARY : AppColors.TEXT_SECONDARY)
        .backgroundColor(!this.proMode ? AppColors.PRIMARY : AppColors.SURFACE)
        .textAlign(TextAlign.Center).width(90).height(34)
        .borderRadius({ topLeft: 17, bottomLeft: 17 })
        .onClick(() => { this.proMode = false })
      Text('专业版')
        .fontSize(AppDimens.FONT_CAPTION)
        .fontColor(this.proMode ? AppColors.ON_PRIMARY : AppColors.TEXT_SECONDARY)
        .backgroundColor(this.proMode ? AppColors.PRIMARY : AppColors.SURFACE)
        .textAlign(TextAlign.Center).width(90).height(34)
        .borderRadius({ topRight: 17, bottomRight: 17 })
        .onClick(() => { this.proMode = true })
    }
    .margin({ top: AppDimens.PADDING_MD })
    .border({ width: 1, color: AppColors.DIVIDER, radius: 17 })
  }

  // ───────────────────── 工具 ─────────────────────
  // 去掉 『xxx』 标签前缀
  strip(s: string): string {
    if (!s) return ''
    const i: number = s.indexOf('』')
    return i >= 0 ? s.substring(i + 1).trim() : s.trim()
  }

  // 干支文字 → 五行颜色(专业版)
  charColor(ch: string): string {
    const gi: number = '甲乙丙丁戊己庚辛壬癸'.indexOf(ch)
    if (gi >= 0) return this.elemColors[Math.floor(gi / 2)]
    const zi: number = '子丑寅卯辰巳午未申酉戌亥'.indexOf(ch)
    if (zi >= 0) {
      const zw: number[] = [4, 2, 0, 0, 2, 1, 1, 2, 3, 3, 2, 4]
      return this.elemColors[zw[zi]]
    }
    return AppColors.ON_SURFACE
  }

  yearZhiIndex(): number {
    if (this.arg.result.columns === undefined || this.arg.result.columns.length === 0) return 0
    const idx: number = ZHI.indexOf(this.arg.result.columns[0].zhi)
    return idx >= 0 ? idx : 0
  }

  zodiacRes(): Resource {
    const arr: Resource[] = [
      $r('app.media.bz_rat'), $r('app.media.bz_cow'), $r('app.media.bz_tiger'),
      $r('app.media.bz_rabbit'), $r('app.media.bz_dragon'), $r('app.media.bz_snake'),
      $r('app.media.bz_horse'), $r('app.media.bz_goat'), $r('app.media.bz_monkey'),
      $r('app.media.bz_hen'), $r('app.media.bz_dog'), $r('app.media.bz_pig')
    ]
    return arr[this.yearZhiIndex()]
  }

  curDaYunStart(): number {
    const d: BaziColumn | undefined = this.arg.result.currentDaYun
    return d !== undefined ? d.startYear : -1
  }
  curLiuNianYear(): number {
    const l: BaziColumn | undefined = this.arg.result.currentLiuNian
    return l !== undefined ? l.startYear : -1
  }
  daYunStartAge(startYear: number): number { return startYear - this.arg.birthYear + 1 }
  decadeYears(start: number): number[] {
    const r: number[] = []
    for (let j = 0; j < 10; j++) r.push(start + j)
    return r
  }

  buildColumns(): BaziColumn[] {
    const list: BaziColumn[] = []
    if (this.arg.result.columns !== undefined) {
      for (let i = 0; i < this.arg.result.columns.length; i++) list.push(this.arg.result.columns[i])
    }
    if (this.arg.result.currentDaYun !== undefined) list.push(this.arg.result.currentDaYun)
    if (this.arg.result.currentLiuNian !== undefined) list.push(this.arg.result.currentLiuNian)
    return list
  }

  // ════════════════════ 简洁版(宣纸命书) ════════════════════
  @Builder
  simpleCard() {
    Column() {
      this.simpleReport()
      this.saveShareRow('baziSimpleCard')
    }
    .width('100%').alignItems(HorizontalAlign.Center)
  }

  // 保存为图片 + 分享到 (简洁版/专业版共用)
  //  · 左: 安全控件 SaveButton, 点击后系统授予一次性写相册授权, 无需声明/申请权限;
  //  · 右: 普通按钮调起系统分享面板, 分享到其它 app。
  @Builder
  saveShareRow(snapshotId: string) {
    Row() {
      SaveButton({ icon: SaveIconStyle.FULL_FILLED, text: SaveDescription.SAVE_IMAGE, buttonType: ButtonType.Capsule })
        .layoutWeight(1).height(40).fontSize(AppDimens.FONT_BODY).fontWeight(FontWeight.Medium)
        .fontColor(AppColors.PRIMARY_DARK).iconColor(AppColors.PRIMARY_DARK)
        .backgroundColor(AppColors.ACCENT)
        .borderRadius(AppDimens.BORDER_RADIUS_LG)
        .onClick((_event: ClickEvent, result: SaveButtonOnClickResult) => {
          if (result === SaveButtonOnClickResult.SUCCESS) {
            this.saveImageToGallery(snapshotId)
          } else {
            this.toast('没有相册写入授权')
          }
        })
      Button('分享到')
        .layoutWeight(1).height(40).fontSize(AppDimens.FONT_BODY).fontWeight(FontWeight.Medium)
        .fontColor(AppColors.ON_PRIMARY).backgroundColor(AppColors.PRIMARY)
        .borderRadius(AppDimens.BORDER_RADIUS_LG)
        .margin({ left: 12 })
        .onClick(() => { this.shareImage(snapshotId) })
    }
    .width('100%').padding({ left: 24, right: 24 })
    .margin({ top: AppDimens.PADDING_SM, bottom: AppDimens.PADDING_SM })
  }

  // 单个『键』值字段 (紧凑, 单行不换行; 键固定, 值可尾部省略)
  @Builder
  infoField(k: string, v: string, red: boolean, weight: number) {
    Row() {
      Text(`『${k}』`).fontSize(12).fontColor(this.inkSoft).flexShrink(0)
      Text(v).fontSize(12).fontColor(red ? this.redColor : this.inkColor)
        .margin({ left: 1 }).layoutWeight(1).maxLines(1)
        .textOverflow({ overflow: TextOverflow.Ellipsis })
    }.layoutWeight(weight).alignItems(VerticalAlign.Top)
  }

  // 四字段一行 (姓名/性别/生肖/年龄), 对齐命盘图紧凑排版
  @Builder
  infoQuad(k1: string, v1: string, k2: string, v2: string,
    k3: string, v3: string, k4: string, v4: string) {
    Row() {
      this.infoField(k1, v1, false, 13)
      this.infoField(k2, v2, false, 9)
      this.infoField(k3, v3, false, 9)
      this.infoField(k4, v4, false, 11)
    }.width('100%').margin({ bottom: 1 }).alignItems(VerticalAlign.Top)
  }

  // 两字段一行 (出生地/时间标准)
  @Builder
  infoDuo(k1: string, v1: string, red1: boolean, k2: string, v2: string, red2: boolean) {
    Row() {
      this.infoField(k1, v1, red1, 11)
      this.infoField(k2, v2, red2, 10)
    }.width('100%').margin({ bottom: 1 }).alignItems(VerticalAlign.Top)
  }

  @Builder
  simpleReport() {
    // 结构 (完全等价于"改动前"稳定版, 只在其内部嵌入生肖印图层):
    //   外层 Column: 决定卡片高度 + 挂宣纸背景 (.backgroundImage, 不占布局)
    //   内层 Stack: 生肖印 (position 居中, 脱离布局) 与 文字 Column 并列
    // 层级: 宣纸(Column.background) → 生肖(Stack 第一子) → 文字(Stack 第二子, 在上)
    // 高度由文字 Column 自然决定, 不会被宣纸/生肖撑高.
    // 生肖 PNG 本身即水墨风格 (大量透明像素), opacity=1.0 用图自身 alpha 混合,
    // 视觉与 Python img.paste(zodiac, (28,150), mask=zodiac) 一致.
    Column() {
      Stack({ alignContent: Alignment.TopStart }) {
        // ① 生肖印 (以自身中心为锚点, 定位到 Stack 中心)
        //   ArkUI 里 .align() 只影响子内部内容对齐, 不影响子在父 Stack 中的位置;
        //   必须用 .position({x:'50%',y:'50%'}) + .markAnchor({x:'50%',y:'50%'}).
        Image(this.zodiacRes())
          .width(200).height(200)
          .objectFit(ImageFit.Contain)
          .opacity(1.0)
          .position({ x: '50%', y: '50%' })
          .markAnchor({ x: '50%', y: '50%' })

        // ② 命书正文 (紧凑排布, 力求单屏完整显示; 覆盖生肖印之上)
        Column() {
          // 抬头品牌 (单行, 不占额外竖向空间)
          Row() {
            Blank()
            Text('☯ 八字命书').fontSize(12).fontColor(this.redColor).fontWeight(FontWeight.Bold)
          }.width('100%').margin({ bottom: 4 })

          // ── 命主信息 (姓名/性别/生肖/年龄 合并为一行) ──
          this.infoQuad('姓名', this.strip(this.arg.result.userName),
            '性别', this.arg.result.gender !== undefined && this.arg.result.gender.indexOf('女') >= 0 ? '女' : '男',
            '生肖', this.strip(this.arg.result.shengXiao),
            '年龄', this.strip(this.arg.result.age))
          this.infoDuo('出生地', `东经${this.arg.longitude.toFixed(2)}° 北纬${this.arg.latitude.toFixed(2)}°`, false,
            '时间标准', this.arg.astEnabled ? '真太阳时' : '北京时间(120°E)', false)
          // 公历 / 农历各占整行 (内容较长)
          Row() {
            Text('『公历』').fontSize(12).fontColor(this.inkSoft).flexShrink(0)
            Text(this.strip(this.arg.result.solarBirth)).fontSize(12).fontColor(this.inkColor)
              .margin({ left: 1 }).layoutWeight(1).maxLines(1).textOverflow({ overflow: TextOverflow.Ellipsis })
          }.width('100%').margin({ bottom: 1 })
          Row() {
            Text('『农历』').fontSize(12).fontColor(this.inkSoft).flexShrink(0)
            Text(this.strip(this.arg.result.lunarBirth)).fontSize(12).fontColor(this.redColor)
              .margin({ left: 1 }).layoutWeight(1).maxLines(1).textOverflow({ overflow: TextOverflow.Ellipsis })
          }.width('100%').margin({ bottom: 1 })

          // ── 节气/司令 段落 ──
          Text(this.jieQiParagraph())
            .fontSize(11).fontColor(this.inkSoft).width('100%').lineHeight(15).margin({ top: 3 })

          Divider().color(this.goldLine).opacity(0.6).margin({ top: 4, bottom: 4 })

          // ── 四柱(主星/天干/地支/藏干) ──
          this.siZhuBlock()

          Divider().color(this.goldLine).opacity(0.6).margin({ top: 4, bottom: 4 })

          // ── 起运/交运 段落 ──
          Text(`${this.strip(this.arg.result.qiYun)}　${this.strip(this.arg.result.jiaoYun)}`)
            .fontSize(11).fontColor(this.inkSoft).width('100%').lineHeight(15)

          // ── 大运 / 流年 网格 (大運标签内嵌于表格左侧) ──
          this.daYunGrid()
        }.width('100%')
      }
      .width('100%')
    }
    .width('100%').padding({ left: 12, right: 12, top: 10, bottom: 10 })
    .backgroundImage($r('app.media.bz_paper'), ImageRepeat.NoRepeat)
    .backgroundImageSize(ImageSize.Cover)
    .id('baziSimpleCard')
    .width('100%')
    .margin({ left: AppDimens.PADDING_SM, right: AppDimens.PADDING_SM, top: AppDimens.PADDING_SM })
    .borderRadius(10).clip(true)
    .border({ width: 1, color: this.goldLine })
    .shadow({ radius: 8, color: AppColors.CARD_SHADOW_COLOR, offsetY: 2 })
  }

  jieQiParagraph(): string {
    let s: string = ''
    const dq: string = this.strip(this.arg.result.dingQiType)
    if (dq.length > 0) s += `依据${dq}。`
    const jq: string = this.arg.result.jieQi.trim().replace(/\n+/g, '；')
    if (jq.length > 0) s += jq + '。'
    if (this.arg.result.siLing !== undefined && this.arg.result.siLing.length > 0) s += `命主月令司令：${this.arg.result.siLing}`
    return s
  }

  @Builder
  siZhuBlock() {
    Row() {
      // 乾造/坤造 (顶部留白, 与天干行对齐, 跳过主星行)
      Text(this.arg.result.gender !== undefined && this.arg.result.gender.indexOf('女') >= 0 ? '坤造' : '乾造')
        .fontSize(12).fontFamily('WenYue').fontColor(this.inkColor).width(24).margin({ top: 16 })
      ForEach(this.arg.result.columns, (c: BaziColumn, i: number) => {
        Column() {
          Text(i === 2 ? '日元' : c.ganShiShen)
            .fontSize(11).fontColor(i === 2 ? this.redColor : this.inkSoft).fontWeight(FontWeight.Medium)
          Text(c.gan).fontSize(18).fontFamily('WenYue')
            .fontColor(i === 2 ? this.redColor : this.inkColor).margin({ top: 1 })
          Text(c.zhi).fontSize(18).fontFamily('WenYue').fontColor(this.inkColor).margin({ top: 1 })
          Column() {
            ForEach(c.cangGan, (cg: CangGanItem) => {
              Row() {
                Text(cg.gan).fontSize(10).fontFamily('WenYue').fontColor(this.inkColor)
                Text(cg.shiShen).fontSize(10).fontColor(this.inkSoft).margin({ left: 2 })
              }
            })
          }.margin({ top: 3 }).alignItems(HorizontalAlign.Center)
        }
        .layoutWeight(1).alignItems(HorizontalAlign.Center)
      })
    }
    .width('100%').alignItems(VerticalAlign.Top)
  }

  // 仅取前 8 步大运(与参考命书一致)
  daYunColumns8(): BaziColumn[] {
    const all: BaziColumn[] = this.arg.result.daYunColumns !== undefined ? this.arg.result.daYunColumns : []
    return all.length > 8 ? all.slice(0, 8) : all
  }

  // 大运 / 流年: 左侧"大運"标签 + 8 列铺满(不滚动), 当年大运/流年仅标红不加框, 流年横排
  @Builder
  daYunGrid() {
    Row() {
      // 大運 标签 (顶部留白, 与干支行对齐, 跳过起运岁数行)
      Text('大運').fontSize(12).fontFamily('WenYue').fontColor(this.inkColor)
        .width(22).margin({ top: 18 })
      ForEach(this.daYunColumns8(), (c: BaziColumn) => {
        this.daYunColumn(c)
      })
    }
    .width('100%').alignItems(VerticalAlign.Top).margin({ top: 3 })
  }

  @Builder
  daYunColumn(c: BaziColumn) {
    Column() {
      // 起运虚岁
      Text(`${this.daYunStartAge(c.startYear)}`).fontSize(10)
        .fontColor(c.startYear === this.curDaYunStart() ? this.redColor : this.inkSoft)
        .fontWeight(FontWeight.Medium)
      // 大运干支(上下堆叠; 当年大运标红+底色, 无边框)
      Column() {
        Text(c.gan).fontSize(16).fontFamily('WenYue')
          .fontColor(c.startYear === this.curDaYunStart() ? this.redColor : this.inkColor)
        Text(c.zhi).fontSize(16).fontFamily('WenYue')
          .fontColor(c.startYear === this.curDaYunStart() ? this.redColor : this.inkColor)
      }
      .padding({ left: 2, right: 2, top: 1, bottom: 1 })
      .backgroundColor(c.startYear === this.curDaYunStart() ? this.redBg : Color.Transparent)
      .borderRadius(3)
      // 起始年
      Text(`${c.startYear}`).fontSize(8).fontColor(this.redColor).margin({ top: 1 })
      Divider().color(this.goldLine).opacity(0.5).margin({ top: 2, bottom: 2 }).width('86%')
      // 10 个流年(横排: 干支同一行; 当年流年标红)
      Column() {
        ForEach(this.decadeYears(c.startYear), (y: number) => {
          Text(GAN[BaziCalc.lnGan(y)] + ZHI[BaziCalc.lnZhi(y)])
            .fontSize(12).fontFamily('WenYue').maxLines(1)
            .fontColor(y === this.curLiuNianYear() ? this.redColor : this.inkColor)
            .padding({ left: 1, right: 1 })
            .backgroundColor(y === this.curLiuNianYear() ? this.redBg : Color.Transparent)
            .borderRadius(3)
        })
      }
      // 止于年
      Text(`${c.startYear + 9}`).fontSize(8).fontColor(this.redColor).margin({ top: 1 })
    }
    .layoutWeight(1).alignItems(HorizontalAlign.Center).padding({ left: 1, right: 1 })
  }

  toast(msg: string, duration: number = 2000) {
    try {
      this.getUIContext().getPromptAction().showToast({ message: msg, duration: duration })
    } catch (e) {
      console.error('showToast failed: ' + JSON.stringify(e))
    }
  }

  // 截图并保存到系统相册。
  //  · 必须在 SaveButton(安全控件)的 onClick 回调内、且授权成功后调用, 由系统授予一次性写相册权限;
  //  · ImagePacker / PixelMap / 文件描述符均确保释放, 避免 native 资源泄漏。
  async saveImageToGallery(snapshotId: string) {
    let pm: image.PixelMap | undefined
    const packer: image.ImagePacker = image.createImagePacker()
    try {
      pm = await this.getUIContext().getComponentSnapshot().get(snapshotId)
      const opts: image.PackingOption = { format: 'image/png', quality: 100 }
      const buf: ArrayBuffer = await packer.packToData(pm, opts)

      const ctx = this.getUIContext().getHostContext() as common.UIAbilityContext
      const helper = photoAccessHelper.getPhotoAccessHelper(ctx)
      // createAsset(photoType, extension): 配合安全控件时无需权限即可创建相册资源
      const uri: string = await helper.createAsset(photoAccessHelper.PhotoType.IMAGE, 'png')
      const f = await fs.open(uri, fs.OpenMode.READ_WRITE | fs.OpenMode.CREATE)
      try {
        await fs.write(f.fd, buf)
      } finally {
        await fs.close(f.fd)
      }
      this.toast('已保存到相册')
    } catch (e) {
      console.error('saveImageToGallery failed: ' + JSON.stringify(e))
      this.toast('保存失败')
    } finally {
      packer.release()
      if (pm !== undefined) {
        await pm.release()
      }
    }
  }

  // 截图并调起系统分享面板, 分享到其它 app。
  //  · 快照 → PNG → 写入沙箱缓存 → fileUri 转可分享 uri → systemShare.ShareController;
  //  · 分享无需相册授权, 用普通按钮触发即可。
  async shareImage(snapshotId: string) {
    let pm: image.PixelMap | undefined
    const packer: image.ImagePacker = image.createImagePacker()
    try {
      pm = await this.getUIContext().getComponentSnapshot().get(snapshotId)
      const opts: image.PackingOption = { format: 'image/png', quality: 100 }
      const buf: ArrayBuffer = await packer.packToData(pm, opts)

      const ctx = this.getUIContext().getHostContext() as common.UIAbilityContext
      const dir: string = ctx.cacheDir + '/shares'
      try {
        fs.mkdirSync(dir)
      } catch (_) {
        // 目录已存在时忽略
      }
      const path: string = dir + '/bazi_' + Date.now() + '.png'
      const file = fs.openSync(path, fs.OpenMode.READ_WRITE | fs.OpenMode.CREATE | fs.OpenMode.TRUNC)
      try {
        fs.writeSync(file.fd, buf)
      } finally {
        fs.closeSync(file)
      }

      const uri: string = fileUri.getUriFromPath(path)
      const shareData: systemShare.SharedData = new systemShare.SharedData({
        utd: utd.UniformDataType.PNG,
        uri: uri,
        title: '八字命盘'
      })
      const controller: systemShare.ShareController = new systemShare.ShareController(shareData)
      await controller.show(ctx, {
        selectionMode: systemShare.SelectionMode.SINGLE,
        previewMode: systemShare.SharePreviewMode.DETAIL
      })
    } catch (e) {
      console.error('shareImage failed: ' + JSON.stringify(e))
      this.toast('分享失败')
    } finally {
      packer.release()
      if (pm !== undefined) {
        await pm.release()
      }
    }
  }

  // ════════════════════ 专业版 ════════════════════
  @Builder
  infoBar() {
    Column() {
      Text(`${this.arg.result.gender !== undefined && this.arg.result.gender.indexOf('女') >= 0 ? '坤' : '乾'}  ${this.strip(this.arg.result.userName)}`)
        .fontSize(AppDimens.FONT_SUBTITLE).fontWeight(FontWeight.Bold).fontColor(AppColors.PRIMARY)
      Text(`${this.strip(this.arg.result.shengXiao)}  ${this.strip(this.arg.result.age)}`)
        .fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.TEXT_SECONDARY).margin({ top: 2 })
      Text(`公历 ${this.strip(this.arg.result.solarBirth)}`)
        .fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.TEXT_SECONDARY).margin({ top: 2 })
      Text(`农历 ${this.strip(this.arg.result.lunarBirth)}`)
        .fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.TEXT_SECONDARY).margin({ top: 2 })
    }
    .width('100%').alignItems(HorizontalAlign.Start)
    .padding(AppDimens.PADDING_MD)
    .margin({ left: AppDimens.PADDING_MD, right: AppDimens.PADDING_MD, top: AppDimens.PADDING_MD })
    .backgroundColor(AppColors.SURFACE).borderRadius(AppDimens.BORDER_RADIUS_MD)
  }

  @Builder
  paiPanTable() {
    // 四柱 + 大运 + 流年 = 6 列, 取消横向滚动, 各列 layoutWeight 等分填满
    Column() {
      this.headerRow()
      this.starRow()
      this.ganRow()
      this.zhiRow()
      this.cangGanRow()
      this.textRow('纳音', this.rowValues(0))
      this.textRow('星运', this.rowValues(1))
      this.textRow('自坐', this.rowValues(2))
      this.textRow('空亡', this.rowValues(3))
      this.shenShaRow()
    }
    .width('100%').alignItems(HorizontalAlign.Center)
    .padding(AppDimens.PADDING_SM)
    .margin({ left: AppDimens.PADDING_MD, right: AppDimens.PADDING_MD, top: AppDimens.PADDING_MD })
    .backgroundColor(AppColors.SURFACE).borderRadius(AppDimens.BORDER_RADIUS_MD)
  }

  @Builder
  cell(label: string) {
    Text(label).fontSize(AppDimens.FONT_SMALL).fontColor(AppColors.TEXT_SECONDARY)
      .textAlign(TextAlign.Center).width(this.LABEL_W)
  }

  @Builder
  headerRow() {
    Row() {
      Text('').width(this.LABEL_W)
      ForEach(this.colTitles, (t: string, i: number) => {
        Text(t).fontSize(AppDimens.FONT_CAPTION)
          .fontColor(i >= 4 ? AppColors.ACCENT : AppColors.PRIMARY)
          .fontWeight(FontWeight.Bold).textAlign(TextAlign.Center).layoutWeight(1)
      })
    }
    .width('100%').height(28).alignItems(VerticalAlign.Center)
    .border({ width: { bottom: 1 }, color: AppColors.DIVIDER })
  }

  @Builder
  starRow() {
    Row() {
      this.cell('主星')
      ForEach(this.buildColumns(), (c: BaziColumn, i: number) => {
        Text(c.ganShiShen).fontSize(AppDimens.FONT_CAPTION)
          .fontColor(i === 2 ? AppColors.PRIMARY : AppColors.TEXT_SECONDARY)
          .fontWeight(i === 2 ? FontWeight.Bold : FontWeight.Normal)
          .textAlign(TextAlign.Center).layoutWeight(1)
      })
    }
    .width('100%').height(26).alignItems(VerticalAlign.Center)
  }

  @Builder
  ganRow() {
    Row() {
      this.cell('天干')
      ForEach(this.buildColumns(), (c: BaziColumn) => {
        Text(c.gan).fontSize(24).fontWeight(FontWeight.Bold)
          .fontColor(this.charColor(c.gan)).textAlign(TextAlign.Center).layoutWeight(1)
      })
    }
    .width('100%').height(40).alignItems(VerticalAlign.Center)
  }

  @Builder
  zhiRow() {
    Row() {
      this.cell('地支')
      ForEach(this.buildColumns(), (c: BaziColumn) => {
        Text(c.zhi).fontSize(24).fontWeight(FontWeight.Bold)
          .fontColor(this.charColor(c.zhi)).textAlign(TextAlign.Center).layoutWeight(1)
      })
    }
    .width('100%').height(40).alignItems(VerticalAlign.Center)
  }

  @Builder
  cangGanRow() {
    Row() {
      this.cell('藏干')
      ForEach(this.buildColumns(), (c: BaziColumn) => {
        Column() {
          ForEach(c.cangGan, (cg: CangGanItem) => {
            Row() {
              Text(cg.gan).fontSize(AppDimens.FONT_SMALL).fontColor(this.charColor(cg.gan))
              Text(cg.shiShen).fontSize(AppDimens.FONT_SMALL).fontColor(AppColors.TEXT_SECONDARY).margin({ left: 2 })
            }
          })
        }
        .layoutWeight(1).alignItems(HorizontalAlign.Center).justifyContent(FlexAlign.Center)
      })
    }
    .width('100%').alignItems(VerticalAlign.Top).padding({ top: 4, bottom: 4 })
  }

  rowValues(field: number): string[] {
    const out: string[] = []
    const cols: BaziColumn[] = this.buildColumns()
    for (let i = 0; i < cols.length; i++) {
      const c: BaziColumn = cols[i]
      if (field === 0) out.push(c.nayin)
      else if (field === 1) out.push(c.xingYun)
      else if (field === 2) out.push(c.ziZuo)
      else out.push(c.kongWang)
    }
    return out
  }

  @Builder
  textRow(label: string, values: string[]) {
    Row() {
      this.cell(label)
      ForEach(values, (v: string) => {
        Text(v).fontSize(AppDimens.FONT_SMALL).fontColor(AppColors.ON_SURFACE)
          .textAlign(TextAlign.Center).layoutWeight(1)
      })
    }
    .width('100%').height(24).alignItems(VerticalAlign.Center)
    .border({ width: { top: 0.5 }, color: AppColors.DIVIDER })
  }

  @Builder
  shenShaRow() {
    Row() {
      this.cell('神煞')
      ForEach(this.buildColumns(), (c: BaziColumn) => {
        Column() {
          ForEach(c.shenSha, (s: string) => {
            Text(s).fontSize(9).fontColor(AppColors.TEXT_SECONDARY).margin({ bottom: 1 })
          })
        }
        .layoutWeight(1).alignItems(HorizontalAlign.Center)
      })
    }
    .width('100%').alignItems(VerticalAlign.Top).padding({ top: 4, bottom: 4 })
    .border({ width: { top: 0.5 }, color: AppColors.DIVIDER })
  }

  @Builder
  summaryBar() {
    Column() {
      Row() {
        Text('司令').fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.TEXT_SECONDARY)
        Text(this.arg.result.siLing).fontSize(AppDimens.FONT_BODY).fontWeight(FontWeight.Bold)
          .fontColor(AppColors.PRIMARY).margin({ left: 4 })
        Blank()
        Text(this.wuXingSummary()).fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.ON_SURFACE)
      }
      .width('100%').margin({ bottom: AppDimens.PADDING_SM })

      Row() {
        ForEach(this.statusChips(), (item: string, i: number) => {
          Text(item).fontSize(AppDimens.FONT_CAPTION).fontColor(AppColors.ON_PRIMARY)
            .textAlign(TextAlign.Center).layoutWeight(1).height(28)
            .backgroundColor(this.elemColors[i]).borderRadius(AppDimens.BORDER_RADIUS_SM)
            .margin({ right: i < 4 ? 4 : 0 })
        })
      }
      .width('100%')
    }
    .width('100%').padding(AppDimens.PADDING_MD)
    .margin({ left: AppDimens.PADDING_MD, right: AppDimens.PADDING_MD, top: AppDimens.PADDING_MD })
    .backgroundColor(AppColors.SURFACE).borderRadius(AppDimens.BORDER_RADIUS_MD)
  }

  wuXingSummary(): string {
    if (this.arg.result.wuXingCount === undefined) return ''
    const names: string[] = ['木', '火', '土', '金', '水']
    let s: string = '含藏气: '
    for (let i = 0; i < 5; i++) s += `${this.arg.result.wuXingCount[i]}${names[i]} `
    return s.trim()
  }

  statusChips(): string[] {
    if (this.arg.result.wuXingStatus === undefined) return ['', '', '', '', '']
    const names: string[] = ['木', '火', '土', '金', '水']
    const out: string[] = []
    for (let i = 0; i < 5; i++) out.push(`${names[i]}${this.arg.result.wuXingStatus[i]}`)
    return out
  }

  // ── 当前选中大运 (默认为命主当前大运) ────────────────────
  curDaYunIdx(): number {
    if (this.selectedDaYunIdx >= 0) return this.selectedDaYunIdx
    const cur: number = this.curDaYunStart()
    if (this.arg.result.daYunColumns !== undefined) {
      for (let i = 0; i < this.arg.result.daYunColumns.length; i++) {
        if (this.arg.result.daYunColumns[i].startYear === cur) return i
      }
    }
    return 0
  }

  // 当前选中大运对应的 10 年流年/小运
  curLiuNianList(): LiuNianItem[] {
    const idx: number = this.curDaYunIdx()
    const all: LiuNianItem[][] | undefined = this.arg.result.liuNianAll
    if (all !== undefined && idx >= 0 && idx < all.length) return all[idx]
    return this.arg.result.liuNian !== undefined ? this.arg.result.liuNian : []
  }

  // 当前选中大运是否为命主"实际"大运 (用于流年红框)
  isSelectedCurDaYun(): boolean {
    const dy = this.arg.result.daYunColumns
    if (dy === undefined || dy.length === 0) return false
    return dy[this.curDaYunIdx()].startYear === this.curDaYunStart()
  }

  @Builder
  liuNianCard() {
    Column() {
      Text('大运 (点击切换查看流年 / 小运)')
        .fontSize(AppDimens.FONT_SUBTITLE).fontWeight(FontWeight.Bold)
        .fontColor(AppColors.PRIMARY).width('100%').margin({ bottom: AppDimens.PADDING_SM })
      Scroll() {
        Row() {
          ForEach(this.arg.result.daYunColumns, (c: BaziColumn, i: number) => {
            Column() {
              Text(c.ganShiShen).fontSize(9).fontColor(AppColors.TEXT_SECONDARY)
              Text(c.gan).fontSize(18).fontWeight(FontWeight.Bold).fontColor(this.charColor(c.gan))
              Text(c.zhi).fontSize(18).fontWeight(FontWeight.Bold).fontColor(this.charColor(c.zhi))
              Text(`${c.startYear}`).fontSize(9).fontColor(AppColors.TEXT_SECONDARY).margin({ top: 2 })
            }
            .width(48).padding({ top: 4, bottom: 4 }).margin({ right: 4 })
            .alignItems(HorizontalAlign.Center)
            .backgroundColor(i === this.curDaYunIdx() ? this.redBg : AppColors.BACKGROUND)
            .borderRadius(AppDimens.BORDER_RADIUS_SM)
            .border(c.startYear === this.curDaYunStart()
              ? { width: 1, color: this.redColor }
              : { width: 0 })
            .onClick(() => { this.selectedDaYunIdx = i })
          })
        }
      }
      .scrollable(ScrollDirection.Horizontal).scrollBar(BarState.Off).width('100%')

      Text(this.liuNianTitle())
        .fontSize(AppDimens.FONT_SUBTITLE).fontWeight(FontWeight.Bold)
        .fontColor(AppColors.PRIMARY).width('100%')
        .margin({ top: AppDimens.PADDING_MD, bottom: AppDimens.PADDING_SM })
      // 10 年流年: 等分宽度展开, 无需横向滚动
      Row() {
        ForEach(this.curLiuNianList(), (item: LiuNianItem) => {
          Column() {
            Text(`${item.year}`).fontSize(9).fontColor(AppColors.TEXT_SECONDARY)
            Text(`${item.age}岁`).fontSize(9).fontColor(AppColors.LUNAR_TEXT)
            Column() {
              Text(item.ganZhi).fontSize(14).fontWeight(FontWeight.Bold).fontColor(AppColors.ON_SURFACE)
              Text(item.ganShiShen).fontSize(9).fontColor(AppColors.ACCENT)
            }.margin({ top: 2 }).alignItems(HorizontalAlign.Center)
            Column() {
              Text(item.xiaoYun).fontSize(11).fontColor(AppColors.TEXT_SECONDARY)
              Text(item.xiaoYunShiShen).fontSize(9).fontColor(AppColors.LUNAR_TEXT)
            }.margin({ top: 1 }).alignItems(HorizontalAlign.Center)
          }
          .layoutWeight(1).padding({ top: 4, bottom: 4 }).margin({ right: 2 })
          .alignItems(HorizontalAlign.Center)
          .backgroundColor((this.isSelectedCurDaYun() && item.year === this.curLiuNianYear())
            ? this.redBg : AppColors.BACKGROUND)
          .borderRadius(AppDimens.BORDER_RADIUS_SM)
        })
      }
      .width('100%').alignItems(VerticalAlign.Top)
    }
    .width('100%').padding(AppDimens.PADDING_LG)
    .margin({ left: AppDimens.PADDING_MD, right: AppDimens.PADDING_MD, top: AppDimens.PADDING_MD })
    .backgroundColor(AppColors.SURFACE).borderRadius(AppDimens.BORDER_RADIUS_LG)
    .shadow({ radius: 8, color: AppColors.CARD_SHADOW_COLOR, offsetY: 2 })
  }

  liuNianTitle(): string {
    const dy = this.arg.result.daYunColumns
    if (dy === undefined || dy.length === 0) return '流年 / 小运'
    const c: BaziColumn = dy[this.curDaYunIdx()]
    const tag: string = c.startYear === this.curDaYunStart() ? ' · 当前大运' : ''
    return `流年 / 小运 (${c.gan}${c.zhi} ${c.startYear}-${c.startYear + 9}${tag})`
  }
}
