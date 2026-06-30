import SwiftUI
import UIKit

// ════════════════════════════════════════════════════════════════
//  BaziResultView — 与鸿蒙 components/BaziResultView.ets 对齐
//  顶部切换: 简洁版(宣纸命书) / 专业版(命盘表)
// ════════════════════════════════════════════════════════════════

/// 宣纸命书专用毛笔字体 (文悦古典明朝体, 与鸿蒙 WenYue.otf 同源).
/// 字体已在 SxwnlCalendarApp.init() 中通过 CTFontManagerRegisterFontsForURL 注册.
/// PostScript Name = "WenYue-GuDianMingChaoTi-NC-W5"; 未加载时 SwiftUI 自动回落系统字体.
private func wenYue(_ size: CGFloat) -> Font {
    Font.custom("WenYue-GuDianMingChaoTi-NC-W5", size: size)
}

struct BaziResultArg {
    let result: BaziResult
    let birthYear: Int
    let astEnabled: Bool
    let longitude: Double
    let latitude: Double
    let lifaLabel: String
}

struct BaziResultView: View {
    let arg: BaziResultArg

    @State private var proMode: Bool = false
    @State private var saveToast: String? = nil
    @Environment(\.dismiss) private var dismiss

    private let elemColors: [Color] = [
        Color(hex: 0x2E7D32),   // 木
        Color(hex: 0xD32F2F),   // 火
        Color(hex: 0x8D6E40),   // 土
        Color(hex: 0xC9A227),   // 金
        Color(hex: 0x1565C0)    // 水
    ]
    // 宣纸命书色板 (与鸿蒙 BaziResultView.ets 对齐)
    private let inkColor = Color(hex: 0x3A2A1C)
    private let inkSoft  = Color(hex: 0x6B5640)
    private let redColor = Color(hex: 0xC0392B)
    private let redBg    = Color(hex: 0xF6DAD2)
    private let goldLine = Color(hex: 0xC9A86A)
    // 兜底纸色: 图片加载前/失败时使用
    private let paperBg  = Color(hex: 0xFAF3E0)

    // 12 生肖水印资源名 (0=子鼠..11=亥猪), 与 Assets.xcassets / 鸿蒙端对齐
    private let zodiacNames: [String] = [
        "bz_rat",    "bz_cow",    "bz_tiger",
        "bz_rabbit", "bz_dragon", "bz_snake",
        "bz_horse",  "bz_goat",   "bz_monkey",
        "bz_hen",    "bz_dog",    "bz_pig"
    ]

    /// 根据年支返回生肖图片名
    private func zodiacName() -> String {
        let zhi = arg.result.columns.first?.zhi ?? "子"
        let order = "子丑寅卯辰巳午未申酉戌亥"
        let idx = order.firstIndex(of: zhi.first ?? "子").map {
            order.distance(from: order.startIndex, to: $0)
        } ?? 0
        return zodiacNames[max(0, min(11, idx))]
    }

    private let colTitles = ["年柱", "月柱", "日柱", "时柱", "大运", "流年"]
    private let labelW: CGFloat = 40
    private let colW: CGFloat = 54

    var body: some View {
        NavigationStack {
            ScrollView {
                VStack(spacing: 0) {
                    versionToggle.padding(.top, AppDimens.paddingMd)
                    if proMode {
                        infoBar
                        paiPanTable
                        summaryBar
                        liuNianCard
                    } else {
                        simpleCard
                    }
                }
                .padding(.bottom, 24)
            }
            .background(AppColors.background)
            .navigationTitle(isFemale ? "排盘 · 坤造" : "排盘 · 乾造")
            .navigationBarTitleDisplayMode(.inline)
            .toolbar {
                ToolbarItem(placement: .navigationBarLeading) {
                    Button("关闭") { dismiss() }
                }
            }
            .overlay(alignment: .bottom) {
                if let toast = saveToast {
                    Text(toast)
                        .font(.system(size: AppDimens.fontCaption))
                        .foregroundColor(.white)
                        .padding(.horizontal, 14).padding(.vertical, 8)
                        .background(Color.black.opacity(0.7))
                        .cornerRadius(8)
                        .padding(.bottom, 40)
                        .transition(.opacity)
                }
            }
        }
    }

    private var isFemale: Bool {
        arg.result.gender.contains("女")
    }

    // ───────────────────── 版本切换 ─────────────────────
    private var versionToggle: some View {
        HStack(spacing: 0) {
            Button { proMode = false } label: {
                Text("简洁版")
                    .font(.system(size: AppDimens.fontCaption,
                                  weight: !proMode ? .bold : .regular))
                    .foregroundColor(!proMode ? AppColors.onPrimary : AppColors.textSecondary)
                    .frame(width: 90, height: 34)
                    .background(!proMode ? AppColors.primary : AppColors.surface)
                    .clipShape(.rect(topLeadingRadius: 17, bottomLeadingRadius: 17))
            }
            Button { proMode = true } label: {
                Text("专业版")
                    .font(.system(size: AppDimens.fontCaption,
                                  weight: proMode ? .bold : .regular))
                    .foregroundColor(proMode ? AppColors.onPrimary : AppColors.textSecondary)
                    .frame(width: 90, height: 34)
                    .background(proMode ? AppColors.primary : AppColors.surface)
                    .clipShape(.rect(bottomTrailingRadius: 17, topTrailingRadius: 17))
            }
        }
        .overlay(RoundedRectangle(cornerRadius: 17)
            .stroke(AppColors.divider, lineWidth: 1))
        .clipShape(RoundedRectangle(cornerRadius: 17))
    }

    // ════════════════════ 简洁版(宣纸命书) ════════════════════
    private var simpleCard: some View {
        VStack(spacing: 0) {
            simpleReport
                .id("baziSimpleCard")

            Button { saveImage() } label: {
                Text("保存为图片")
                    .font(.system(size: AppDimens.fontBody, weight: .medium))
                    .foregroundColor(AppColors.primaryDark)
                    .frame(maxWidth: .infinity)
                    .frame(height: 40)
                    .background(AppColors.accent)
                    .cornerRadius(AppDimens.radiusLg)
            }
            .padding(.horizontal, 64)
            .padding(.vertical, AppDimens.paddingSm)
        }
    }

    private var simpleReport: some View {
        VStack(alignment: .leading, spacing: 0) {
            // 抬头
            HStack {
                Spacer()
                Text("☯ 八字命书")
                    .font(.system(size: 13, weight: .bold))
                    .foregroundColor(redColor)
            }
            .padding(.bottom, 10)

            // ── 命主信息 ──
            infoPair("姓名", strip(arg.result.userName), false,
                     "性别", isFemale ? "女" : "男")
            infoPair("生肖", strip(arg.result.shengXiao), false,
                     "年龄", strip(arg.result.age))
            infoPair("出生地",
                     String(format: "东经%.2f° 北纬%.2f°",
                            arg.longitude, arg.latitude), false,
                     "时间标准",
                     arg.astEnabled ? "真太阳时" : "北京时间(120°E)")
            singlePair("公历生日", strip(arg.result.solarBirth), false)
                .padding(.bottom, 4)
            singlePair("农历生日", strip(arg.result.lunarBirth), true)
                .padding(.bottom, 8)

            // 节气/司令段落
            Text(jieQiParagraph())
                .font(.system(size: 12))
                .foregroundColor(inkSoft)
                .lineSpacing(4)
                .padding(.bottom, 10)

            Rectangle()
                .fill(goldLine.opacity(0.6))
                .frame(height: 1)
                .padding(.bottom, 10)

            // 四柱
            siZhuBlock

            Rectangle()
                .fill(goldLine.opacity(0.6))
                .frame(height: 1)
                .padding(.top, 6).padding(.bottom, 10)

            // 起运/交运
            Text("\(strip(arg.result.qiYun))　\(strip(arg.result.jiaoYun))")
                .font(.system(size: 12))
                .foregroundColor(inkSoft)
                .lineSpacing(4)
                .padding(.bottom, 10)

            // 大运 / 流年
            Text("大运 / 流年")
                .font(.system(size: 13, weight: .bold))
                .foregroundColor(inkColor)
                .padding(.bottom, 6)
            daYunGrid
        }
        .padding(.horizontal, 16).padding(.vertical, 14)
        // 背景: 宣纸纹理图 (兜底纯色); 中部叠 12 生肖水印
        .background(
            ZStack {
                paperBg
                Image("bz_paper")
                    .resizable()
                    .scaledToFill()
                    .clipped()
                Image(zodiacName())
                    .resizable()
                    .scaledToFit()
                    .frame(width: 220, height: 220)
                    .opacity(0.07)
            }
        )
        .overlay(RoundedRectangle(cornerRadius: 10)
            .stroke(goldLine, lineWidth: 1))
        .cornerRadius(10)
        .padding(.horizontal, AppDimens.paddingMd)
        .padding(.top, AppDimens.paddingMd)
        .shadow(color: Color.black.opacity(0.1), radius: 4, y: 2)
    }

    private func infoPair(_ k1: String, _ v1: String, _ redV1: Bool,
                          _ k2: String, _ v2: String) -> some View {
        HStack(alignment: .top, spacing: 10) {
            HStack(alignment: .top, spacing: 2) {
                Text("『\(k1)』")
                    .font(.system(size: 13))
                    .foregroundColor(inkSoft)
                Text(v1)
                    .font(.system(size: 13))
                    .foregroundColor(redV1 ? redColor : inkColor)
                Spacer(minLength: 0)
            }
            .frame(maxWidth: .infinity, alignment: .leading)
            if !k2.isEmpty {
                HStack(alignment: .top, spacing: 2) {
                    Text("『\(k2)』")
                        .font(.system(size: 13))
                        .foregroundColor(inkSoft)
                    Text(v2)
                        .font(.system(size: 13))
                        .foregroundColor(inkColor)
                    Spacer(minLength: 0)
                }
                .frame(maxWidth: .infinity, alignment: .leading)
            }
        }
        .padding(.bottom, 4)
    }

    private func singlePair(_ k: String, _ v: String, _ red: Bool) -> some View {
        HStack(alignment: .top, spacing: 2) {
            Text("『\(k)』")
                .font(.system(size: 13))
                .foregroundColor(inkSoft)
            Text(v)
                .font(.system(size: 13))
                .foregroundColor(red ? redColor : inkColor)
            Spacer(minLength: 0)
        }
    }

    private var siZhuBlock: some View {
        HStack(alignment: .top, spacing: 0) {
            Text(isFemale ? "坤造" : "乾造")
                .font(wenYue(13))
                .foregroundColor(inkColor)
                .frame(width: 26, alignment: .leading)
            ForEach(Array(arg.result.columns.enumerated()), id: \.offset) { i, c in
                VStack(spacing: 2) {
                    Text(i == 2 ? "日元" : c.ganShiShen)
                        .font(.system(size: 11, weight: .medium))
                        .foregroundColor(i == 2 ? redColor : inkSoft)
                    Text(c.gan)
                        .font(wenYue(18))
                        .foregroundColor(i == 2 ? redColor : inkColor)
                    Text(c.zhi)
                        .font(wenYue(18))
                        .foregroundColor(inkColor)
                    VStack(spacing: 1) {
                        ForEach(Array(c.cangGan.enumerated()), id: \.offset) { _, cg in
                            HStack(spacing: 2) {
                                Text(cg.gan)
                                    .font(wenYue(11))
                                    .foregroundColor(inkColor)
                                Text(cg.shiShen)
                                    .font(.system(size: 10))
                                    .foregroundColor(inkSoft)
                            }
                        }
                    }
                    .padding(.top, 6)
                }
                .frame(maxWidth: .infinity)
            }
        }
    }

    private func jieQiParagraph() -> String {
        var s = ""
        let dq = strip(arg.result.dingQiType)
        if !dq.isEmpty { s += "依据\(dq)。" }
        let jq = strip(arg.result.jieQi).replacingOccurrences(
            of: "\n+", with: "；", options: .regularExpression)
        if !jq.isEmpty { s += "\(jq)。" }
        if !arg.result.siLing.isEmpty {
            s += "命主月令司令：\(arg.result.siLing)"
        }
        return s
    }

    private func daYunColumns8() -> [BaziColumn] {
        let all = arg.result.daYunColumns
        return all.count > 8 ? Array(all.prefix(8)) : all
    }

    private var daYunGrid: some View {
        HStack(alignment: .top, spacing: 0) {
            ForEach(Array(daYunColumns8().enumerated()), id: \.offset) { _, c in
                daYunColumnView(c)
            }
        }
    }

    private func daYunColumnView(_ c: BaziColumn) -> some View {
        let isCur = c.startYear == curDaYunStart
        return VStack(spacing: 0) {
            Text("\(daYunStartAge(c.startYear))")
                .font(.system(size: 11, weight: .medium))
                .foregroundColor(isCur ? redColor : inkSoft)
            VStack(spacing: 0) {
                Text(c.gan).font(wenYue(17))
                    .foregroundColor(isCur ? redColor : inkColor)
                Text(c.zhi).font(wenYue(17))
                    .foregroundColor(isCur ? redColor : inkColor)
            }
            .padding(.horizontal, 2).padding(.vertical, 1)
            .background(isCur ? redBg : .clear)
            .cornerRadius(3)
            Text("\(c.startYear)")
                .font(.system(size: 9))
                .foregroundColor(redColor)
                .padding(.top, 2)
            Rectangle().fill(goldLine.opacity(0.5))
                .frame(height: 0.5).padding(.vertical, 3)
            VStack(spacing: 0) {
                ForEach(decadeYears(c.startYear), id: \.self) { y in
                    let curLN = y == curLiuNianYear
                    Text(BaziCalc.GAN[BaziCalc.lnGan(y)]
                         + BaziCalc.ZHI[BaziCalc.lnZhi(y)])
                        .font(wenYue(13))
                        .foregroundColor(curLN ? redColor : inkColor)
                        .padding(.horizontal, 1)
                        .background(curLN ? redBg : .clear)
                        .cornerRadius(3)
                        .padding(.bottom, 2)
                }
            }
            Text("\(c.startYear + 9)")
                .font(.system(size: 9))
                .foregroundColor(redColor)
                .padding(.top, 2)
        }
        .frame(maxWidth: .infinity)
        .padding(.horizontal, 1)
    }

    private var curDaYunStart: Int { arg.result.currentDaYun?.startYear ?? -1 }
    private var curLiuNianYear: Int { arg.result.currentLiuNian?.startYear ?? -1 }
    private func daYunStartAge(_ y: Int) -> Int { y - arg.birthYear + 1 }
    private func decadeYears(_ s: Int) -> [Int] { (0..<10).map { s + $0 } }

    // ════════════════════ 专业版 ════════════════════
    private var infoBar: some View {
        VStack(alignment: .leading, spacing: 2) {
            Text("\(isFemale ? "坤" : "乾")  \(strip(arg.result.userName))")
                .font(.system(size: AppDimens.fontSubtitle, weight: .bold))
                .foregroundColor(AppColors.primary)
            Text("\(strip(arg.result.shengXiao))  \(strip(arg.result.age))")
                .font(.system(size: AppDimens.fontCaption))
                .foregroundColor(AppColors.textSecondary)
            Text("公历 \(strip(arg.result.solarBirth))")
                .font(.system(size: AppDimens.fontCaption))
                .foregroundColor(AppColors.textSecondary)
            Text("农历 \(strip(arg.result.lunarBirth))")
                .font(.system(size: AppDimens.fontCaption))
                .foregroundColor(AppColors.textSecondary)
        }
        .frame(maxWidth: .infinity, alignment: .leading)
        .padding(AppDimens.paddingMd)
        .background(AppColors.surface)
        .cornerRadius(AppDimens.radiusMd)
        .padding(.horizontal, AppDimens.paddingMd)
        .padding(.top, AppDimens.paddingMd)
    }

    private var paiPanTable: some View {
        ScrollView(.horizontal, showsIndicators: false) {
            VStack(alignment: .leading, spacing: 0) {
                headerRow
                starRow
                ganRow
                zhiRow
                cangGanRow
                textRow("纳音", rowValues(0))
                textRow("星运", rowValues(1))
                textRow("自坐", rowValues(2))
                textRow("空亡", rowValues(3))
                shenShaRow
            }
        }
        .padding(AppDimens.paddingSm)
        .background(AppColors.surface)
        .cornerRadius(AppDimens.radiusMd)
        .padding(.horizontal, AppDimens.paddingMd)
        .padding(.top, AppDimens.paddingMd)
    }

    private var headerRow: some View {
        HStack(spacing: 0) {
            Text("").frame(width: labelW)
            ForEach(Array(colTitles.enumerated()), id: \.offset) { i, t in
                Text(t)
                    .font(.system(size: AppDimens.fontCaption, weight: .bold))
                    .foregroundColor(i >= 4 ? AppColors.accent : AppColors.primary)
                    .frame(width: colW, alignment: .center)
            }
        }
        .frame(height: 28)
        .overlay(Rectangle()
            .frame(height: 1)
            .foregroundColor(AppColors.divider),
                 alignment: .bottom)
    }

    private var starRow: some View {
        HStack(spacing: 0) {
            labelCell("主星")
            ForEach(Array(allColumns.enumerated()), id: \.offset) { i, c in
                Text(c.ganShiShen)
                    .font(.system(size: AppDimens.fontCaption,
                                  weight: i == 2 ? .bold : .regular))
                    .foregroundColor(i == 2 ? AppColors.primary : AppColors.textSecondary)
                    .frame(width: colW, alignment: .center)
            }
        }
        .frame(height: 26)
    }

    private var ganRow: some View {
        HStack(spacing: 0) {
            labelCell("天干")
            ForEach(Array(allColumns.enumerated()), id: \.offset) { _, c in
                Text(c.gan)
                    .font(.system(size: 26, weight: .bold))
                    .foregroundColor(charColor(c.gan))
                    .frame(width: colW, alignment: .center)
            }
        }
        .frame(height: 44)
    }

    private var zhiRow: some View {
        HStack(spacing: 0) {
            labelCell("地支")
            ForEach(Array(allColumns.enumerated()), id: \.offset) { _, c in
                Text(c.zhi)
                    .font(.system(size: 26, weight: .bold))
                    .foregroundColor(charColor(c.zhi))
                    .frame(width: colW, alignment: .center)
            }
        }
        .frame(height: 44)
    }

    private var cangGanRow: some View {
        HStack(alignment: .top, spacing: 0) {
            labelCell("藏干").frame(maxHeight: .infinity, alignment: .top)
            ForEach(Array(allColumns.enumerated()), id: \.offset) { _, c in
                VStack(spacing: 0) {
                    ForEach(Array(c.cangGan.enumerated()), id: \.offset) { _, cg in
                        HStack(spacing: 2) {
                            Text(cg.gan)
                                .font(.system(size: AppDimens.fontSmall))
                                .foregroundColor(charColor(cg.gan))
                            Text(cg.shiShen)
                                .font(.system(size: AppDimens.fontSmall))
                                .foregroundColor(AppColors.textSecondary)
                        }
                    }
                }
                .frame(width: colW)
            }
        }
        .padding(.vertical, 4)
    }

    private func textRow(_ label: String, _ values: [String]) -> some View {
        HStack(spacing: 0) {
            labelCell(label)
            ForEach(Array(values.enumerated()), id: \.offset) { _, v in
                Text(v)
                    .font(.system(size: AppDimens.fontSmall))
                    .foregroundColor(AppColors.onSurface)
                    .frame(width: colW, alignment: .center)
            }
        }
        .frame(height: 24)
        .overlay(Rectangle()
            .frame(height: 0.5)
            .foregroundColor(AppColors.divider),
                 alignment: .top)
    }

    private var shenShaRow: some View {
        HStack(alignment: .top, spacing: 0) {
            labelCell("神煞").frame(maxHeight: .infinity, alignment: .top)
            ForEach(Array(allColumns.enumerated()), id: \.offset) { _, c in
                VStack(spacing: 1) {
                    ForEach(Array(c.shenSha.enumerated()), id: \.offset) { _, s in
                        Text(s)
                            .font(.system(size: 9))
                            .foregroundColor(AppColors.textSecondary)
                    }
                }
                .frame(width: colW)
            }
        }
        .padding(.vertical, 4)
        .overlay(Rectangle()
            .frame(height: 0.5)
            .foregroundColor(AppColors.divider),
                 alignment: .top)
    }

    private func labelCell(_ s: String) -> some View {
        Text(s)
            .font(.system(size: AppDimens.fontSmall))
            .foregroundColor(AppColors.textSecondary)
            .frame(width: labelW, alignment: .center)
    }

    private var summaryBar: some View {
        VStack(spacing: AppDimens.paddingSm) {
            HStack {
                Text("司令")
                    .font(.system(size: AppDimens.fontCaption))
                    .foregroundColor(AppColors.textSecondary)
                Text(arg.result.siLing)
                    .font(.system(size: AppDimens.fontBody, weight: .bold))
                    .foregroundColor(AppColors.primary)
                Spacer()
                Text(wuXingSummary())
                    .font(.system(size: AppDimens.fontCaption))
                    .foregroundColor(AppColors.onSurface)
            }
            HStack(spacing: 4) {
                let chips = statusChips()
                ForEach(Array(chips.enumerated()), id: \.offset) { i, c in
                    Text(c)
                        .font(.system(size: AppDimens.fontCaption))
                        .foregroundColor(AppColors.onPrimary)
                        .frame(maxWidth: .infinity)
                        .frame(height: 28)
                        .background(i < elemColors.count ? elemColors[i] : AppColors.primary)
                        .cornerRadius(AppDimens.radiusSm)
                }
            }
        }
        .padding(AppDimens.paddingMd)
        .background(AppColors.surface)
        .cornerRadius(AppDimens.radiusMd)
        .padding(.horizontal, AppDimens.paddingMd)
        .padding(.top, AppDimens.paddingMd)
    }

    private var liuNianCard: some View {
        VStack(alignment: .leading, spacing: 0) {
            Text("大运")
                .font(.system(size: AppDimens.fontSubtitle, weight: .bold))
                .foregroundColor(AppColors.primary)
                .padding(.bottom, AppDimens.paddingSm)
            ScrollView(.horizontal, showsIndicators: false) {
                HStack(spacing: 4) {
                    ForEach(Array(arg.result.daYunColumns.enumerated()), id: \.offset) { _, c in
                        proDaYunCell(c)
                    }
                }
            }

            Text("流年 / 小运（当前大运）")
                .font(.system(size: AppDimens.fontSubtitle, weight: .bold))
                .foregroundColor(AppColors.primary)
                .padding(.top, AppDimens.paddingMd)
                .padding(.bottom, AppDimens.paddingSm)
            ScrollView(.horizontal, showsIndicators: false) {
                HStack(spacing: 4) {
                    ForEach(Array(arg.result.liuNian.enumerated()), id: \.offset) { _, it in
                        liuNianCell(it)
                    }
                }
            }
        }
        .padding(AppDimens.paddingLg)
        .frame(maxWidth: .infinity, alignment: .leading)
        .background(AppColors.surface)
        .cornerRadius(AppDimens.radiusLg)
        .padding(.horizontal, AppDimens.paddingMd)
        .padding(.top, AppDimens.paddingMd)
        .shadow(color: Color.black.opacity(0.06), radius: 4, y: 2)
    }

    private func proDaYunCell(_ c: BaziColumn) -> some View {
        let isCur = c.startYear == curDaYunStart
        return VStack(spacing: 0) {
            Text(c.ganShiShen)
                .font(.system(size: 9))
                .foregroundColor(AppColors.textSecondary)
            Text(c.gan)
                .font(.system(size: 18, weight: .bold))
                .foregroundColor(charColor(c.gan))
            Text(c.zhi)
                .font(.system(size: 18, weight: .bold))
                .foregroundColor(charColor(c.zhi))
            Text("\(c.startYear)")
                .font(.system(size: 9))
                .foregroundColor(AppColors.textSecondary)
                .padding(.top, 2)
        }
        .frame(width: 48)
        .padding(.vertical, 4)
        .background(isCur ? redBg : AppColors.background)
        .cornerRadius(AppDimens.radiusSm)
    }

    private func liuNianCell(_ it: LiuNianItem) -> some View {
        let isCur = it.year == curLiuNianYear
        return VStack(spacing: 0) {
            Text("\(it.year)")
                .font(.system(size: 9))
                .foregroundColor(AppColors.textSecondary)
            Text("\(it.age)岁")
                .font(.system(size: 9))
                .foregroundColor(AppColors.lunarText)
            HStack(spacing: 1) {
                Text(it.ganZhi)
                    .font(.system(size: 13, weight: .bold))
                    .foregroundColor(AppColors.onSurface)
                Text(it.ganShiShen)
                    .font(.system(size: 9))
                    .foregroundColor(AppColors.accent)
            }
            .padding(.top, 2)
            HStack(spacing: 1) {
                Text(it.xiaoYun)
                    .font(.system(size: 11))
                    .foregroundColor(AppColors.textSecondary)
                Text(it.xiaoYunShiShen)
                    .font(.system(size: 9))
                    .foregroundColor(AppColors.lunarText)
            }
        }
        .frame(width: 54)
        .padding(.vertical, 4)
        .background(isCur ? redBg : AppColors.background)
        .cornerRadius(AppDimens.radiusSm)
    }

    // ── 工具 ────────────────────────────────────────────────
    private var allColumns: [BaziColumn] {
        var list = arg.result.columns
        if let d = arg.result.currentDaYun { list.append(d) }
        if let l = arg.result.currentLiuNian { list.append(l) }
        return list
    }

    private func rowValues(_ field: Int) -> [String] {
        allColumns.map { c in
            switch field {
            case 0: return c.nayin
            case 1: return c.xingYun
            case 2: return c.ziZuo
            default: return c.kongWang
            }
        }
    }

    private func strip(_ s: String) -> String {
        if let i = s.firstIndex(of: "』") {
            return String(s[s.index(after: i)...])
                .trimmingCharacters(in: .whitespaces)
        }
        return s.trimmingCharacters(in: .whitespaces)
    }

    private func charColor(_ ch: String) -> Color {
        let ganRange = "甲乙丙丁戊己庚辛壬癸"
        if let i = ganRange.firstIndex(of: Character(ch)) {
            return elemColors[ganRange.distance(from: ganRange.startIndex, to: i) / 2]
        }
        let zhiRange = "子丑寅卯辰巳午未申酉戌亥"
        if let i = zhiRange.firstIndex(of: Character(ch)) {
            let zw = [4, 2, 0, 0, 2, 1, 1, 2, 3, 3, 2, 4]
            return elemColors[zw[zhiRange.distance(from: zhiRange.startIndex, to: i)]]
        }
        return AppColors.onSurface
    }

    private func wuXingSummary() -> String {
        let names = ["木", "火", "土", "金", "水"]
        var s = "含藏气: "
        for i in 0..<5 {
            if i < arg.result.wuXingCount.count {
                s += "\(arg.result.wuXingCount[i])\(names[i]) "
            }
        }
        return s.trimmingCharacters(in: .whitespaces)
    }

    private func statusChips() -> [String] {
        let names = ["木", "火", "土", "金", "水"]
        return (0..<5).map { i in
            let stat = i < arg.result.wuXingStatus.count
                ? arg.result.wuXingStatus[i] : ""
            return "\(names[i])\(stat)"
        }
    }

    // ── 截图保存 ────────────────────────────────────────────
    private func saveImage() {
        let view = simpleReport
            .frame(width: UIScreen.main.bounds.width - 2 * AppDimens.paddingMd)
        let renderer = ImageRenderer(content: view)
        renderer.scale = UIScreen.main.scale
        guard let ui = renderer.uiImage else {
            showToast("截图失败"); return
        }
        UIImageWriteToSavedPhotosAlbum(ui, nil, nil, nil)
        showToast("已保存到相册")
    }

    private func showToast(_ msg: String) {
        withAnimation { saveToast = msg }
        DispatchQueue.main.asyncAfter(deadline: .now() + 2.5) {
            withAnimation { saveToast = nil }
        }
    }
}
