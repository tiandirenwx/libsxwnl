import Foundation

// ════════════════════════════════════════════════════════════════
//  BaziCalc — 与鸿蒙 model/BaziInput.ets 对齐的纯计算工具
//
//  仅依赖 SxwnlBridge 提供的数据接口, 不引入任何 UI。
// ════════════════════════════════════════════════════════════════

enum BaziCalc {

    static let GAN: [String] = ["甲", "乙", "丙", "丁", "戊",
                                "己", "庚", "辛", "壬", "癸"]
    static let ZHI: [String] = ["子", "丑", "寅", "卯", "辰", "巳",
                                "午", "未", "申", "酉", "戌", "亥"]

    /// 六十甲子序号 → 名称
    static func jiaZiName(_ idx: Int) -> String {
        let i = ((idx % 60) + 60) % 60
        return GAN[i % 10] + ZHI[i % 12]
    }

    /// 干支索引 → 六十甲子序号
    static func gz60(gan: Int, zhi: Int) -> Int {
        for n in 0..<60 where n % 10 == gan && n % 12 == zhi {
            return n
        }
        return 0
    }

    /// 五虎遁: 由年干得 12 个合法月柱(寅月起)的六十甲子序号
    static func validMonths(yearIdx: Int) -> [Int] {
        let yg = ((yearIdx % 60) + 60) % 60 % 10
        let base = ((yg % 5) * 2 + 2) % 10
        return (0..<12).map { k in
            gz60(gan: (base + k) % 10, zhi: (2 + k) % 12)
        }
    }

    /// 五鼠遁: 由日干得 12 个合法时柱(子时起)的六十甲子序号
    static func validHours(dayIdx: Int) -> [Int] {
        let dg = ((dayIdx % 60) + 60) % 60 % 10
        let base = ((dg % 5) * 2) % 10
        return (0..<12).map { h in
            gz60(gan: (base + h) % 10, zhi: h)
        }
    }

    /// 在合法列表中按地支保留, 返回新的六十甲子序号
    static func remapByZhi(validList: [Int], oldIdx: Int) -> Int {
        let z = ((oldIdx % 60) + 60) % 60 % 12
        for v in validList where v % 12 == z { return v }
        return validList.first ?? 0
    }

    // ── 流年(立春界)干支 ──
    static func lnGan(_ year: Int) -> Int { (((year - 4) % 10) + 10) % 10 }
    static func lnZhi(_ year: Int) -> Int { (((year - 4) % 12) + 12) % 12 }

    // ── 年份解析 (B212 / *212 / -211) ──
    static func parseYear(_ s: String) -> Int? {
        let astro = YearUtil.yearStrToAstro(s)
        return astro == Int.min ? nil : astro
    }

    /// 天文年 → 显示文本: 2024 → "2024年", -210 → "公元前211年"
    static func formatYear(_ y: Int) -> String {
        if y <= 0 { return "公元前\(1 - y)年" }
        return "\(y)年"
    }

    // ── 月长 ──
    static func solarDaysInMonth(year: Int, month: Int) -> Int {
        if month == 2 {
            let leap = (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0)
            return leap ? 29 : 28
        }
        if month == 4 || month == 6 || month == 9 || month == 11 { return 30 }
        return 31
    }

    static func lunarDaysInMonth(year: Int, month: Int,
                                 isLeap: Bool, isSpec: Bool) -> Int {
        let d = SxwnlBridge.getLunarMonthDays(year: year, month: month,
                                              isLeap: isLeap, isSpec: isSpec)
        return (d == 29 || d == 30) ? d : 30
    }

    static func pad2(_ n: Int) -> String {
        n < 10 ? "0\(n)" : "\(n)"
    }

    /// 出生时间一条记录文本: "公历 2026年06月02日 13点16分"
    static func formatRecord(_ sel: DateSelection,
                             lunarMonths: [LunarMonth]) -> String {
        let ymd = "\(formatYear(sel.year))\(pad2(sel.month))月\(pad2(sel.day))日"
        let hm = "\(sel.hour)点\(pad2(sel.minute))分"
        if sel.inputMode == .lunar {
            var mName = "\(sel.month)月"
            for m in lunarMonths where m.month == sel.month && m.isLeap == sel.isLeap {
                mName = m.name
                break
            }
            return "农历 \(formatYear(sel.year))\(mName)\(sel.day)日 \(hm)"
        }
        return "公历 \(ymd) \(hm)"
    }

    /// 表单 → SxwnlBridge.calcBazi 的入参 (调用时直接展开传入)
    struct BaziCallParams {
        let name: String
        let gender: Bool
        let year: Int
        let month: Int
        let day: Int
        let hour: Int
        let minute: Int
        let isLunar: Bool
        let isLeap: Bool
        let isSpec: Bool
        let astEnabled: Bool
        let longitude: Double
        let latitude: Double
        let lifa: Int
    }

    static func toParams(_ form: BaziFormData) -> BaziCallParams {
        let s = form.selection
        return BaziCallParams(
            name: form.name.isEmpty ? "匿名" : form.name,
            gender: form.gender,
            year: s.year, month: s.month, day: s.day,
            hour: s.hour, minute: s.minute,
            isLunar: s.inputMode == .lunar,
            isLeap: s.inputMode == .lunar ? s.isLeap : false,
            isSpec: s.inputMode == .lunar ? s.isSpec : false,
            astEnabled: form.astEnabled,
            longitude: form.longitude,
            latitude: form.latitude,
            lifa: form.lifaModern.rawValue
        )
    }
}
