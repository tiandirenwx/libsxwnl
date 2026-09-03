import Foundation

// ════════════════════════════════════════════════════════════════
//  SxwnlBridge.swift — 与 capi/sxwnl_capi.h 一一对应的 Swift 封装
//  对齐鸿蒙端 NativeBridge.ets 的全部接口与字段
// ════════════════════════════════════════════════════════════════

// MARK: - 单日完整信息 (含节日 / 回历 / 纳音 / JD 扩展字段)

struct DayInfo {
    let solarYear: Int
    let solarMonth: Int
    let solarDay: Int

    let lunarYear: Int
    let lunarMonth: Int
    let lunarDay: Int
    let isLeapMonth: Bool

    let weekDay: Int

    let yearGan: Int
    let yearZhi: Int
    let monthGan: Int
    let monthZhi: Int
    let dayGan: Int
    let dayZhi: Int

    // 与鸿蒙 NAPI/Android 一致, 数值索引名为 jieQi/yueXiang;
    // 与字符串 jieQiName/yueXiangName 共存, 不冲突.
    let jieQi: Int        // -1 表示无, 0..23
    let yueXiang: Int     // -1 表示无, 0..3
    let constellation: Int
    let jian12: Int

    let yearGZ: String
    let monthGZ: String
    let dayGZ: String
    let lunarMonthName: String
    let lunarDayName: String
    let jieQiName: String
    let jieQiTime: String
    // 历谱口径节气(整日表+QB, 对齐 sxwnl 网页版): 日历格子标签建议用此值,
    // 古代(1645年前)与天文口径 jieQi 可能差 1 天; 精确时刻仍用 jieQiTime.
    let lipuJieQi: Int    // -1 表示无, 0..23
    let lipuJieQiName: String
    let shengXiao: String
    let constellationName: String
    let weekName: String
    let yueXiangName: String
    let jian12Name: String

    // ── 节日 & 杂节 ─────────────────────────
    let holiday: String   // A 类
    let major: String     // B 类
    let minor: String     // C 类
    let misc: String      // 杂节(数九/三伏...)
    let isOffDay: Bool

    // ── 纪年扩展 ───────────────────────────
    let lunarJunYear: Int
    let lunarLichunYear: Int
    let huangdiYear: Int

    // ── 回历 ───────────────────────────────
    let moslemYear: Int
    let moslemMonth: Int
    let moslemDay: Int

    // ── 纳音 ───────────────────────────────
    let yearNaYin: String
    let monthNaYin: String
    let dayNaYin: String

    // ── 月相时刻 ───────────────────────────
    let yueXiangTime: String

    // ── 儒略日 (对应当日 12:00) ─────────────
    let julianDay: Double
}

// MARK: - 节气 & 农历月

struct JieQiItem {
    let idx: Int
    let solarMonth: Int
    let solarDay: Int
    let name: String
    let time: String
}

struct LunarMonth {
    let month: Int       // 1..12
    let isLeap: Bool
    let isSpec: Bool     // 古代特殊月
    let name: String     // "正月"/"闰四月"/"后九月"
}

struct SolarDate {
    let year: Int
    let month: Int
    let day: Int
}

// MARK: - 年历 (按农历月聚合)

struct YearCalJieQi {
    let idx: Int            // 0..23
    let name: String
    let gz: String
    let solarMonth: Int
    let solarDay: Int
    let time: String        // "HH:MM:SS" 精确交气时刻(天文定气)
    let dayOffset: Int
    let dayName: String     // 月内日名 "初一"/"十五"...
    let accMonth: Int       // 精确交气公历月 (古代可能与 solarMonth 差1天)
    let accDay: Int         // 精确交气公历日
}

struct YearCalMonth {
    let monthIdx: Int
    let monthName: String   // "正月"/"闰二月"/"后九月"
    let isLeap: Bool
    let isSpec: Bool
    let dayCount: Int       // 29/30
    let solarYear: Int
    let solarMonth: Int
    let solarDay: Int
    let shuoGZ: String
    let jieQi: [YearCalJieQi]
}

// MARK: - 日月升降 (RTS)

struct DayRTS {
    let sunRise: String
    let sunSet: String
    let sunMeridian: String
    let moonRise: String
    let moonSet: String
    let moonMeridian: String
    let civilDawn: String
    let civilDusk: String
    let dayLength: String
    let lightLength: String
}

// MARK: - 八字

struct CangGanItem {
    let gan: String
    let shiShen: String
}

struct BaziColumn {
    let gan: String
    let zhi: String
    let ganShiShen: String   // 主星 (日柱为 "男"/"女")
    let nayin: String
    let xingYun: String
    let ziZuo: String
    let kongWang: String
    let cangGan: [CangGanItem]
    let shenSha: [String]
    let startYear: Int       // 大运/流年起年; 四柱为 0
}

struct LiuNianItem {
    let year: Int
    let age: Int
    let ganZhi: String
    let ganShiShen: String
    let zhiShiShen: String
    let xiaoYun: String
    let xiaoYunShiShen: String
}

struct ReverseItem {
    let year: Int
    let month: Int
    let day: Int
    let hour: Int            // -1 表示时柱未知
    let yearGZ: String
    let monthGZ: String
    let dayGZ: String
    let hourGZ: String
}

struct BaziResult {
    let userName: String
    let gender: String
    let solarBirth: String
    let lunarBirth: String
    let dateOfBirth: String
    let shengXiao: String
    let age: String
    let lifa: String
    let dingQiType: String
    let ast: String
    let jieQi: String
    let qiYun: String
    let jiaoYun: String

    // ── 完整盘面 ────────────────────────────
    let columns: [BaziColumn]            // 年/月/日/时 完整 4 列
    let currentDaYun: BaziColumn?        // 当前大运
    let currentLiuNian: BaziColumn?      // 当前流年
    let daYunColumns: [BaziColumn]       // 大运完整列
    let wuXingCount: [Int]               // [木火土金水] (含藏干)
    let wuXingStatus: [String]           // [旺/相/休/囚/死]
    let siLing: String                   // 司令
    let liuNian: [LiuNianItem]           // 流年+小运列表 (从出生年起 120 年)

    // 按大运分组的流年: 与 daYunColumns 一一对应, 每桶 10 个.
    //   对齐鸿蒙 NAPI 返回字段 (napi_bazi.cpp: result.v("liuNianAll", ...)).
    let liuNianAll: [[LiuNianItem]]
}

// MARK: - SxwnlBridge

final class SxwnlBridge {

    // ═══ Calendar / Day ═══════════════════════════════════════

    static func getDayInfo(year: Int, month: Int, day: Int) -> DayInfo? {
        var raw = SxwnlDayInfo()
        guard sxwnl_get_day_info(Int32(year), Int32(month), Int32(day), &raw) == 0 else {
            return nil
        }
        return convertDayInfo(raw)
    }

    static func getMonthDays(year: Int, month: Int) -> [DayInfo] {
        var buffer = [SxwnlDayInfo](repeating: SxwnlDayInfo(), count: 42)
        let n = sxwnl_get_month_days(Int32(year), Int32(month), &buffer, 42)
        guard n > 0 else { return [] }
        return buffer.prefix(Int(n)).map(convertDayInfo)
    }

    /// 与鸿蒙端 NativeBridge.getMonthData 一致(更直观的名字)
    static func getMonthData(year: Int, month: Int) -> [DayInfo] {
        return getMonthDays(year: year, month: month)
    }

    static func lunarToSolar(year: Int, month: Int, day: Int, isLeap: Bool) -> SolarDate? {
        var oy: Int32 = 0, om: Int32 = 0, od: Int32 = 0
        guard sxwnl_lunar_to_solar(Int32(year), Int32(month), Int32(day),
                                   isLeap, &oy, &om, &od) == 0 else { return nil }
        return SolarDate(year: Int(oy), month: Int(om), day: Int(od))
    }

    static func solarToLunar(year: Int, month: Int, day: Int) -> DayInfo? {
        var raw = SxwnlDayInfo()
        guard sxwnl_solar_to_lunar(Int32(year), Int32(month), Int32(day), &raw) == 0 else {
            return nil
        }
        return convertDayInfo(raw)
    }

    static func getYearLeapMonth(year: Int) -> Int {
        Int(sxwnl_get_year_leap_month(Int32(year)))
    }

    static func getLunarMonths(year: Int) -> [LunarMonth] {
        var buffer = [SxwnlLunarMonth](repeating: SxwnlLunarMonth(), count: 16)
        let n = sxwnl_get_lunar_months(Int32(year), &buffer, 16)
        guard n > 0 else { return [] }
        return buffer.prefix(Int(n)).map { m in
            LunarMonth(
                month: Int(m.month),
                isLeap: m.is_leap != 0,
                isSpec: m.is_spec != 0,
                name: extractString(from: m.name, capacity: 16)
            )
        }
    }

    static func getLunarMonthDays(year: Int, month: Int, isLeap: Bool, isSpec: Bool) -> Int {
        Int(sxwnl_get_lunar_month_days(Int32(year), Int32(month), isLeap, isSpec))
    }

    /// 公历某年月真实存在的日号(升序); 自动跳过 1582-10 改历缺口。
    static func getSolarMonthValidDays(year: Int, month: Int) -> [Int] {
        var buffer = [Int32](repeating: 0, count: 31)
        let n = sxwnl_get_solar_month_valid_days(Int32(year), Int32(month), &buffer, 31)
        guard n > 0 else { return [] }
        return buffer.prefix(Int(n)).map { Int($0) }
    }

    // ═══ Year Calendar ════════════════════════════════════════

    static func getYearCalendar(year: Int) -> [YearCalMonth] {
        var buffer = [SxwnlYearCalMonth](repeating: SxwnlYearCalMonth(), count: 16)
        let n = sxwnl_get_year_calendar(Int32(year), &buffer, 16)
        guard n > 0 else { return [] }
        return buffer.prefix(Int(n)).map { m in
            var jqList: [YearCalJieQi] = []
            let jqCount = Int(m.jq_count)
            var jqTuple = m.jieqi
            withUnsafePointer(to: &jqTuple) { ptr in
                ptr.withMemoryRebound(to: SxwnlYearCalJieQi.self, capacity: 4) { base in
                    for i in 0..<min(jqCount, 4) {
                        let j = base[i]
                        jqList.append(YearCalJieQi(
                            idx: Int(j.idx),
                            name: extractString(from: j.name, capacity: 12),
                            gz: extractString(from: j.gz, capacity: 8),
                            solarMonth: Int(j.solar_month),
                            solarDay: Int(j.solar_day),
                            time: extractString(from: j.time, capacity: 32),
                            dayOffset: Int(j.day_offset),
                            dayName: extractString(from: j.day_name, capacity: 12),
                            accMonth: Int(j.acc_month),
                            accDay: Int(j.acc_day)
                        ))
                    }
                }
            }
            return YearCalMonth(
                monthIdx: Int(m.month_idx),
                monthName: extractString(from: m.month_name, capacity: 16),
                isLeap: m.is_leap != 0,
                isSpec: m.is_spec != 0,
                dayCount: Int(m.day_count),
                solarYear: Int(m.solar_year),
                solarMonth: Int(m.solar_month),
                solarDay: Int(m.solar_day),
                shuoGZ: extractString(from: m.shuo_gz, capacity: 8),
                jieQi: jqList
            )
        }
    }

    // ═══ JieQi ════════════════════════════════════════════════

    static func getJieQiList(year: Int) -> [JieQiItem] {
        var buffer = [SxwnlJieQiItem](repeating: SxwnlJieQiItem(), count: 30)
        let n = sxwnl_get_jieqi_list(Int32(year), &buffer, 30)
        guard n > 0 else { return [] }
        return buffer.prefix(Int(n)).map { item in
            JieQiItem(
                idx: Int(item.idx),
                solarMonth: Int(item.solar_month),
                solarDay: Int(item.solar_day),
                name: extractString(from: item.name, capacity: 12),
                time: extractString(from: item.time, capacity: 32)
            )
        }
    }

    // ═══ RTS (日月升降/中天) ═══════════════════════════════════

    static func calcDayRTS(year: Int, month: Int, day: Int,
                           longitude: Double, latitude: Double,
                           tzHours: Double) -> DayRTS? {
        var raw = SxwnlDayRTS()
        guard sxwnl_calc_day_rts(Int32(year), Int32(month), Int32(day),
                                 longitude, latitude, tzHours, &raw) == 0 else {
            return nil
        }
        return DayRTS(
            sunRise: extractString(from: raw.sun_rise, capacity: 16),
            sunSet: extractString(from: raw.sun_set, capacity: 16),
            sunMeridian: extractString(from: raw.sun_meridian, capacity: 16),
            moonRise: extractString(from: raw.moon_rise, capacity: 16),
            moonSet: extractString(from: raw.moon_set, capacity: 16),
            moonMeridian: extractString(from: raw.moon_meridian, capacity: 16),
            civilDawn: extractString(from: raw.civil_dawn, capacity: 16),
            civilDusk: extractString(from: raw.civil_dusk, capacity: 16),
            dayLength: extractString(from: raw.day_length, capacity: 16),
            lightLength: extractString(from: raw.light_length, capacity: 16)
        )
    }

    // ═══ Bazi ════════════════════════════════════════════════

    /// 与鸿蒙端对齐的完整八字排盘 — 新接口
    static func calcBazi(
        name: String,
        gender: Bool,
        year: Int, month: Int, day: Int,
        hour: Int, minute: Int, second: Double = 0,
        isLunar: Bool,
        isLeap: Bool = false,
        isSpec: Bool = false,
        astEnabled: Bool,
        longitude: Double,
        latitude: Double,
        lifa: Int
    ) -> BaziResult? {
        var params = SxwnlBaziParams()
        writeCString(name, into: &params.name, capacity: 64)
        params.gender = gender
        params.is_ast = astEnabled
        params.longitude = longitude
        params.latitude = latitude
        params.lifa = Int32(lifa)
        params.is_lunar = isLunar
        params.is_leap = isLeap
        params.is_spec = isSpec
        params.year = Int32(year)
        params.month = Int32(month)
        params.day = Int32(day)
        params.hour = Int32(hour)
        params.minute = Int32(minute)
        params.second = second

        guard let bazi = sxwnl_bazi_create(&params) else { return nil }
        defer { sxwnl_bazi_free(bazi) }

        return buildBaziResult(bazi)
    }

    /// 四柱反推 — 时柱未知时 hg=hz=-1
    static func baziReverse(yg: Int, yz: Int, mg: Int, mz: Int,
                            dg: Int, dz: Int, hg: Int, hz: Int,
                            startYear: Int, endYear: Int) -> [ReverseItem] {
        let sz: [Int32] = [Int32(yg), Int32(yz), Int32(mg), Int32(mz),
                           Int32(dg), Int32(dz), Int32(hg), Int32(hz)]
        var buffer = [SxwnlReverseItem](repeating: SxwnlReverseItem(), count: 64)
        let n = sz.withUnsafeBufferPointer { sptr -> Int32 in
            sxwnl_bazi_reverse(sptr.baseAddress, Int32(startYear), Int32(endYear),
                               &buffer, 64)
        }
        guard n > 0 else { return [] }
        return buffer.prefix(Int(n)).map { r in
            let g: (Int) -> String = { i in
                var t = r.ganzhi
                return withUnsafePointer(to: &t) { ptr -> String in
                    ptr.withMemoryRebound(to: CChar.self, capacity: 32) { base in
                        String(cString: base.advanced(by: i * 8))
                    }
                }
            }
            return ReverseItem(
                year: Int(r.year),
                month: Int(r.month),
                day: Int(r.day),
                hour: Int(r.hour),
                yearGZ: g(0), monthGZ: g(1),
                dayGZ: g(2), hourGZ: g(3)
            )
        }
    }

    // ═══ Internal helpers ═════════════════════════════════════

    private static func buildBaziResult(_ bazi: SxwnlBazi) -> BaziResult {
        // ── 文本字段 ──
        let userName = cstr(sxwnl_bazi_get_user_name(bazi))
        let genderStr = cstr(sxwnl_bazi_get_gender(bazi))
        let solarBirth = cstr(sxwnl_bazi_get_solar_birth(bazi))
        let lunarBirth = cstr(sxwnl_bazi_get_lunar_birth(bazi))
        let dateOfBirth = cstr(sxwnl_bazi_get_date_of_birth(bazi))
        let shengXiao = cstr(sxwnl_bazi_get_sheng_xiao(bazi))
        let age = cstr(sxwnl_bazi_get_age(bazi))
        let lifaStr = cstr(sxwnl_bazi_get_lifa(bazi))
        let dingQiType = cstr(sxwnl_bazi_get_ding_qi_type(bazi))
        let ast = cstr(sxwnl_bazi_get_ast(bazi))
        let jieQi = cstr(sxwnl_bazi_get_jie_qi(bazi))
        let qiYun = cstr(sxwnl_bazi_get_qi_yun(bazi))
        let jiaoYun = cstr(sxwnl_bazi_get_jiao_yun(bazi))

        // ── 大运计数 (用于后续大运完整列) ──
        let dyCount = sxwnl_bazi_get_da_yun_count(bazi)

        // ── 完整 4 柱 ──
        var colBuf = [SxwnlBaziColumn](repeating: SxwnlBaziColumn(), count: 4)
        let colCount = sxwnl_bazi_get_columns(bazi, &colBuf)
        let columns = colBuf.prefix(Int(colCount)).map { convertColumn($0) }

        // ── 当前大运 ──
        var curDy = SxwnlBaziColumn()
        let curDyOK = sxwnl_bazi_get_current_da_yun(bazi, &curDy) == 0
        let currentDaYun: BaziColumn? = curDyOK ? convertColumn(curDy) : nil

        // ── 当前流年 ──
        var curLn = SxwnlBaziColumn()
        let curLnOK = sxwnl_bazi_get_current_liu_nian(bazi, &curLn) == 0
        let currentLiuNian: BaziColumn? = curLnOK ? convertColumn(curLn) : nil

        // ── 大运完整列 ──
        var dyColBuf = [SxwnlBaziColumn](repeating: SxwnlBaziColumn(), count: max(Int(dyCount), 1))
        let dyColN = sxwnl_bazi_get_da_yun_columns(bazi, &dyColBuf, dyCount)
        let daYunColumns = dyColBuf.prefix(Int(dyColN)).map { convertColumn($0) }

        // ── 五行 ──
        var wxBuf = [Int32](repeating: 0, count: 5)
        sxwnl_bazi_get_wuxing_count(bazi, &wxBuf, true)
        let wuXingCount = wxBuf.map { Int($0) }

        typealias C8 = (CChar, CChar, CChar, CChar, CChar, CChar, CChar, CChar)
        let emptyC8: C8 = (0, 0, 0, 0, 0, 0, 0, 0)
        var statusBuf: (
            (CChar, CChar, CChar, CChar, CChar, CChar, CChar, CChar),
            (CChar, CChar, CChar, CChar, CChar, CChar, CChar, CChar),
            (CChar, CChar, CChar, CChar, CChar, CChar, CChar, CChar),
            (CChar, CChar, CChar, CChar, CChar, CChar, CChar, CChar),
            (CChar, CChar, CChar, CChar, CChar, CChar, CChar, CChar)
        ) = (emptyC8, emptyC8, emptyC8, emptyC8, emptyC8)
        withUnsafeMutablePointer(to: &statusBuf) { ptr in
            ptr.withMemoryRebound(to: (CChar, CChar, CChar, CChar, CChar, CChar, CChar, CChar).self, capacity: 5) { base in
                sxwnl_bazi_get_wuxing_status(bazi, base)
            }
        }
        var wuXingStatus: [String] = []
        withUnsafePointer(to: &statusBuf) { ptr in
            ptr.withMemoryRebound(to: (CChar, CChar, CChar, CChar, CChar, CChar, CChar, CChar).self, capacity: 5) { base in
                for i in 0..<5 {
                    wuXingStatus.append(extractString(from: base[i], capacity: 8))
                }
            }
        }

        // ── 司令 ──
        let siLing = cstr(sxwnl_bazi_get_si_ling(bazi))

        // ── 流年 + 小运 列表 (默认起始为出生年, 取 120 年) ──
        let lnStart: Int32 = parseYear(fromBirth: solarBirth) ?? 1900
        let liuNian = fetchLiuNian(bazi, startYear: lnStart, count: 120)

        // ── 按大运分组的流年: 与 daYunColumns 一一对应, 每桶 10 个 ──
        //   对齐鸿蒙 NAPI 行为. 出错或 startYear 无效的桶用空列表占位,
        //   保证下游可直接用 daYunColumns[i] ↔ liuNianAll[i] 一一访问.
        let liuNianAll: [[LiuNianItem]] = dyColBuf.prefix(Int(dyColN)).map { col in
            col.start_year > 0
                ? fetchLiuNian(bazi, startYear: col.start_year, count: 10)
                : []
        }

        return BaziResult(
            userName: userName, gender: genderStr,
            solarBirth: solarBirth, lunarBirth: lunarBirth,
            dateOfBirth: dateOfBirth, shengXiao: shengXiao,
            age: age, lifa: lifaStr, dingQiType: dingQiType,
            ast: ast, jieQi: jieQi, qiYun: qiYun, jiaoYun: jiaoYun,
            columns: columns,
            currentDaYun: currentDaYun,
            currentLiuNian: currentLiuNian,
            daYunColumns: daYunColumns,
            wuXingCount: wuXingCount,
            wuXingStatus: wuXingStatus,
            siLing: siLing,
            liuNian: liuNian,
            liuNianAll: liuNianAll
        )
    }

    /// 读取 [startYear, startYear+count) 区间的流年+小运. 越界/错误均返回空列表.
    private static func fetchLiuNian(
        _ bazi: SxwnlBazi, startYear: Int32, count: Int32
    ) -> [LiuNianItem] {
        let cap = max(Int(count), 1)
        var buf = [SxwnlLiuNianItem](repeating: SxwnlLiuNianItem(), count: cap)
        let n = sxwnl_bazi_get_liu_nian(bazi, startYear, &buf, count)
        guard n > 0 else { return [] }
        return buf.prefix(Int(n)).map { l in
            LiuNianItem(
                year: Int(l.year),
                age: Int(l.age),
                ganZhi: extractString(from: l.gan_zhi, capacity: 8),
                ganShiShen: extractString(from: l.gan_shi_shen, capacity: 8),
                zhiShiShen: extractString(from: l.zhi_shi_shen, capacity: 8),
                xiaoYun: extractString(from: l.xiao_yun, capacity: 8),
                xiaoYunShiShen: extractString(from: l.xiao_yun_shi_shen, capacity: 8)
            )
        }
    }

    /// "公历 2024年01月01日 ..." → 2024
    private static func parseYear(fromBirth s: String) -> Int32? {
        for tok in s.split(separator: " ") {
            let str = String(tok)
            if str.contains("年"),
               let y = Int(str.split(separator: "年").first ?? "") {
                return Int32(y)
            }
        }
        return nil
    }

    private static func convertColumn(_ c: SxwnlBaziColumn) -> BaziColumn {
        // 藏干
        var cgList: [CangGanItem] = []
        var cgTuple = c.cang_gan
        var ssTuple = c.cang_gan_shi_shen
        let cnt = Int(c.cang_gan_count)
        withUnsafePointer(to: &cgTuple) { gp in
            gp.withMemoryRebound(to: (CChar, CChar, CChar, CChar, CChar, CChar, CChar, CChar).self,
                                 capacity: 3) { gbase in
                withUnsafePointer(to: &ssTuple) { sp in
                    sp.withMemoryRebound(to: (CChar, CChar, CChar, CChar, CChar, CChar, CChar, CChar).self,
                                         capacity: 3) { sbase in
                        for i in 0..<min(cnt, 3) {
                            cgList.append(CangGanItem(
                                gan: extractString(from: gbase[i], capacity: 8),
                                shiShen: extractString(from: sbase[i], capacity: 8)
                            ))
                        }
                    }
                }
            }
        }
        // 神煞
        var ssList: [String] = []
        let ssCnt = Int(c.shen_sha_count)
        var ssArr = c.shen_sha
        withUnsafePointer(to: &ssArr) { ptr in
            ptr.withMemoryRebound(to: CChar.self, capacity: 12 * 20) { base in
                for i in 0..<min(ssCnt, 12) {
                    ssList.append(String(cString: base.advanced(by: i * 20)))
                }
            }
        }
        return BaziColumn(
            gan: extractString(from: c.gan, capacity: 8),
            zhi: extractString(from: c.zhi, capacity: 8),
            ganShiShen: extractString(from: c.gan_shi_shen, capacity: 12),
            nayin: extractString(from: c.nayin, capacity: 16),
            xingYun: extractString(from: c.xing_yun, capacity: 8),
            ziZuo: extractString(from: c.zi_zuo, capacity: 8),
            kongWang: extractString(from: c.kong_wang, capacity: 12),
            cangGan: cgList,
            shenSha: ssList,
            startYear: Int(c.start_year)
        )
    }

    private static func convertDayInfo(_ raw: SxwnlDayInfo) -> DayInfo {
        DayInfo(
            solarYear: Int(raw.solar_year),
            solarMonth: Int(raw.solar_month),
            solarDay: Int(raw.solar_day),
            lunarYear: Int(raw.lunar_year),
            lunarMonth: Int(raw.lunar_month),
            lunarDay: Int(raw.lunar_day),
            isLeapMonth: raw.is_leap_month,
            weekDay: Int(raw.week_day),
            yearGan: Int(raw.year_gan),
            yearZhi: Int(raw.year_zhi),
            monthGan: Int(raw.month_gan),
            monthZhi: Int(raw.month_zhi),
            dayGan: Int(raw.day_gan),
            dayZhi: Int(raw.day_zhi),
            jieQi: Int(raw.jie_qi),
            yueXiang: Int(raw.yue_xiang),
            constellation: Int(raw.constellation),
            jian12: Int(raw.jian12),
            yearGZ: extractString(from: raw.year_gz, capacity: 8),
            monthGZ: extractString(from: raw.month_gz, capacity: 8),
            dayGZ: extractString(from: raw.day_gz, capacity: 8),
            lunarMonthName: extractString(from: raw.lunar_month_name, capacity: 16),
            lunarDayName: extractString(from: raw.lunar_day_name, capacity: 12),
            jieQiName: extractString(from: raw.jie_qi_name, capacity: 12),
            jieQiTime: extractString(from: raw.jie_qi_time, capacity: 32),
            lipuJieQi: Int(raw.lipu_jie_qi),
            lipuJieQiName: extractString(from: raw.lipu_jie_qi_name, capacity: 12),
            shengXiao: extractString(from: raw.sheng_xiao, capacity: 8),
            constellationName: extractString(from: raw.constellation_name, capacity: 12),
            weekName: extractString(from: raw.week_name, capacity: 16),
            yueXiangName: extractString(from: raw.yue_xiang_name, capacity: 12),
            jian12Name: extractString(from: raw.jian12_name, capacity: 8),
            holiday: extractString(from: raw.holiday, capacity: 128),
            major: extractString(from: raw.major, capacity: 256),
            minor: extractString(from: raw.minor, capacity: 256),
            misc: extractString(from: raw.misc, capacity: 64),
            isOffDay: raw.is_off_day,
            lunarJunYear: Int(raw.lunar_jun_year),
            lunarLichunYear: Int(raw.lunar_lichun_year),
            huangdiYear: Int(raw.huangdi_year),
            moslemYear: Int(raw.moslem_year),
            moslemMonth: Int(raw.moslem_month),
            moslemDay: Int(raw.moslem_day),
            yearNaYin: extractString(from: raw.year_nayin, capacity: 16),
            monthNaYin: extractString(from: raw.month_nayin, capacity: 16),
            dayNaYin: extractString(from: raw.day_nayin, capacity: 16),
            yueXiangTime: extractString(from: raw.yue_xiang_time, capacity: 32),
            julianDay: raw.julian_day
        )
    }

    /// 从定长 C 字符串元组(tuple)读取 UTF-8 字符串
    private static func extractString<T>(from tuple: T, capacity: Int) -> String {
        withUnsafePointer(to: tuple) { ptr in
            ptr.withMemoryRebound(to: CChar.self, capacity: capacity) {
                String(cString: $0)
            }
        }
    }

    /// 写入 Swift 字符串到定长 C 字符串 tuple
    private static func writeCString<T>(_ s: String, into tuple: inout T, capacity: Int) {
        let bytes = s.utf8CString
        bytes.withUnsafeBufferPointer { buf in
            let count = min(buf.count, capacity - 1)
            withUnsafeMutablePointer(to: &tuple) { ptr in
                ptr.withMemoryRebound(to: CChar.self, capacity: capacity) { dest in
                    if count > 0 {
                        buf.baseAddress!.withMemoryRebound(to: CChar.self, capacity: count) { src in
                            _ = memcpy(dest, src, count)
                        }
                    }
                    dest[count] = 0
                }
            }
        }
    }

    /// 安全 C 字符串
    private static func cstr(_ p: UnsafePointer<CChar>?) -> String {
        guard let p = p else { return "" }
        return String(cString: p)
    }

    // ═══ Almanac (老黄历) ═════════════════════════════════════
    //
    //  SxwnlAlmanac 是 ~13KB POD struct; 在 leaf 函数中栈分配安全
    //  (iOS 主线程默认 1MB 栈, 子线程 512KB). 不要在深递归路径中调用.

    /// 取公历某日老黄历完整数据 (二十八宿/黄道黑道/冲煞/方位/彭祖/神煞/宜忌/吉时/用事).
    static func getAlmanac(year: Int, month: Int, day: Int) -> DayAlmanac? {
        var raw = SxwnlAlmanac()
        guard sxwnl_get_almanac(Int32(year), Int32(month), Int32(day), &raw) == 0 else {
            return nil
        }
        return convertAlmanac(&raw)
    }

    /// 老黄历静态知识 (董公总论/口诀/方位 等, 全局只取一次).
    static func getAlmanacTopics() -> [AlmanacTopic] {
        let kMax: Int32 = 64
        var buffer = [SxwnlAlmanacTopic](repeating: SxwnlAlmanacTopic(),
                                         count: Int(kMax))
        let n = sxwnl_get_almanac_topics(&buffer, kMax)
        guard n > 0 else { return [] }
        // 严格截断, 防越界
        let nn = min(Int(n), Int(kMax))
        return buffer.prefix(nn).map { t in
            AlmanacTopic(
                category: extractString(from: t.category, capacity: 32),
                title: extractString(from: t.title, capacity: 64),
                body: extractString(from: t.body, capacity: 1024)
            )
        }
    }

    // ═══ Geo (地理目录) ═══════════════════════════════════════

    static func geoListProvinces() -> [GeoProvince] {
        let kMax: Int32 = 64
        var buffer = [SxwnlGeoProvince](repeating: SxwnlGeoProvince(),
                                        count: Int(kMax))
        let n = sxwnl_geo_list_provinces(&buffer, kMax)
        guard n > 0 else { return [] }
        let nn = min(Int(n), Int(kMax))
        return buffer.prefix(nn).map { p in
            GeoProvince(
                province: extractString(from: p.province, capacity: 64),
                cityCount: Int(p.city_count)
            )
        }
    }

    static func geoListCities(province: String) -> [GeoCity] {
        // 单省最多约 300 城, 给 512 余量. 改用 Array 让 Swift 接管内存
        // 生命周期, 避免手动 allocate/deinitialize/deallocate 顺序出错 (POD
        // 类型 deinit 是 no-op, 但 Swift 严格模式要求成对调用).
        let kMax: Int32 = 512
        var buffer = [SxwnlGeoCity](repeating: SxwnlGeoCity(), count: Int(kMax))
        let n = province.withCString { cstr in
            sxwnl_geo_list_cities(cstr, &buffer, kMax)
        }
        guard n > 0 else { return [] }
        let nn = min(Int(n), Int(kMax))
        return buffer.prefix(nn).map { convertGeoCity($0) }
    }

    static func geoSearch(keyword: String, limit: Int = 64) -> [GeoCity] {
        let lim = Int32(max(1, min(limit, 1024)))
        var buffer = [SxwnlGeoCity](repeating: SxwnlGeoCity(), count: Int(lim))
        let n = keyword.withCString { cstr in
            sxwnl_geo_search(cstr, &buffer, lim)
        }
        guard n > 0 else { return [] }
        let nn = min(Int(n), Int(lim))
        return buffer.prefix(nn).map { convertGeoCity($0) }
    }

    static func geoListTimezones() -> [TimezoneGroup] {
        let kMax: Int32 = 128
        var buffer = [SxwnlTimezoneGroup](repeating: SxwnlTimezoneGroup(),
                                          count: Int(kMax))
        let n = sxwnl_geo_list_timezone_groups(&buffer, kMax)
        guard n > 0 else { return [] }
        let nn = min(Int(n), Int(kMax))
        return buffer.prefix(nn).map { g -> TimezoneGroup in
            // cities 是 [8][32] 二维数组, 取前 city_count 项
            var citiesList: [String] = []
            let cc = max(0, min(Int(g.city_count), 8))
            var citiesTuple = g.cities
            withUnsafePointer(to: &citiesTuple) { ptr in
                ptr.withMemoryRebound(to: CChar.self, capacity: 8 * 32) { base in
                    for k in 0..<cc {
                        citiesList.append(String(cString: base.advanced(by: k * 32)))
                    }
                }
            }
            return TimezoneGroup(
                continent: extractString(from: g.continent, capacity: 32),
                country: extractString(from: g.country, capacity: 96),
                timezone: g.timezone,
                cities: citiesList
            )
        }
    }

    static func geoDefault() -> GeoCity? {
        var c = SxwnlGeoCity()
        guard sxwnl_geo_default(&c) == 0 else { return nil }
        return convertGeoCity(c)
    }

    // ── Almanac / Geo 内部转换器 ──────────────────────────────

    private static func convertGeoCity(_ c: SxwnlGeoCity) -> GeoCity {
        GeoCity(
            province: extractString(from: c.province, capacity: 64),
            district: extractString(from: c.district, capacity: 64),
            longitude: c.longitude,
            latitude: c.latitude,
            timezone: c.timezone
        )
    }

    private static func convertAlmanac(_ a: UnsafePointer<SxwnlAlmanac>) -> DayAlmanac {
        let p = a.pointee

        // ── quotes ──
        var quotes: [AlmanacQuote] = []
        let qn = max(0, min(Int(p.quote_count), 4))
        var qTuple = p.quotes
        withUnsafePointer(to: &qTuple) { ptr in
            ptr.withMemoryRebound(to: SxwnlAlmanacQuote.self, capacity: 4) { base in
                for i in 0..<qn {
                    let q = base[i]
                    quotes.append(AlmanacQuote(
                        source: extractString(from: q.source, capacity: 24),
                        title: extractString(from: q.title, capacity: 40),
                        luck: extractString(from: q.luck, capacity: 8),
                        body: extractString(from: q.body, capacity: 1024)
                    ))
                }
            }
        }

        // ── shenSha ──
        var shen: [ShenSha] = []
        let sn = max(0, min(Int(p.shen_sha_count), 24))
        var sTuple = p.shen_sha
        withUnsafePointer(to: &sTuple) { ptr in
            ptr.withMemoryRebound(to: SxwnlShenSha.self, capacity: 24) { base in
                for i in 0..<sn {
                    let s = base[i]
                    shen.append(ShenSha(
                        name: extractString(from: s.name, capacity: 32),
                        isLucky: s.is_lucky,
                        weight: Int(s.weight)
                    ))
                }
            }
        }

        // ── yi / ji (定长 16 字节字符串 × 5) ──
        let yiList = extractStringList(tuple: p.yi,
                                       capacity: 5,
                                       itemLen: 16,
                                       count: max(0, min(Int(p.yi_count), 5)))
        let jiList = extractStringList(tuple: p.ji,
                                       capacity: 5,
                                       itemLen: 16,
                                       count: max(0, min(Int(p.ji_count), 5)))

        // ── luckyHours ──
        var hours: [LuckyHour] = []
        let hn = max(0, min(Int(p.lucky_hour_count), 8))
        var hTuple = p.lucky_hours
        withUnsafePointer(to: &hTuple) { ptr in
            ptr.withMemoryRebound(to: SxwnlLuckyHour.self, capacity: 8) { base in
                for i in 0..<hn {
                    let h = base[i]
                    hours.append(LuckyHour(
                        name: extractString(from: h.name, capacity: 12),
                        zhi: Int(h.zhi)
                    ))
                }
            }
        }

        // ── events (用事择吉) ──
        var events: [EventAdvice] = []
        let en = max(0, min(Int(p.event_count), 8))
        var eTuple = p.events
        withUnsafePointer(to: &eTuple) { ptr in
            ptr.withMemoryRebound(to: SxwnlEventAdvice.self, capacity: 8) { base in
                for i in 0..<en {
                    let e = base[i]
                    events.append(EventAdvice(
                        event: extractString(from: e.event, capacity: 16),
                        suitable: e.suitable,
                        reason: extractString(from: e.reason, capacity: 40)
                    ))
                }
            }
        }

        // ── notes (定长 80 字节字符串 × 4) ──
        let notes = extractStringList(tuple: p.notes,
                                      capacity: 4,
                                      itemLen: 80,
                                      count: max(0, min(Int(p.note_count), 4)))

        return DayAlmanac(
            xiu: extractString(from: p.xiu, capacity: 8),
            xiuZheng: extractString(from: p.xiu_zheng, capacity: 4),
            xiuAnimal: extractString(from: p.xiu_animal, capacity: 8),
            xiuLuck: extractString(from: p.xiu_luck, capacity: 4),
            xiuGong: extractString(from: p.xiu_gong, capacity: 4),
            twelveGod: extractString(from: p.twelve_god, capacity: 8),
            huangHei: extractString(from: p.huang_hei, capacity: 8),
            isHuangDao: p.is_huang_dao,
            chongShengXiao: extractString(from: p.chong_sheng_xiao, capacity: 8),
            chongGanZhi: extractString(from: p.chong_gan_zhi, capacity: 8),
            sha: extractString(from: p.sha, capacity: 8),
            xiShenFang: extractString(from: p.xi_shen_fang, capacity: 8),
            yangGuiFang: extractString(from: p.yang_gui_fang, capacity: 8),
            yinGuiFang: extractString(from: p.yin_gui_fang, capacity: 8),
            fuShenFang: extractString(from: p.fu_shen_fang, capacity: 8),
            caiShenFang: extractString(from: p.cai_shen_fang, capacity: 8),
            pengZuGan: extractString(from: p.peng_zu_gan, capacity: 64),
            pengZuZhi: extractString(from: p.peng_zu_zhi, capacity: 64),
            quotes: quotes,
            shenSha: shen,
            yi: yiList,
            ji: jiList,
            luckyHours: hours,
            events: events,
            notes: notes
        )
    }

    /// 从 Swift tuple 形式的"定长 char[itemLen][capacity]"中读取前 count 个字符串
    private static func extractStringList<T>(tuple: T, capacity: Int,
                                              itemLen: Int, count: Int) -> [String] {
        var tup = tuple
        var out: [String] = []
        out.reserveCapacity(count)
        withUnsafePointer(to: &tup) { ptr in
            ptr.withMemoryRebound(to: CChar.self, capacity: capacity * itemLen) { base in
                for i in 0..<count {
                    out.append(String(cString: base.advanced(by: i * itemLen)))
                }
            }
        }
        return out
    }
}
