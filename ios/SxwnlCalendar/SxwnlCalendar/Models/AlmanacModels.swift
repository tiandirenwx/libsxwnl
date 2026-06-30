import Foundation

// ════════════════════════════════════════════════════════════════
//  AlmanacModels — 老黄历 + 地理目录数据模型
//
//  与鸿蒙 NativeBridge.ets / Android Models.kt 字段一一对齐, 全部为
//  纯 value type. 数组按 *_count 严格截断 (Bridge 层保证), 永不暴露
//  上界外的未初始化数据.
// ════════════════════════════════════════════════════════════════

// MARK: - 老黄历

/// 择日典籍语录 (董公择日要诀 / 玉匣记 / 通胜经 ...)
struct AlmanacQuote {
    let source: String      // 典籍来源
    let title: String       // 段标题
    let luck: String        // "吉"/"凶"/"平"/"混"/""
    let body: String        // 原文
}

/// 神煞 (天德/月厌大祸/三合 ...)
struct ShenSha {
    let name: String
    let isLucky: Bool       // 吉神=true / 凶神=false
    let weight: Int         // 1一般 2中 3大煞
}

/// 吉时
struct LuckyHour {
    let name: String        // "福德"/"凤辇"/"贵人(阳)"
    let zhi: Int            // 0..11
}

/// 用事择吉建议
struct EventAdvice {
    let event: String       // "动土"/"上梁"/"安床" 等
    let suitable: Bool
    let reason: String
}

/// 单日老黄历完整数据
struct DayAlmanac {
    // 二十八宿
    let xiu: String
    let xiuZheng: String
    let xiuAnimal: String
    let xiuLuck: String
    let xiuGong: String

    // 黄道黑道
    let twelveGod: String
    let huangHei: String
    let isHuangDao: Bool

    // 冲煞
    let chongShengXiao: String
    let chongGanZhi: String
    let sha: String

    // 五吉神方位
    let xiShenFang: String
    let yangGuiFang: String
    let yinGuiFang: String
    let fuShenFang: String
    let caiShenFang: String

    // 彭祖百忌
    let pengZuGan: String
    let pengZuZhi: String

    // 典籍语录 / 神煞 / 宜忌 / 吉时 / 用事 / 备注
    let quotes: [AlmanacQuote]
    let shenSha: [ShenSha]
    let yi: [String]
    let ji: [String]
    let luckyHours: [LuckyHour]
    let events: [EventAdvice]
    let notes: [String]
}

/// 静态知识 (董公总论/口诀/方位 等)
struct AlmanacTopic: Identifiable {
    let id = UUID()
    let category: String    // "总论"/"基础理论"/"建筑"/"口诀"/"方位"
    let title: String
    let body: String
}

// MARK: - 地理目录 (省/市 + 国际时区)
//
//   数据来源: libsxwnl 内部 src/geo.cpp 的 JWv (国内 4000+ 市县) +
//   SQv (国际时区). 上层无需重复维护. Bridge 层一次性懒加载即可.

struct GeoProvince: Hashable, Identifiable {
    let id = UUID()
    let province: String
    let cityCount: Int
}

struct GeoCity: Hashable, Identifiable {
    let id = UUID()
    let province: String
    let district: String
    let longitude: Double   // 度, 东+ 西-
    let latitude: Double    // 度, 北+ 南-
    let timezone: Double    // 小时, 东+ 西- (国内统一 8)
}

struct TimezoneGroup: Hashable, Identifiable {
    let id = UUID()
    let continent: String   // "亚洲"/"欧洲"/...
    let country: String     // "日本"/"加拿大东部时区"/...
    let timezone: Double    // 标准时偏移 (不含 DST)
    let cities: [String]    // 代表城市名, 至多 8 个
}
