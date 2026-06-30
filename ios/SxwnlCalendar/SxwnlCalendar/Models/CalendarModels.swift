import Foundation
import SwiftUI

// ════════════════════════════════════════════════════════════════
//  对应鸿蒙 BaziInput.ets / 月历页面数据模型
// ════════════════════════════════════════════════════════════════

// MARK: - 月历 / 日单元

struct CalendarDay: Identifiable {
    let id = UUID()
    let info: DayInfo?
    let isCurrentMonth: Bool
    let isToday: Bool
}

struct MonthData {
    let year: Int
    let month: Int
    let days: [CalendarDay]
    let leapMonth: Int

    var title: String {
        "\(year)年\(month)月"
    }
}

// MARK: - 月相 / 节气事件 (用于底栏摘要)

struct MonthEvent: Identifiable {
    let id = UUID()
    let day: Int        // 公历日
    let time: String    // "HH:MM:SS"
    let name: String    // 朔月/望月/上弦/下弦 或 节气名
}

// MARK: - 历法类型 (对应 C++ LiFaType)
//   旧 BaziView 使用此枚举; 新 UI 使用更精简的 LiFaModern。

enum LiFaType: Int, CaseIterable {
    case shixian = 0
    case daming = 1
    case dayan = 2
    case xuanming = 3
    case yitian = 4
    case chongzhen = 5
    case shoushi = 6
    case tongtian = 7
    case chongtian = 8
    case guantian = 9
    case jiyuan = 10
    case jingchu = 11
    case pingqi = 12

    var name: String {
        switch self {
        case .shixian:    return "时宪历"
        case .daming:     return "大明历"
        case .dayan:      return "大衍历"
        case .xuanming:   return "宣明历"
        case .yitian:     return "仪天历"
        case .chongzhen:  return "崇祯历"
        case .shoushi:    return "授时历"
        case .tongtian:   return "统天历"
        case .chongtian:  return "崇天历"
        case .guantian:   return "观天历"
        case .jiyuan:     return "纪元历"
        case .jingchu:    return "景初历"
        case .pingqi:     return "平气"
        }
    }
}

// MARK: - 现代历法 (与鸿蒙 BaziInput.ets 一致)
//   仅 3 项, 用于新版八字表单。

enum LiFaModern: Int, CaseIterable, Identifiable {
    case dingQi = 11        // 现代农历定气
    case pingXiaZhi = 13    // 平气定夏至
    case pingDongZhi = 12   // 平气定冬至

    var id: Int { rawValue }
    var label: String {
        switch self {
        case .dingQi:      return "定气"
        case .pingXiaZhi:  return "平气定夏至"
        case .pingDongZhi: return "平气定冬至"
        }
    }
}

// MARK: - 出生时间选择 (公历 / 农历 / 反推)

enum BirthInputMode: Int {
    case solar = 0       // 公历
    case lunar = 1       // 农历
    case reverse = 2     // 四柱反推 (结果统一归一为公历)
}

struct DateSelection {
    var inputMode: BirthInputMode = .solar
    var year: Int = 2000
    var month: Int = 1
    var day: Int = 1
    var hour: Int = 12
    var minute: Int = 0
    var isLeap: Bool = false
    var isSpec: Bool = false
}

// MARK: - 八字表单
//   字段同时兼容旧 BaziView 的平铺访问与新 UI 的结构化访问

struct BaziFormData {
    var name: String = ""
    var gender: Bool = false
    var year: Int = 1990
    var month: Int = 1
    var day: Int = 1
    var hour: Int = 12
    var minute: Int = 0
    var isLunar: Bool = false
    var isLeap: Bool = false
    var lifa: LiFaType = .shixian
    var isAst: Bool = false
    var longitude: Double = 116.4
    var latitude: Double = 39.9

    // ── 新 UI 字段 (与鸿蒙 BaziInput.ets 对齐) ──
    var lifaModern: LiFaModern = .dingQi
    var astEnabled: Bool { get { isAst } set { isAst = newValue } }
    var inputMode: BirthInputMode = .solar
    var isSpec: Bool = false

    /// 派生: 当前是否反推
    var isReverse: Bool { inputMode == .reverse }

    /// 同步派生 DateSelection (用于新 UI / BaziCalc)
    var selection: DateSelection {
        get {
            DateSelection(
                inputMode: inputMode,
                year: year, month: month, day: day,
                hour: hour, minute: minute,
                isLeap: isLeap, isSpec: isSpec
            )
        }
        set {
            inputMode = newValue.inputMode
            year = newValue.year
            month = newValue.month
            day = newValue.day
            hour = newValue.hour
            minute = newValue.minute
            isLeap = newValue.isLeap
            isSpec = newValue.isSpec
            isLunar = newValue.inputMode == .lunar
        }
    }
}
