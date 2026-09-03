package com.sxwnl.calendar.data

// ════════════════════════════════════════════════════════════════
//  数据模型 — 与鸿蒙端 NativeBridge.ets 对齐
//
//  所有 data class 均通过 @JvmOverloads + 全默认值参数, 生成对 JNI
//  友好的无参构造器, 供 jni_bridge.cpp 用 GetMethodID("<init>", "()V")
//  快速创建空对象, 再用 SetXxxField 填充字段。
// ════════════════════════════════════════════════════════════════

// ─── 单日完整信息 ────────────────────────────────────────────────

data class DayInfo @JvmOverloads constructor(
    var solarYear: Int = 0,
    var solarMonth: Int = 0,
    var solarDay: Int = 0,

    var lunarYear: Int = 0,
    var lunarMonth: Int = 0,
    var lunarDay: Int = 0,
    var isLeapMonth: Boolean = false,

    var weekDay: Int = 0,

    var yearGan: Int = 0,
    var yearZhi: Int = 0,
    var monthGan: Int = 0,
    var monthZhi: Int = 0,
    var dayGan: Int = 0,
    var dayZhi: Int = 0,

    var jieQi: Int = -1,
    var yueXiang: Int = -1,
    var constellation: Int = 0,
    var jian12: Int = 0,

    // 命名: yearGZ/monthGZ/dayGZ (大写 GZ) — 与 NAPI/iOS 字段保持一致, JNI
    //   会同名设置. 历史曾用 yearGz, 已统一.
    var yearGZ: String = "",
    var monthGZ: String = "",
    var dayGZ: String = "",
    var lunarMonthName: String = "",
    var lunarDayName: String = "",
    var jieQiName: String = "",
    var jieQiTime: String = "",
    // 历谱口径节气(整日表+QB, 对齐 sxwnl 网页版): 日历格子标签建议用此值,
    // 古代(1645年前)与天文口径 jieQi 可能差 1 天; 精确时刻仍用 jieQiTime。
    var lipuJieQi: Int = -1,
    var lipuJieQiName: String = "",
    var shengXiao: String = "",
    var constellationName: String = "",
    var weekName: String = "",
    var yueXiangName: String = "",
    var jian12Name: String = "",

    // ── 节日 & 杂节 ──────────────────────────
    var holiday: String = "",
    var major: String = "",
    var minor: String = "",
    var misc: String = "",
    var isOffDay: Boolean = false,

    // ── 纪年扩展 ─────────────────────────────
    var lunarJunYear: Int = 0,
    var lunarLichunYear: Int = 0,
    var huangdiYear: Int = 0,

    // ── 回历 ─────────────────────────────────
    var moslemYear: Int = 0,
    var moslemMonth: Int = 0,
    var moslemDay: Int = 0,

    // ── 纳音 ─────────────────────────────────
    var yearNaYin: String = "",
    var monthNaYin: String = "",
    var dayNaYin: String = "",

    // ── 月相极值时刻 ─────────────────────────
    var yueXiangTime: String = "",

    // ── 儒略日 (当日 12:00) ──────────────────
    var julianDay: Double = 0.0
)

// ─── 节气 / 农历月 ──────────────────────────────────────────────

data class JieQiItem @JvmOverloads constructor(
    var idx: Int = 0,
    var solarMonth: Int = 0,
    var solarDay: Int = 0,
    var name: String = "",
    var time: String = ""
)

data class LunarMonth @JvmOverloads constructor(
    var month: Int = 0,
    var isLeap: Boolean = false,
    var isSpec: Boolean = false,
    var name: String = ""
)

data class SolarDate @JvmOverloads constructor(
    var year: Int = 0,
    var month: Int = 0,
    var day: Int = 0
)

// ─── 年历: 按农历月聚合 ─────────────────────────────────────────

data class YearCalJieQi @JvmOverloads constructor(
    var idx: Int = 0,
    var name: String = "",
    var gz: String = "",
    var solarMonth: Int = 0,
    var solarDay: Int = 0,
    var time: String = "",
    var dayOffset: Int = 0,
    var dayName: String = "",
    // 精确交气(天文定气)所在公历日期; 古代可能与历谱 solarMonth/solarDay 差 1 天
    var accMonth: Int = 0,
    var accDay: Int = 0
)

class YearCalMonth @JvmOverloads constructor(
    var monthIdx: Int = 0,
    var monthName: String = "",
    var isLeap: Boolean = false,
    var isSpec: Boolean = false,
    var dayCount: Int = 0,
    var solarYear: Int = 0,
    var solarMonth: Int = 0,
    var solarDay: Int = 0,
    var shuoGZ: String = "",
    var jieQi: Array<YearCalJieQi> = emptyArray()
)

// ─── 日月升降 (RTS) ─────────────────────────────────────────────

data class DayRTS @JvmOverloads constructor(
    var sunRise: String = "",
    var sunSet: String = "",
    var sunMeridian: String = "",
    var moonRise: String = "",
    var moonSet: String = "",
    var moonMeridian: String = "",
    var civilDawn: String = "",
    var civilDusk: String = "",
    var dayLength: String = "",
    var lightLength: String = ""
)

// ─── 八字 ───────────────────────────────────────────────────────

data class CangGanItem @JvmOverloads constructor(
    var gan: String = "",
    var shiShen: String = ""
)

class BaziColumn @JvmOverloads constructor(
    var gan: String = "",
    var zhi: String = "",
    var ganShiShen: String = "",
    var nayin: String = "",
    var xingYun: String = "",
    var ziZuo: String = "",
    var kongWang: String = "",
    var cangGan: Array<CangGanItem> = emptyArray(),
    var shenSha: Array<String> = emptyArray(),
    var startYear: Int = 0
)

data class LiuNianItem @JvmOverloads constructor(
    var year: Int = 0,
    var age: Int = 0,
    var ganZhi: String = "",
    var ganShiShen: String = "",
    var zhiShiShen: String = "",
    var xiaoYun: String = "",
    var xiaoYunShiShen: String = ""
)

data class ReverseItem @JvmOverloads constructor(
    var year: Int = 0,
    var month: Int = 0,
    var day: Int = 0,
    var hour: Int = -1,
    var yearGZ: String = "",
    var monthGZ: String = "",
    var dayGZ: String = "",
    var hourGZ: String = ""
)

class BaziResult @JvmOverloads constructor(
    var userName: String = "",
    var gender: String = "",
    var solarBirth: String = "",
    var lunarBirth: String = "",
    var dateOfBirth: String = "",
    var shengXiao: String = "",
    var age: String = "",
    var lifa: String = "",
    var dingQiType: String = "",
    var ast: String = "",
    var jieQi: String = "",
    var qiYun: String = "",
    var jiaoYun: String = "",

    // ── 完整盘面 ─────────────────────────────
    var columns: List<BaziColumn> = emptyList(),
    var currentDaYun: BaziColumn? = null,
    var currentLiuNian: BaziColumn? = null,
    var daYunColumns: List<BaziColumn> = emptyList(),
    var wuXingCount: IntArray = IntArray(5),
    var     wuXingStatus: Array<String> = arrayOf("", "", "", "", ""),
    var siLing: String = "",
    var liuNian: List<LiuNianItem> = emptyList(),

    // 按大运分组的流年: 与 daYunColumns 一一对应, 每桶 10 个流年.
    //  对齐鸿蒙 NAPI 返回字段 (napi_bazi.cpp: result.v("liuNianAll", allLn)).
    //  rendering 时直接 daYunColumns[i] ↔ liuNianAll[i] 对应展开即可.
    var liuNianAll: List<List<LiuNianItem>> = emptyList()
)

// ─── 八字输入参数 ───────────────────────────────────────────────

data class BaziParams @JvmOverloads constructor(
    var name: String = "",
    var gender: Boolean = false,
    var lifa: Int = 11,
    var astEnabled: Boolean = false,
    var longitude: Double = 116.3833,
    var latitude: Double = 39.9,
    var inputMode: Int = 0,             // 0公历 1农历 2反推
    var year: Int = 2000,
    var month: Int = 1,
    var day: Int = 1,
    var hour: Int = 12,
    var minute: Int = 0,
    var isLeap: Boolean = false,
    var isSpec: Boolean = false
)

// ─── 老黄历 (Almanac) — 与鸿蒙 DayAlmanac 一致 ──────────────────

/** 择日典籍语录 (董公择日要诀 / 玉匣记 / 通胜经 ...) */
data class AlmanacQuote @JvmOverloads constructor(
    var source: String = "",            // 典籍来源
    var title: String = "",             // 段标题
    var luck: String = "",              // "吉"/"凶"/"平"/"混"/""
    var body: String = ""               // 原文
)

/** 神煞 (天德/月厌大祸/三合 ...) */
data class ShenSha @JvmOverloads constructor(
    var name: String = "",
    var isLucky: Boolean = true,
    var weight: Int = 1                 // 1一般 2中 3大煞
)

/** 吉时 */
data class LuckyHour @JvmOverloads constructor(
    var name: String = "",              // "福德"/"凤辇"/"贵人(阳)"
    var zhi: Int = 0                    // 0..11
)

/** 用事择吉建议 */
data class EventAdvice @JvmOverloads constructor(
    var event: String = "",             // "动土"/"上梁"/"安床" ...
    var suitable: Boolean = false,
    var reason: String = ""
)

/** 单日老黄历完整数据 */
data class DayAlmanac @JvmOverloads constructor(
    var xiu: String = "",
    var xiuZheng: String = "",
    var xiuAnimal: String = "",
    var xiuLuck: String = "",
    var xiuGong: String = "",

    var twelveGod: String = "",
    var huangHei: String = "",
    var isHuangDao: Boolean = false,

    var chongShengXiao: String = "",
    var chongGanZhi: String = "",
    var sha: String = "",

    var xiShenFang: String = "",
    var yangGuiFang: String = "",
    var yinGuiFang: String = "",
    var fuShenFang: String = "",
    var caiShenFang: String = "",

    var pengZuGan: String = "",
    var pengZuZhi: String = "",

    var quotes: Array<AlmanacQuote> = emptyArray(),
    var shenSha: Array<ShenSha> = emptyArray(),
    var yi: Array<String> = emptyArray(),
    var ji: Array<String> = emptyArray(),
    var luckyHours: Array<LuckyHour> = emptyArray(),
    var events: Array<EventAdvice> = emptyArray(),
    var notes: Array<String> = emptyArray()
)

/** 静态知识 (董公总论/口诀/方位 等) */
data class AlmanacTopic @JvmOverloads constructor(
    var category: String = "",          // "总论"/"基础理论"/"建筑"/"口诀"/"方位"
    var title: String = "",
    var body: String = ""
)

// ─── 地理目录 (省/市 + 国际时区) — 与鸿蒙 GeoCity 一致 ───────────

data class GeoProvince @JvmOverloads constructor(
    var province: String = "",
    var cityCount: Int = 0
)

data class GeoCity @JvmOverloads constructor(
    var province: String = "",
    var district: String = "",
    var longitude: Double = 0.0,
    var latitude: Double = 0.0,
    var timezone: Double = 8.0
)

data class TimezoneGroup @JvmOverloads constructor(
    var continent: String = "",
    var country: String = "",
    var timezone: Double = 0.0,
    var cities: Array<String> = emptyArray()
)
