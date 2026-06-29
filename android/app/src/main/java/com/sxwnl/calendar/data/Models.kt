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

    var yearGz: String = "",
    var monthGz: String = "",
    var dayGz: String = "",
    var lunarMonthName: String = "",
    var lunarDayName: String = "",
    var jieQiName: String = "",
    var jieQiTime: String = "",
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
    var dayName: String = ""
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
    var wuXingStatus: Array<String> = arrayOf("", "", "", "", ""),
    var siLing: String = "",
    var liuNian: List<LiuNianItem> = emptyList()
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
