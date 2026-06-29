import sxwnl from 'libsxwnl_napi.so'

export interface DayInfo {
  solarYear: number
  solarMonth: number
  solarDay: number
  lunarYear: number
  lunarMonth: number
  lunarDay: number
  isLeapMonth: boolean
  weekDay: number
  yearGan: number
  yearZhi: number
  monthGan: number
  monthZhi: number
  dayGan: number
  dayZhi: number
  jieQi: number
  yueXiang: number
  constellation: number
  jian12: number
  yearGZ: string
  monthGZ: string
  dayGZ: string
  lunarMonthName: string
  lunarDayName: string
  jieQiName: string
  jieQiTime: string
  shengXiao: string
  constellationName: string
  weekName: string
  yueXiangName: string
  yueXiangTime: string
  jian12Name: string

  // 节日 & 杂节(数九/三伏等)
  holiday: string      // A 类: 法定/重要节日
  major: string        // B 类: 中等节日
  minor: string        // C 类: 一般节日
  misc: string         // 杂节
  isOffDay: boolean    // 法定/周末放假

  // 纪年扩展
  lunarJunYear: number     // 农历纪年(春节为界)
  lunarLichunYear: number  // 农历纪年(立春为界)
  huangdiYear: number      // 黄帝纪年

  // 回历(Hijri)
  moslemYear: number
  moslemMonth: number
  moslemDay: number

  // 纳音
  yearNaYin: string
  monthNaYin: string
  dayNaYin: string

  // 儒略日(对应当日12点, 与 sxwnl/lunar.js 一致)
  julianDay: number
}

// ── 年历: 一年的农历月+节气信息 (对应 JS nianLi2HTML) ───────
export interface YearCalJieQi {
  idx: number          // 节气索引 0..23 (0=冬至)
  name: string         // 节气名称
  gz: string           // 节气日干支
  solarMonth: number
  solarDay: number
  time: string         // "HH:MM:SS" 精确时刻
  dayOffset: number    // 距本月初一的天数
  dayName: string      // 月内日名 (初一/十五...)
}

export interface YearCalMonth {
  monthIdx: number     // 农历月名索引 0..11 (0=正)
  monthName: string    // "正月"/"闰二月"/"后九月"
  isLeap: boolean
  isSpec: boolean
  dayCount: number     // 月大小 29/30
  solarYear: number    // 朔日的公历年
  solarMonth: number
  solarDay: number
  shuoGZ: string       // 朔日干支
  jieQi: YearCalJieQi[]
}

export interface CangGanItem {
  gan: string
  shiShen: string
}

export interface BaziColumn {
  gan: string
  zhi: string
  ganShiShen: string
  nayin: string
  xingYun: string
  ziZuo: string
  kongWang: string
  cangGan: CangGanItem[]
  shenSha: string[]
  startYear: number
}

export interface LiuNianItem {
  year: number
  age: number
  ganZhi: string
  ganShiShen: string
  zhiShiShen: string
  xiaoYun: string
  xiaoYunShiShen: string
}

export interface BaziResult {
  userName: string
  gender: string
  solarBirth: string
  lunarBirth: string
  dateOfBirth: string
  shengXiao: string
  age: string
  lifa: string
  dingQiType: string
  ast: string
  jieQi: string
  qiYun: string
  jiaoYun: string
  columns: BaziColumn[]
  currentDaYun?: BaziColumn
  currentLiuNian?: BaziColumn
  daYunColumns: BaziColumn[]
  wuXingCount: number[]
  wuXingStatus: string[]
  siLing: string
  liuNian: LiuNianItem[]
  // 每个大运对应的 10 年流年/小运 (与 daYunColumns 等长, 按相同顺序)
  // 用于专业版"切换大运 → 展示对应流年"
  liuNianAll?: LiuNianItem[][]
}

export interface BaziParams {
  year: number
  month: number
  day: number
  hour: number
  minute: number
  second: number
  gender: number
  isLunar: boolean
  isLeapMonth: boolean
  isSpecMonth: boolean
  name: string
  astEnabled: boolean
  longitude: number
  latitude: number
  lifa: number
}

export interface LunarMonth {
  month: number
  isLeap: boolean
  isSpec: boolean
  name: string
}

export interface ReverseItem {
  year: number
  month: number
  day: number
  hour: number
  yearGZ: string
  monthGZ: string
  dayGZ: string
  hourGZ: string
}

export interface JieQiItem {
  name: string
  time: string
  month: number
  day: number
}

export interface SolarDate {
  year: number
  month: number
  day: number
}

// ── 日月食可视化数据 ─────────────────────────────────────

export interface EclipseGeoPoint {
  longitude: number  // 度, 东正
  latitude: number   // 度, 北正
  jd: number         // 力学时 JD
}

export interface SolarEclipseItem {
  jd: number
  year: number
  month: number
  day: number
  hour: number
  minute: number
  type: string       // T/A/P/H/T0/A0/T1/A1/H2/H3
  typeName: string   // 全食/环食/偏食 等
  gamma: number
  magnitude: number  // 食分
  centerLon: number  // 100 表示无中心
  centerLat: number
  pathWidth: number  // km
  duration: number   // s
  saros: number      // Saros 周期序号
  sarosMember: number
  season: string     // 食季
}

export interface SolarEclipsePath {
  type: string
  typeName: string
  gamma: number
  magnitude: number
  pathWidthKm: number
  durationSeconds: number
  maxEclipseJd: number
  centerLine: EclipseGeoPoint[]
  umbraNorth: EclipseGeoPoint[]
  umbraSouth: EclipseGeoPoint[]
  penumbraNorth: EclipseGeoPoint[]
  penumbraSouth: EclipseGeoPoint[]
  partialStart: EclipseGeoPoint
  centralStart: EclipseGeoPoint
  greatestEclipse: EclipseGeoPoint
  centralEnd: EclipseGeoPoint
  partialEnd: EclipseGeoPoint
}

export interface EclipseFrame {
  sunRadius: number   // 弧度
  moonX: number       // 月心相对日心的 x 偏移 (弧度)
  moonY: number
  moonRadius: number  // 弧度
  jd: number
  magnitude: number
}

export interface LocalSolarEclipse {
  type: string         // 偏/全/环
  maxMagnitude: number
  moonSunRatio: number
  frames: EclipseFrame[]
  tC1: number   // 初亏 JD (0 = N/A)
  tC2: number   // 食既/环食始
  tMax: number  // 食甚
  tC3: number   // 生光/环食终
  tC4: number   // 复圆
  tSunrise: number
  tSunset: number
}

export interface LunarEclipseItem {
  jd: number
  year: number
  month: number
  day: number
  hour: number
  minute: number
  type: string         // N/B/P/T
  typeName: string     // 半影食/偏食/全食
  magnitude: number
  saros: number        // Saros 周期序号
  sarosMember: number
  season: string       // 食季
}

export interface LunarEclipseFrame {
  moonX: number          // 月心相对地影中心 (弧度)
  moonY: number
  moonRadius: number
  umbraRadius: number    // 本影半径
  penumbraRadius: number // 半影半径
  jd: number
  coverage: number       // 月亮直径被覆盖比例
}

export interface LunarEclipseDetail {
  type: string
  typeName: string
  maxMagnitude: number
  frames: LunarEclipseFrame[]
  tP1: number   // 半影食始
  tU1: number   // 初亏
  tU2: number   // 食既
  tMax: number
  tU3: number   // 生光
  tU4: number   // 复圆
  tP4: number   // 半影食终
}

// ── 日月升降/中天 (单天, 对应 JS SZJ.calcRTS) ─────────────
export interface DayRTS {
  sunRise: string       // 日出 HH:MM:SS
  sunSet: string        // 日落
  sunMeridian: string   // 日上中天
  moonRise: string      // 月出, "--:--:--" 表示当日无升起
  moonSet: string       // 月落
  moonMeridian: string  // 月中天
  civilDawn: string     // 晨起天亮 (民用晨)
  civilDusk: string     // 晚上天黑 (民用昏)
  dayLength: string     // 日照时间 (日出→日落)
  lightLength: string   // 白天时间 (民用晨→民用昏)
}

export class NativeBridge {
  static getDayInfo(year: number, month: number, day: number): DayInfo {
    return sxwnl.getDayInfo(year, month, day) as DayInfo
  }

  static getMonthData(year: number, month: number): DayInfo[] {
    return sxwnl.getMonthData(year, month) as DayInfo[]
  }

  static lunarToSolar(year: number, month: number, day: number, isLeap: boolean): SolarDate {
    return sxwnl.lunarToSolar(year, month, day, isLeap) as SolarDate
  }

  static solarToLunar(year: number, month: number, day: number): DayInfo {
    return sxwnl.solarToLunar(year, month, day) as DayInfo
  }

  static getJieQiList(year: number): JieQiItem[] {
    return sxwnl.getJieQiList(year) as JieQiItem[]
  }

  static getYearLeapMonth(year: number): number {
    return sxwnl.getYearLeapMonth(year) as number
  }

  static getLunarMonths(year: number): LunarMonth[] {
    return (sxwnl.getLunarMonths(year) ?? []) as LunarMonth[]
  }

  static getLunarMonthDays(year: number, month: number, isLeap: boolean, isSpec: boolean): number {
    return sxwnl.getLunarMonthDays(year, month, isLeap, isSpec) as number
  }

  // 公历某年月真实存在的日号(升序); 自动跳过 1582-10 改历缺口, 并给出正确月长
  static getSolarMonthValidDays(year: number, month: number): number[] {
    return (sxwnl.getSolarMonthValidDays(year, month) ?? []) as number[]
  }

  // 年历: 按农历月聚合, 含朔日干支/月大小/本月节气, 用于年视图展示
  static getYearCalendar(year: number): YearCalMonth[] {
    return (sxwnl.getYearCalendar(year) ?? []) as YearCalMonth[]
  }

  static calcBazi(params: BaziParams): BaziResult {
    return sxwnl.calcBazi(params) as BaziResult
  }

  // 四柱反推: 干支用索引(天干0-9, 地支0-11); 时柱未知传 hg=hz=-1
  static baziReverse(yg: number, yz: number, mg: number, mz: number,
    dg: number, dz: number, hg: number, hz: number,
    startYear: number, endYear: number): ReverseItem[] {
    const r = sxwnl.baziReverse(yg, yz, mg, mz, dg, dz, hg, hz, startYear, endYear)
    return (r ?? []) as ReverseItem[]
  }

  // 日月升降/中天 (单天, 对应 JS SZJ.calcRTS)
  //   longitude/latitude: 度, 东经/北纬为正
  //   tzHours:            所在时区偏移(小时, 东向为正), 北京 = 8
  static calcDayRTS(year: number, month: number, day: number,
    longitude: number, latitude: number, tzHours: number): DayRTS | null {
    const r = sxwnl.calcDayRTS(year, month, day, longitude, latitude, tzHours)
    return r === null ? null : r as DayRTS
  }

  // ── 日月食 (结构化, 用于可视化) ─────────────────────────

  // 从 (year, month) 开始, 最多返回 count 条日食 (含偏食以上)
  static searchSolarEclipses(year: number, month: number, count: number): SolarEclipseItem[] {
    return (sxwnl.searchSolarEclipses(year, month, count) ?? []) as SolarEclipseItem[]
  }

  static getSolarEclipsePath(year: number, month: number, day: number): SolarEclipsePath | null {
    const r = sxwnl.getSolarEclipsePath(year, month, day)
    return r === null ? null : r as SolarEclipsePath
  }

  static getLocalSolarEclipse(year: number, month: number, day: number,
    longitude: number, latitude: number, frameInterval: number): LocalSolarEclipse | null {
    const r = sxwnl.getLocalSolarEclipse(year, month, day, longitude, latitude, frameInterval)
    return r === null ? null : r as LocalSolarEclipse
  }

  // 从 (year, month) 开始, 最多返回 count 条月食 (含半影食以上)
  static searchLunarEclipses(year: number, month: number, count: number): LunarEclipseItem[] {
    return (sxwnl.searchLunarEclipses(year, month, count) ?? []) as LunarEclipseItem[]
  }

  static getLunarEclipseDetail(year: number, month: number, day: number, frameInterval: number): LunarEclipseDetail | null {
    const r = sxwnl.getLunarEclipseDetail(year, month, day, frameInterval)
    return r === null ? null : r as LunarEclipseDetail
  }

  // 世界地图海岸线轮廓; number[] 普通数组, 经纬度(弧度)交替, 段间用 1e7 分隔.
  // 注意: HarmonyOS Ark 引擎不认 NAPI 创建的 Float64Array (会抛
  // "does not have [[TypedArrayName]] internal slot"), 必须用 number[].
  static getWorldMapDitu0(): number[] {
    const r = sxwnl.getWorldMapDitu0()
    return r === null ? [] : (r as number[])
  }

  // idx=1 → ditu1 大图海岸; idx=2 → ditu2 国界. 首次调用自动加载内建数据.
  static getWorldMapData(idx: number): number[] {
    const r = sxwnl.getWorldMapData(idx)
    return r === null ? [] : (r as number[])
  }
}
