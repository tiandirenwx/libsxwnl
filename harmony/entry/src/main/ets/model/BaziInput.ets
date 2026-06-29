/**
 * 八字输入模型与纯计算逻辑(无 UI 依赖, 便于后续跨平台复用)。
 *
 * 仅依赖 NativeBridge 提供的数据接口, 不引入任何 ArkUI 组件。
 */
import { NativeBridge, BaziParams, LunarMonth } from '../services/NativeBridge'

// ─── 历法类型(对应 C++ LiFaType 枚举值) ───
export const LIFA_DING_QI: number = 11        // 现代农历定气
export const LIFA_PING_DONG_ZHI: number = 12  // 平气·定冬至
export const LIFA_PING_XIA_ZHI: number = 13   // 平气·定夏至

export interface LiFaOption {
  label: string
  value: number
}

export const LIFA_OPTIONS: LiFaOption[] = [
  { label: '定气', value: LIFA_DING_QI },
  { label: '平气定夏至', value: LIFA_PING_XIA_ZHI },
  { label: '平气定冬至', value: LIFA_PING_DONG_ZHI }
]

// ─── 干支常量 ───
export const GAN: string[] = ['甲', '乙', '丙', '丁', '戊', '己', '庚', '辛', '壬', '癸']
export const ZHI: string[] = ['子', '丑', '寅', '卯', '辰', '巳', '午', '未', '申', '酉', '戌', '亥']

// 输入模式
export const MODE_SOLAR: number = 0   // 公历
export const MODE_LUNAR: number = 1   // 农历
export const MODE_REVERSE: number = 2 // 四柱反推

/** 出生时间选择结果(选择器回填给表单) */
export interface DateSelection {
  inputMode: number   // 0公历 1农历 (反推结果统一归一为公历)
  year: number
  month: number
  day: number
  hour: number
  minute: number
  isLeap: boolean
  isSpec: boolean
}

/** 八字表单数据(可序列化, 跨平台共享语义) */
export interface BaziForm {
  name: string
  gender: number        // 0男 1女
  lifa: number          // LiFaType
  astEnabled: boolean   // 真太阳时
  longitude: number
  latitude: number
  sel: DateSelection
}

/** 八字相关纯逻辑工具 */
export class BaziCalc {
  // 六十甲子序号 → 名称
  static jiaZiName(idx: number): string {
    const i: number = ((idx % 60) + 60) % 60
    return GAN[i % 10] + ZHI[i % 12]
  }

  // 干支索引 → 六十甲子序号
  static gz60(gan: number, zhi: number): number {
    for (let n = 0; n < 60; n++) {
      if (n % 10 === gan && n % 12 === zhi) return n
    }
    return 0
  }

  // 五虎遁: 由年干得 12 个合法月柱(寅月起)的六十甲子序号
  static validMonths(yearIdx: number): number[] {
    const yg: number = ((yearIdx % 60) + 60) % 60 % 10
    const base: number = ((yg % 5) * 2 + 2) % 10
    const r: number[] = []
    for (let k = 0; k < 12; k++) r.push(BaziCalc.gz60((base + k) % 10, (2 + k) % 12))
    return r
  }

  // 五鼠遁: 由日干得 12 个合法时柱(子时起)的六十甲子序号
  static validHours(dayIdx: number): number[] {
    const dg: number = ((dayIdx % 60) + 60) % 60 % 10
    const base: number = ((dg % 5) * 2) % 10
    const r: number[] = []
    for (let h = 0; h < 12; h++) r.push(BaziCalc.gz60((base + h) % 10, h))
    return r
  }

  // 在合法列表中按地支保留, 返回新的六十甲子序号
  static remapByZhi(validList: number[], oldIdx: number): number {
    const z: number = ((oldIdx % 60) + 60) % 60 % 12
    for (let i = 0; i < validList.length; i++) {
      if (validList[i] % 12 === z) return validList[i]
    }
    return validList[0]
  }

  // 流年(立春界)干支
  static lnGan(year: number): number { return (((year - 4) % 10) + 10) % 10 }
  static lnZhi(year: number): number { return (((year - 4) % 12) + 12) % 12 }

  // 年份解析: 支持 B212 / *212 / -211 表示公元前(返回天文年)
  static parseYear(s: string): number {
    const t: string = s.trim()
    if (t.length === 0) return NaN
    const c0: string = t.charAt(0)
    if (c0 === 'B' || c0 === 'b' || c0 === '*') {
      const n: number = Number(t.substring(1))
      if (isNaN(n)) return NaN
      return -(Math.floor(n) - 1)
    }
    const n: number = Number(t)
    if (isNaN(n)) return NaN
    return Math.floor(n)
  }

  // 天文年 → 显示文本
  static formatYear(y: number): string {
    if (y <= 0) return `公元前${1 - y}年`
    return `${y}年`
  }

  // 公历某月天数
  static solarDaysInMonth(year: number, month: number): number {
    if (month === 2) {
      const leap: boolean = (year % 4 === 0 && year % 100 !== 0) || (year % 400 === 0)
      return leap ? 29 : 28
    }
    if (month === 4 || month === 6 || month === 9 || month === 11) return 30
    return 31
  }

  // 农历某月天数(29/30)
  static lunarDaysInMonth(year: number, month: number, isLeap: boolean, isSpec: boolean): number {
    const d: number = NativeBridge.getLunarMonthDays(year, month, isLeap, isSpec)
    return (d === 29 || d === 30) ? d : 30
  }

  // 两位补零
  static pad2(n: number): string {
    return n < 10 ? `0${n}` : `${n}`
  }

  // 出生时间 → 一条记录文本(参考截图: "公历 2026年06月02日 13点16分")
  static formatRecord(sel: DateSelection, lunarMonths: LunarMonth[]): string {
    const ymd: string = `${BaziCalc.formatYear(sel.year)}${BaziCalc.pad2(sel.month)}月${BaziCalc.pad2(sel.day)}日`
    const hm: string = `${sel.hour}点${BaziCalc.pad2(sel.minute)}分`
    if (sel.inputMode === MODE_LUNAR) {
      let mName: string = `${sel.month}月`
      for (let i = 0; i < lunarMonths.length; i++) {
        if (lunarMonths[i].month === sel.month && lunarMonths[i].isLeap === sel.isLeap) {
          mName = lunarMonths[i].name
          break
        }
      }
      return `农历 ${BaziCalc.formatYear(sel.year)}${mName}${sel.day}日 ${hm}`
    }
    return `公历 ${ymd} ${hm}`
  }

  // 表单 → NativeBridge 入参
  static toParams(form: BaziForm): BaziParams {
    const s: DateSelection = form.sel
    return {
      year: s.year, month: s.month, day: s.day,
      hour: s.hour, minute: s.minute, second: 0,
      gender: form.gender,
      isLunar: s.inputMode === MODE_LUNAR,
      isLeapMonth: s.inputMode === MODE_LUNAR ? s.isLeap : false,
      isSpecMonth: s.inputMode === MODE_LUNAR ? s.isSpec : false,
      name: form.name.length > 0 ? form.name : '匿名',
      astEnabled: form.astEnabled,
      longitude: form.longitude,
      latitude: form.latitude,
      lifa: form.lifa
    } as BaziParams
  }
}
