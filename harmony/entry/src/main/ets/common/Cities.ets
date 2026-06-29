// 出生地经纬度数据 (源自寿星天文历 JWv 压缩字符串)
//
// 数据规模: 32 个省级行政区, 每个 30-200 个地级市/区县, 合计 ~3500 条
// 第一级 = 省/直辖市/特别行政区, 第二级 = 城市/区县
// 八字计算只取经纬度, 时区一律按北京时间(UT+8) 处理
//
// 与 Android `Cities.kt` / iOS `Cities.swift` 字段结构对齐

import { JWV_RAW } from './JWVRaw'

export interface ObserverCity {
  name: string
  longitude: number   // 度, 东正
  latitude: number    // 度, 北正
  // tzHours 仅供日食/RTS 等需要观测点时区的模块使用; 八字模块固定按北京时间(8)处理
  tzHours: number
}

export interface Region {
  name: string        // 省 / 直辖市 / 特别行政区
  cities: ObserverCity[]
}

interface JWPair {
  lon: number
  lat: number
}

// 单字符 base62 解码: '0'-'9' = 0-9, 'A'-'Z' = 10-35, 'a'-'z' = 36-61
function decodeCh(c: number): number {
  if (c > 96) return c - (97 - 36)
  if (c > 64) return c - (65 - 10)
  return c - 48
}

// 4 字符压缩串 → 经纬度 (度): [纬度度, 纬度分, 经度度-73, 经度分]
function decodeJW(code: string): JWPair {
  const c0: number = decodeCh(code.charCodeAt(0))
  const c1: number = decodeCh(code.charCodeAt(1))
  const c2: number = decodeCh(code.charCodeAt(2))
  const c3: number = decodeCh(code.charCodeAt(3))
  const pair: JWPair = {
    lon: c2 + c3 / 60.0 + 73.0,
    lat: c0 + c1 / 60.0
  }
  return pair
}

function buildRegions(): Region[] {
  const out: Region[] = []
  for (const raw of JWV_RAW) {
    const parts: string[] = raw.split(' ')
    if (parts.length < 2) continue
    const regionName: string = parts[0]
    const cities: ObserverCity[] = []
    for (let i = 1; i < parts.length; i++) {
      const entry: string = parts[i]
      if (entry.length < 5) continue
      const pair: JWPair = decodeJW(entry.substring(0, 4))
      const city: ObserverCity = {
        name: entry.substring(4),
        longitude: pair.lon,
        latitude: pair.lat,
        tzHours: 8
      }
      cities.push(city)
    }
    if (cities.length > 0) {
      const region: Region = { name: regionName, cities: cities }
      out.push(region)
    }
  }
  return out
}

export class Cities {
  static readonly REGIONS: Region[] = buildRegions()

  // 默认北京天安门
  static readonly DEFAULT: ObserverCity = Cities.REGIONS[0].cities[0]

  // 旧版兼容: 扁平预设(供日食观测点快速选择条等模块使用), 取首批 30 个 JWv 内置项
  static readonly PRESET: ObserverCity[] = Cities.REGIONS[0].cities.slice(0, 30)
}
