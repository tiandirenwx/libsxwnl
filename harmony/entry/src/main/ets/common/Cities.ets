// 出生地经纬度数据 — 薄包装 NativeBridge.geoXxx()
//
// 数据来源: libsxwnl C++ 层的 GeoPostion 单例 (src/geo.cpp 中的 JWv).
// 本文件不再持有任何硬编码城市表; 上层 (BaziPage / CityPickerDialog /
// EclipsePage / EclipseComponents) 可继续使用 Cities.REGIONS / DEFAULT /
// PRESET 等既有静态 API, 与跨端 native 库的目录数据自动保持同步.
//
// 与 Android `Cities.kt` / iOS `Cities.swift` 字段结构对齐.

import { GeoCity, GeoProvince, NativeBridge } from '../services/NativeBridge'

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

function toObserverCity(c: GeoCity): ObserverCity {
  return {
    name: c.district,
    longitude: c.longitude,
    latitude: c.latitude,
    tzHours: c.timezone,
  }
}

function loadRegions(): Region[] {
  const out: Region[] = []
  const provinces: GeoProvince[] = NativeBridge.geoListProvinces()
  for (let i = 0; i < provinces.length; i++) {
    const pname: string = provinces[i].province
    const raw: GeoCity[] = NativeBridge.geoListCities(pname)
    if (raw.length === 0) continue
    const cities: ObserverCity[] = []
    for (let k = 0; k < raw.length; k++) cities.push(toObserverCity(raw[k]))
    out.push({ name: pname, cities: cities })
  }
  return out
}

// 兜底默认 (NAPI 完全不可用时, 防止下游因空数组解引用崩溃)
const FALLBACK: ObserverCity = {
  name: '天安门', longitude: 116.3833, latitude: 39.9, tzHours: 8,
}

// ── 静态目录 (整个 App 共享一次性加载结果) ──────────────
// libsxwnl 是 in-process native, 一次性拉 32 省 × ~100 城市约 50ms,
// 之后的所有访问都是普通对象索引, 不再触发 NAPI.
const REGIONS: Region[] = loadRegions()
const DEFAULT: ObserverCity =
  REGIONS.length > 0 && REGIONS[0].cities.length > 0 ? REGIONS[0].cities[0] : FALLBACK
const PRESET: ObserverCity[] =
  REGIONS.length > 0 ? REGIONS[0].cities.slice(0, 30) : [DEFAULT]

export class Cities {
  static readonly REGIONS: Region[] = REGIONS
  static readonly DEFAULT: ObserverCity = DEFAULT
  // 旧版兼容: 扁平预设 (供日食观测点快速选择条等模块使用), 取首批 30 个 JWv 内置项
  static readonly PRESET: ObserverCity[] = PRESET
}
