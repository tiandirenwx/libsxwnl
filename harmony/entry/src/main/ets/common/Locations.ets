import { GeoCity, NativeBridge } from '../services/NativeBridge'

// ═══════════════════════════════════════════════════════════════════
//  Locations — 地理位置工具函数 + 默认值
//
//  所有目录数据 (32 省 / 4000+ 市县 / 国际时区) 都驻留在 libsxwnl C++ 层
//  (src/geo.cpp 的 GeoPostion 单例), 本文件不存任何硬编码城市表.
//  通过 NativeBridge.geoXxx() 一次性懒加载即可.
// ═══════════════════════════════════════════════════════════════════

export type Location = GeoCity

// 默认地点 (北京天安门, 由 libsxwnl getDefaultGeoPos 提供).
// NAPI 不可用时也给出兜底值, 避免 UI 显示空数据.
export function getDefaultLocation(): Location {
  const fromNative = NativeBridge.geoDefault()
  if (fromNative !== null) return fromNative
  return {
    province: '北京市', district: '天安门',
    longitude: 116.3833, latitude: 39.9, timezone: 8,
  }
}

// 显示标签 "北京市·天安门"
export function locationLabel(loc: Location): string {
  if (loc.province && loc.district && loc.province !== loc.district) {
    return `${loc.province}·${loc.district}`
  }
  return loc.district || loc.province
}

// 显示时区, e.g. 8 → "UTC+8", -5 → "UTC-5", 5.5 → "UTC+5:30"
export function formatTimezone(tz: number): string {
  const sign = tz >= 0 ? '+' : '-'
  const abs  = Math.abs(tz)
  const h    = Math.floor(abs)
  const m    = Math.round((abs - h) * 60)
  if (m === 0) return `UTC${sign}${h}`
  return `UTC${sign}${h}:${m.toString().padStart(2, '0')}`
}

// 显示经纬度 "116.4°E 39.9°N"
export function formatLonLat(lon: number, lat: number): string {
  const lonDir = lon >= 0 ? 'E' : 'W'
  const latDir = lat >= 0 ? 'N' : 'S'
  return `${Math.abs(lon).toFixed(1)}°${lonDir} ${Math.abs(lat).toFixed(1)}°${latDir}`
}
