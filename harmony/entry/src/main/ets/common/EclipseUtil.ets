// ════════════════════════════════════════════════════════════════
//  EclipseUtil — 日月食类型分类 + JD 格式化
//
//  与 Android / iOS 端 EclipseUtil 行为对齐
// ════════════════════════════════════════════════════════════════

export class EclipseUtil {

  /** 日食类型徽章颜色 */
  static solarBadgeColor(type: string): string {
    if (!type) return '#9E9E9E'
    const c = type.charAt(0)
    if (c === 'T') return '#D32F2F'
    if (c === 'A') return '#F57C00'
    if (c === 'H') return '#7B1FA2'
    if (c === 'P') return '#616161'
    return '#9E9E9E'
  }

  static lunarBadgeColor(type: string): string {
    if (type === 'T') return '#D32F2F'
    if (type === 'P') return '#616161'
    if (type === 'B') return '#1976D2'
    return '#9E9E9E'
  }

  /** 力学时 JD → "YYYY-MM-DD HH:MM:SS" */
  static jdToDateTime(jd: number): string {
    if (!jd || jd <= 0) return '--'
    const jdAdj = jd + 0.5
    const Z = Math.floor(jdAdj)
    const F = jdAdj - Z
    let A = Z
    if (Z >= 2_299_161) {
      const alpha = Math.floor((Z - 1_867_216.25) / 36_524.25)
      A = Z + 1 + alpha - Math.floor(alpha / 4)
    }
    const B = A + 1524
    const C = Math.floor((B - 122.1) / 365.25)
    const D = Math.floor(365.25 * C)
    const E = Math.floor((B - D) / 30.6001)
    const day = (B - D - Math.floor(30.6001 * E)) + F
    const mon = E < 14 ? E - 1 : E - 13
    const yr = mon > 2 ? C - 4716 : C - 4715
    const dInt = Math.floor(day)
    let dFrac = day - dInt
    let hh = Math.floor(dFrac * 24); dFrac = dFrac * 24 - hh
    let mm = Math.floor(dFrac * 60); dFrac = dFrac * 60 - mm
    let ss = Math.floor(dFrac * 60 + 0.5)
    if (ss >= 60) { ss = 0; mm += 1 }
    if (mm >= 60) { mm = 0; hh += 1 }
    return EclipseUtil.pad(yr, 4) + '-' + EclipseUtil.pad(mon, 2) + '-' + EclipseUtil.pad(dInt, 2)
      + ' ' + EclipseUtil.pad(hh, 2) + ':' + EclipseUtil.pad(mm, 2) + ':' + EclipseUtil.pad(ss, 2)
  }

  static jdToTime(jd: number): string {
    if (!jd || jd <= 0) return '--'
    const s = EclipseUtil.jdToDateTime(jd)
    return s.length >= 19 ? s.substring(11) : '--'
  }

  static magnitudeProgress(mag: number): number {
    return Math.max(0, Math.min(1.5, mag)) / 1.5
  }

  static formatDuration(sec: number): string {
    if (sec <= 0) return '—'
    const s = Math.floor(sec)
    const h = Math.floor(s / 3600)
    const m = Math.floor((s % 3600) / 60)
    const ss = s % 60
    return h > 0
      ? h + 'h' + EclipseUtil.pad(m, 2) + 'm' + EclipseUtil.pad(ss, 2) + 's'
      : m + 'm' + EclipseUtil.pad(ss, 2) + 's'
  }

  private static pad(v: number, n: number): string {
    const s = Math.abs(v).toString()
    const sign = v < 0 ? '-' : ''
    return sign + (s.length >= n ? s : '0'.repeat(n - s.length) + s)
  }

  // ─── ΔT (TD − UT1) Espenak/Meeus polynomial ──────────────
  static deltaT(year: number): number {
    const y = year
    if (y < -500 || y > 2150) {
      const u = (y - 1820) / 100
      return -20 + 32 * u * u
    }
    if (y < 500) {
      const u = y / 100
      return 10583.6 + u * (-1014.41 + u * (33.78311 + u * (-5.952053 + u * (-0.1798452 + u * (0.022174192 + u * 0.0090316521)))))
    }
    if (y < 1600) {
      const u = (y - 1000) / 100
      return 1574.2 + u * (-556.01 + u * (71.23472 + u * (0.319781 + u * (-0.8503463 + u * (-0.005050998 + u * 0.0083572073)))))
    }
    if (y < 1700) { const t = y - 1600; return 120 + t * (-0.9808 + t * (-0.01532 + t / 7129)) }
    if (y < 1800) { const t = y - 1700; return 8.83 + t * (0.1603 + t * (-0.0059285 + t * (0.00013336 - t / 1174000))) }
    if (y < 1860) { const t = y - 1800; return 13.72 + t * (-0.332447 + t * (0.0068612 + t * (0.0041116 + t * (-0.00037436 + t * (0.0000121272 + t * (-0.0000001699 + t * 0.000000000875)))))) }
    if (y < 1900) { const t = y - 1860; return 7.62 + t * (0.5737 + t * (-0.251754 + t * (0.01680668 + t * (-0.0004473624 + t / 233174)))) }
    if (y < 1920) { const t = y - 1900; return -2.79 + t * (1.494119 + t * (-0.0598939 + t * (0.0061966 - t * 0.000197))) }
    if (y < 1941) { const t = y - 1920; return 21.20 + t * (0.84493 + t * (-0.076100 + t * 0.0020936)) }
    if (y < 1961) { const t = y - 1950; return 29.07 + t * (0.407 + t * (-1 / 233 + t / 2547)) }
    if (y < 1986) { const t = y - 1975; return 45.45 + t * (1.067 + t * (-1 / 260 - t / 718)) }
    if (y < 2005) { const t = y - 2000; return 63.86 + t * (0.3345 + t * (-0.060374 + t * (0.0017275 + t * (0.000651814 + t * 0.00002373599)))) }
    if (y < 2050) { const t = y - 2000; return 62.92 + t * (0.32217 + t * 0.005589) }
    if (y < 2150) return -20 + 32 * Math.pow((y - 1820) / 100, 2) - 0.5628 * (2150 - y)
    const u = (y - 1820) / 100
    return -20 + 32 * u * u
  }

  /** 力学时 JD → 本地时 "HH:MM:SS" (含 ΔT + 时区) */
  static jdTdToLocal(jd: number, tzHours: number, withDate: boolean = false): string {
    if (!jd || jd <= 0) return '--'
    const s0 = EclipseUtil.jdToDateTime(jd)
    const yr = parseInt(s0.substring(0, 4)) || 2000
    const dt = EclipseUtil.deltaT(yr) / 86400
    const ut = jd - dt + tzHours / 24
    const s = EclipseUtil.jdToDateTime(ut)
    return withDate ? s : (s.length >= 19 ? s.substring(11) : '--')
  }

  static tzLabel(tzHours: number): string {
    const sign = tzHours >= 0 ? '+' : '-'
    const abs = Math.abs(tzHours)
    const h = Math.floor(abs)
    const m = Math.floor((abs - h) * 60)
    return 'UTC' + sign + EclipseUtil.padStatic(h, 2) + ':' + EclipseUtil.padStatic(m, 2)
  }

  /** 经度(弧度) → 屏幕 X 比例 0..1 */
  static projectLonX(lonRad: number): number {
    return (lonRad + Math.PI) / (2 * Math.PI)
  }
  /** 纬度(弧度) → 屏幕 Y 比例 0..1, 北上 */
  static projectLatY(latRad: number): number {
    return 0.5 - latRad / Math.PI
  }

  static lerp(a: number, b: number, t: number): number {
    return a + (b - a) * t
  }

  static solarEmoji(type: string): string {
    if (!type) return '🌗'
    const c = type.charAt(0)
    if (c === 'T') return '🌑'
    if (c === 'A') return '🌒'
    if (c === 'H') return '🌓'
    return '🌗'
  }
  static lunarEmoji(type: string): string {
    if (type === 'T') return '🌕'
    if (type === 'P') return '🌖'
    if (type === 'B') return '🌗'
    return '🌘'
  }

  static lonLabel(lonDeg: number): string {
    const abs = Math.abs(lonDeg)
    return abs.toFixed(2) + '°' + (lonDeg >= 0 ? 'E' : 'W')
  }
  static latLabel(latDeg: number): string {
    const abs = Math.abs(latDeg)
    return abs.toFixed(2) + '°' + (latDeg >= 0 ? 'N' : 'S')
  }

  private static padStatic(v: number, n: number): string {
    const s = Math.abs(v).toString()
    const sign = v < 0 ? '-' : ''
    return sign + (s.length >= n ? s : '0'.repeat(n - s.length) + s)
  }
}
