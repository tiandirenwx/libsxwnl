import fs from '@ohos.file.fs'
import { EclipseUtil } from './EclipseUtil'

// ════════════════════════════════════════════════════════════════
//  EclipseShareUtil — ICS 日程导出 (ArkTS)
//
//  ICS 完全本地构造, 写入沙箱缓存目录.
//  与 Android/iOS 端 EclipseShareUtil 对齐 (截图分享暂用系统拍照/录屏).
// ════════════════════════════════════════════════════════════════

export class EclipseShareUtil {
  static buildICS(
    uid: string,
    summary: string,
    description: string,
    startJd: number,
    endJd: number,
    location: string = ''
  ): string {
    const start = EclipseShareUtil.jdToIcsUtc(startJd)
    const end = EclipseShareUtil.jdToIcsUtc(endJd > startJd ? endJd : startJd + 1 / 24)
    const now = EclipseShareUtil.dateToIcsUtc(new Date())
    const lines: string[] = [
      'BEGIN:VCALENDAR',
      'VERSION:2.0',
      'PRODID:-//libsxwnl//Harmony Eclipse//ZH',
      'CALSCALE:GREGORIAN',
      'BEGIN:VEVENT',
      'UID:' + uid,
      'DTSTAMP:' + now,
      'DTSTART:' + start,
      'DTEND:' + end,
      'SUMMARY:' + EclipseShareUtil.escape(summary),
      'DESCRIPTION:' + EclipseShareUtil.escape(description)
    ]
    if (location && location.length > 0) lines.push('LOCATION:' + EclipseShareUtil.escape(location))
    lines.push('END:VEVENT')
    lines.push('END:VCALENDAR')
    return lines.join('\r\n')
  }

  /** 写入 ICS 到沙箱缓存目录, 返回文件路径; 失败返回空串 */
  static writeICS(cacheDir: string, ics: string, fileName: string): string {
    const path = cacheDir + '/' + fileName
    try {
      const file = fs.openSync(path, fs.OpenMode.READ_WRITE | fs.OpenMode.CREATE | fs.OpenMode.TRUNC)
      try {
        fs.writeSync(file.fd, ics)
      } finally {
        fs.closeSync(file)
      }
      return path
    } catch (e) {
      return ''
    }
  }

  private static jdToIcsUtc(jd: number): string {
    const s = EclipseUtil.jdToDateTime(jd)
    if (s.length < 19) return '19700101T000000Z'
    return s.substring(0, 4) + s.substring(5, 7) + s.substring(8, 10) + 'T' +
           s.substring(11, 13) + s.substring(14, 16) + s.substring(17, 19) + 'Z'
  }

  private static dateToIcsUtc(d: Date): string {
    const y = d.getUTCFullYear(), m = d.getUTCMonth() + 1, day = d.getUTCDate()
    const hh = d.getUTCHours(), mm = d.getUTCMinutes(), ss = d.getUTCSeconds()
    const p = (v: number, n: number) => {
      const s = v.toString()
      return s.length >= n ? s : '0'.repeat(n - s.length) + s
    }
    return p(y, 4) + p(m, 2) + p(day, 2) + 'T' + p(hh, 2) + p(mm, 2) + p(ss, 2) + 'Z'
  }

  private static escape(s: string): string {
    return s.replace(/\\/g, '\\\\').replace(/,/g, '\\,')
            .replace(/;/g, '\\;').replace(/\n/g, '\\n')
  }
}
