import Foundation
import SwiftUI
import UIKit

// ════════════════════════════════════════════════════════════════
//  EclipseShareUtil — 截图分享 + ICS 日历导出
//  与 Android `EclipseShareUtil.kt` / 鸿蒙端对齐
// ════════════════════════════════════════════════════════════════

enum EclipseShareUtil {

    /// 渲染任意 SwiftUI View 到 PNG, 并保存到临时目录, 返回 URL
    @MainActor
    static func snapshot<Content: View>(of view: Content,
                                        size: CGSize? = nil,
                                        fileName: String = "eclipse.png") -> URL? {
        let controller = UIHostingController(rootView: view.frame(maxWidth: 600))
        let target = size ?? controller.sizeThatFits(in: CGSize(width: 600, height: 9999))
        controller.view.bounds = CGRect(origin: .zero, size: target)
        controller.view.backgroundColor = .clear

        let renderer = UIGraphicsImageRenderer(size: target)
        let img = renderer.image { _ in
            controller.view.drawHierarchy(in: controller.view.bounds, afterScreenUpdates: true)
        }
        guard let data = img.pngData() else { return nil }
        let url = FileManager.default.temporaryDirectory.appendingPathComponent(fileName)
        try? data.write(to: url, options: .atomic)
        return url
    }

    /// 构造 ICS 文件 (单事件)
    /// - startJd / endJd: 力学时 JD
    /// - tzHours: 输出 UTC 时间, 由 startJd/endJd 内部已经是 TD
    static func buildICS(uid: String,
                         summary: String,
                         description: String,
                         startJd: Double,
                         endJd: Double,
                         location: String = "") -> String {
        let startTd = jdToUtcCalendar(startJd)
        let endTd   = jdToUtcCalendar(endJd > startJd ? endJd : startJd + 1.0/24)
        let now     = jdToUtcCalendar(jdFromDate(Date()))

        let lines: [String] = [
            "BEGIN:VCALENDAR",
            "VERSION:2.0",
            "PRODID:-//libsxwnl//iOS Eclipse//ZH",
            "CALSCALE:GREGORIAN",
            "BEGIN:VEVENT",
            "UID:\(uid)",
            "DTSTAMP:\(now)",
            "DTSTART:\(startTd)",
            "DTEND:\(endTd)",
            "SUMMARY:\(escape(summary))",
            "DESCRIPTION:\(escape(description))",
            (location.isEmpty ? nil : "LOCATION:\(escape(location))") ?? "",
            "END:VEVENT",
            "END:VCALENDAR"
        ]
        return lines.filter { !$0.isEmpty }.joined(separator: "\r\n")
    }

    /// 保存 ICS 文本到临时目录, 返回 URL
    static func writeICS(_ ics: String, fileName: String) -> URL? {
        let url = FileManager.default.temporaryDirectory.appendingPathComponent(fileName)
        do {
            try ics.write(to: url, atomically: true, encoding: .utf8)
            return url
        } catch { return nil }
    }

    // ─── helpers ──────────────────────────────────────────

    private static func jdFromDate(_ date: Date) -> Double {
        // Unix epoch JD = 2440587.5
        return date.timeIntervalSince1970 / 86400.0 + 2440587.5
    }

    /// JD (assumed TD ~= UTC for ICS purposes) → "YYYYMMDDTHHMMSSZ"
    private static func jdToUtcCalendar(_ jd: Double) -> String {
        let s = EclipseUtil.jdToDateTime(jd)
        // s = "YYYY-MM-DD HH:MM:SS"
        guard s.count >= 19 else { return "19700101T000000Z" }
        let y = s.prefix(4), m = s.dropFirst(5).prefix(2), d = s.dropFirst(8).prefix(2)
        let hh = s.dropFirst(11).prefix(2), mm = s.dropFirst(14).prefix(2)
        let ss = s.dropFirst(17).prefix(2)
        return "\(y)\(m)\(d)T\(hh)\(mm)\(ss)Z"
    }

    private static func escape(_ s: String) -> String {
        s.replacingOccurrences(of: "\\", with: "\\\\")
            .replacingOccurrences(of: ",", with: "\\,")
            .replacingOccurrences(of: ";", with: "\\;")
            .replacingOccurrences(of: "\n", with: "\\n")
    }
}

// MARK: - Share Sheet

struct ShareSheet: UIViewControllerRepresentable {
    let items: [Any]
    func makeUIViewController(context: Context) -> UIActivityViewController {
        UIActivityViewController(activityItems: items, applicationActivities: nil)
    }
    func updateUIViewController(_ uiViewController: UIActivityViewController, context: Context) {}
}
