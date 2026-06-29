package com.sxwnl.calendar.util

import android.content.Context
import android.content.Intent
import android.graphics.Bitmap
import android.graphics.Canvas
import android.net.Uri
import android.view.View
import androidx.core.content.FileProvider
import java.io.File
import java.io.FileOutputStream
import java.text.SimpleDateFormat
import java.util.*

/**
 * EclipseShareUtil — 截图分享 + ICS 日历导出
 *
 * 与 iOS [EclipseShareUtil.swift] / 鸿蒙端对齐.
 */
object EclipseShareUtil {

    // ─── 截图: View → Bitmap → 临时 PNG → 分享 ───────────────

    /** 将任意 [View] 渲染为 Bitmap. view 未布局时返回 1×1 占位, 调用方负责处理. */
    fun captureView(view: View): Bitmap {
        val w = view.width.coerceAtLeast(1)
        val h = view.height.coerceAtLeast(1)
        val bmp = Bitmap.createBitmap(w, h, Bitmap.Config.ARGB_8888)
        val c = Canvas(bmp)
        view.draw(c)
        return bmp
    }

    /** 保存 [Bitmap] 到 cache 目录, 返回 File. 落盘后调用 [Bitmap.recycle] 释放. */
    fun saveBitmap(context: Context, bmp: Bitmap, fileName: String): File {
        val dir = File(context.cacheDir, "shares").apply { mkdirs() }
        val file = File(dir, fileName)
        try {
            FileOutputStream(file).use { os ->
                bmp.compress(Bitmap.CompressFormat.PNG, 92, os)
            }
        } finally {
            if (!bmp.isRecycled) bmp.recycle()
        }
        return file
    }

    /** 调起系统分享 (图片). */
    fun shareImage(context: Context, file: File) {
        val uri: Uri = FileProvider.getUriForFile(
            context, "${context.packageName}.fileprovider", file)
        val intent = Intent(Intent.ACTION_SEND).apply {
            type = "image/png"
            putExtra(Intent.EXTRA_STREAM, uri)
            addFlags(Intent.FLAG_GRANT_READ_URI_PERMISSION)
        }
        context.startActivity(Intent.createChooser(intent, "分享日月食截图")
            .addFlags(Intent.FLAG_ACTIVITY_NEW_TASK))
    }

    /** 调起系统分享 (ICS 日历文件). */
    fun shareFile(context: Context, file: File, mime: String = "text/calendar") {
        val uri: Uri = FileProvider.getUriForFile(
            context, "${context.packageName}.fileprovider", file)
        val intent = Intent(Intent.ACTION_SEND).apply {
            type = mime
            putExtra(Intent.EXTRA_STREAM, uri)
            addFlags(Intent.FLAG_GRANT_READ_URI_PERMISSION)
        }
        context.startActivity(Intent.createChooser(intent, "导出日程")
            .addFlags(Intent.FLAG_ACTIVITY_NEW_TASK))
    }

    // ─── ICS 文件构造 ──────────────────────────────────────

    /**
     * 构造 ICS (单事件).
     *
     * @param startJd 起始时刻 (TD, JD)
     * @param endJd 结束时刻 (TD, JD); 若 ≤ startJd, 自动 +1h
     */
    fun buildICS(
        uid: String,
        summary: String,
        description: String,
        startJd: Double,
        endJd: Double,
        location: String = ""
    ): String {
        val start = jdToIcsUtc(startJd)
        val end = jdToIcsUtc(if (endJd > startJd) endJd else startJd + 1.0 / 24)
        val now = SimpleDateFormat("yyyyMMdd'T'HHmmss'Z'", Locale.US).apply {
            timeZone = TimeZone.getTimeZone("UTC")
        }.format(Date())

        val sb = StringBuilder()
        sb.append("BEGIN:VCALENDAR\r\n")
        sb.append("VERSION:2.0\r\n")
        sb.append("PRODID:-//libsxwnl//Android Eclipse//ZH\r\n")
        sb.append("CALSCALE:GREGORIAN\r\n")
        sb.append("BEGIN:VEVENT\r\n")
        sb.append("UID:$uid\r\n")
        sb.append("DTSTAMP:$now\r\n")
        sb.append("DTSTART:$start\r\n")
        sb.append("DTEND:$end\r\n")
        sb.append("SUMMARY:${escape(summary)}\r\n")
        sb.append("DESCRIPTION:${escape(description)}\r\n")
        if (location.isNotEmpty()) sb.append("LOCATION:${escape(location)}\r\n")
        sb.append("END:VEVENT\r\n")
        sb.append("END:VCALENDAR\r\n")
        return sb.toString()
    }

    fun writeICS(context: Context, ics: String, fileName: String): File {
        val dir = File(context.cacheDir, "shares").apply { mkdirs() }
        val file = File(dir, fileName)
        file.writeText(ics, Charsets.UTF_8)
        return file
    }

    // ─── 内部 ─────────────────────────────────────────────

    /** 力学时 JD → "YYYYMMDDTHHMMSSZ" (作为 UTC 输出). */
    private fun jdToIcsUtc(jd: Double): String {
        val s = EclipseUtil.jdToDateTime(jd) // "yyyy-MM-dd HH:mm:ss"
        if (s.length < 19) return "19700101T000000Z"
        return "${s.substring(0, 4)}${s.substring(5, 7)}${s.substring(8, 10)}T" +
               "${s.substring(11, 13)}${s.substring(14, 16)}${s.substring(17, 19)}Z"
    }

    private fun escape(s: String) = s
        .replace("\\", "\\\\")
        .replace(",", "\\,")
        .replace(";", "\\;")
        .replace("\n", "\\n")
}
