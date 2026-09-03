package com.sxwnl.calendar.ui.components

import android.content.ContentValues
import android.content.Context
import android.graphics.Bitmap
import android.graphics.Picture
import android.os.Build
import android.os.Environment
import android.provider.MediaStore
import android.widget.Toast
import androidx.annotation.DrawableRes
import androidx.compose.foundation.Image
import androidx.compose.foundation.background
import androidx.compose.foundation.border
import androidx.compose.foundation.clickable
import androidx.compose.foundation.horizontalScroll
import androidx.compose.foundation.layout.*
import androidx.compose.foundation.rememberScrollState
import androidx.compose.foundation.shape.RoundedCornerShape
import androidx.compose.foundation.verticalScroll
import androidx.compose.material3.*
import androidx.compose.runtime.*
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.draw.clip
import androidx.compose.ui.draw.drawWithCache
import androidx.compose.ui.graphics.Canvas
import androidx.compose.ui.graphics.Color
import androidx.compose.ui.graphics.drawscope.draw
import androidx.compose.ui.graphics.drawscope.drawIntoCanvas
import androidx.compose.ui.graphics.nativeCanvas
import androidx.compose.ui.layout.ContentScale
import androidx.compose.ui.platform.LocalContext
import androidx.compose.ui.res.painterResource
import androidx.compose.ui.text.font.Font
import androidx.compose.ui.text.font.FontFamily
import androidx.compose.ui.text.font.FontWeight
import androidx.compose.ui.text.style.TextAlign
import androidx.compose.ui.unit.TextUnit
import androidx.compose.ui.unit.dp
import androidx.compose.ui.unit.sp
import com.sxwnl.calendar.R
import com.sxwnl.calendar.data.BaziColumn
import com.sxwnl.calendar.data.BaziResult
import com.sxwnl.calendar.data.LiuNianItem
import com.sxwnl.calendar.ui.theme.*
import com.sxwnl.calendar.util.BaziCalc
import com.sxwnl.calendar.util.EclipseShareUtil
import kotlinx.coroutines.CoroutineScope
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.launch
import kotlinx.coroutines.withContext
import java.io.File
import java.io.FileOutputStream

// ════════════════════════════════════════════════════════════════
//  BaziResultView — 与鸿蒙 components/BaziResultView.ets 对齐
//  简洁版(宣纸命书) + 专业版(命盘表)
// ════════════════════════════════════════════════════════════════

data class BaziResultArg(
    val result: BaziResult,
    val birthYear: Int,
    val astEnabled: Boolean,
    val longitude: Double,
    val latitude: Double,
    val lifaLabel: String
)

private val ElemColors = listOf(
    Color(0xFF2E7D32), // 木
    Color(0xFFD32F2F), // 火
    Color(0xFF8D6E40), // 土
    Color(0xFFC9A227), // 金
    Color(0xFF1565C0)  // 水
)

// 宣纸命书色板 (与鸿蒙 BaziResultView.ets 对齐)
private val InkColor = Color(0xFF3A2A1C)
private val InkSoft  = Color(0xFF6B5640)
private val RedColor = Color(0xFFC0392B)
private val RedBg    = Color(0xFFF6DAD2)
private val GoldLine = Color(0xFFC9A86A)
// 兜底纸色: 加载宣纸图前/失败时使用
private val PaperBg  = Color(0xFFFAF3E0)

// 12 生肖水印 (索引 0=子鼠..11=亥猪), 与鸿蒙 zodiacRes() 对齐
@DrawableRes
private val ZodiacDrawables: IntArray = intArrayOf(
    R.drawable.bz_rat,    R.drawable.bz_cow,    R.drawable.bz_tiger,
    R.drawable.bz_rabbit, R.drawable.bz_dragon, R.drawable.bz_snake,
    R.drawable.bz_horse,  R.drawable.bz_goat,   R.drawable.bz_monkey,
    R.drawable.bz_hen,    R.drawable.bz_dog,    R.drawable.bz_pig
)

private const val ZHI_ORDER = "子丑寅卯辰巳午未申酉戌亥"

@DrawableRes
private fun zodiacDrawable(yearZhi: String): Int {
    val idx = ZHI_ORDER.indexOf(yearZhi.firstOrNull() ?: '子').coerceAtLeast(0)
    return ZodiacDrawables[idx]
}

// 宣纸命书专用毛笔字体 (与鸿蒙 BaziPage.ets 注册的 WenYue.otf 同源)
private val WenYueFont = FontFamily(Font(R.font.wen_yue))

private const val LabelW = 40
private const val ColW = 54
private val ColTitles = listOf("年柱", "月柱", "日柱", "时柱", "大运", "流年")

@OptIn(ExperimentalMaterial3Api::class)
@Composable
fun BaziResultView(
    arg: BaziResultArg,
    onDismiss: () -> Unit
) {
    val sheetState = rememberModalBottomSheetState(skipPartiallyExpanded = true)
    ModalBottomSheet(
        onDismissRequest = onDismiss,
        sheetState = sheetState,
        containerColor = Background,
        shape = RoundedCornerShape(topStart = Dimens.radiusLg, topEnd = Dimens.radiusLg)
    ) {
        ResultBody(arg, onDismiss)
    }
}

@Composable
private fun ResultBody(arg: BaziResultArg, onDismiss: () -> Unit) {
    var proMode by remember { mutableStateOf(false) }
    val isFemale = arg.result.gender.contains("女")
    val ctx = LocalContext.current
    val scope = rememberCoroutineScope()
    // Compose 1.6 兼容方案: 用 android.graphics.Picture 录制 Composable.
    // 在待截图区域上挂 drawWithCache, 截图时直接 picture → Bitmap。
    //   简洁版录制"命书卡片", 专业版录制"命盘表整体"。
    val simpleCardPicture = remember { Picture() }
    val proCardPicture = remember { Picture() }

    Column(Modifier.fillMaxSize()) {
        // 顶部条
        Row(
            Modifier
                .fillMaxWidth()
                .padding(Dimens.paddingLg),
            verticalAlignment = Alignment.CenterVertically
        ) {
            Text(
                "关闭", color = TextSecondary,
                modifier = Modifier.clickable(onClick = onDismiss)
            )
            Spacer(Modifier.weight(1f))
            Text(
                if (isFemale) "排盘 · 坤造" else "排盘 · 乾造",
                fontSize = Dimens.fontSubtitle, fontWeight = FontWeight.Bold,
                color = Primary
            )
            Spacer(Modifier.weight(1f))
            Box(Modifier.width(28.dp))
        }
        HorizontalDivider(color = DividerColor)

        Column(
            Modifier
                .fillMaxSize()
                .verticalScroll(rememberScrollState())
                .padding(bottom = 24.dp)
        ) {
            // 版本切换
            Row(
                Modifier
                    .padding(top = Dimens.paddingMd)
                    .align(Alignment.CenterHorizontally)
                    .clip(RoundedCornerShape(17.dp))
                    .border(1.dp, DividerColor, RoundedCornerShape(17.dp))
            ) {
                ToggleChip(label = "简洁版", active = !proMode,
                    cornerStart = true) { proMode = false }
                ToggleChip(label = "专业版", active = proMode,
                    cornerEnd = true) { proMode = true }
            }

            if (proMode) {
                // 把专业版命盘整体录制到 Picture 供"保存/分享"截图 (含背景)
                Column(
                    Modifier
                        .fillMaxWidth()
                        .recordToPicture(proCardPicture)
                        .background(Background)
                ) {
                    InfoBar(arg, isFemale)
                    PaiPanTable(arg)
                    SummaryBar(arg)
                    LiuNianCard(arg)
                }
                SaveShareRow(
                    onSave = { saveBaziImage(ctx, scope, proCardPicture) },
                    onShare = { shareBaziImage(ctx, scope, proCardPicture) }
                )
            } else {
                SimpleCard(
                    arg = arg,
                    isFemale = isFemale,
                    picture = simpleCardPicture,
                    onSave = { saveBaziImage(ctx, scope, simpleCardPicture) },
                    onShare = { shareBaziImage(ctx, scope, simpleCardPicture) }
                )
            }
        }
    }
}

@Composable
private fun ToggleChip(
    label: String, active: Boolean,
    cornerStart: Boolean = false, cornerEnd: Boolean = false,
    onClick: () -> Unit
) {
    val shape = if (cornerStart)
        RoundedCornerShape(topStart = 17.dp, bottomStart = 17.dp)
    else if (cornerEnd)
        RoundedCornerShape(topEnd = 17.dp, bottomEnd = 17.dp)
    else RoundedCornerShape(0.dp)
    Box(
        Modifier
            .clip(shape)
            .background(if (active) Primary else Surface)
            .clickable(onClick = onClick)
            .width(90.dp)
            .height(34.dp),
        contentAlignment = Alignment.Center
    ) {
        Text(
            label, fontSize = Dimens.fontCaption,
            fontWeight = if (active) FontWeight.Bold else FontWeight.Normal,
            color = if (active) OnPrimary else TextSecondary
        )
    }
}

// ════════════════════ 简洁版(宣纸命书) ════════════════════

@Composable
private fun SimpleCard(
    arg: BaziResultArg, isFemale: Boolean,
    picture: Picture,
    onSave: () -> Unit,
    onShare: () -> Unit
) {
    Column(
        Modifier.fillMaxWidth(),
        horizontalAlignment = Alignment.CenterHorizontally
    ) {
        SimpleReport(
            arg = arg,
            isFemale = isFemale,
            modifier = Modifier.recordToPicture(picture)
        )
        SaveShareRow(onSave = onSave, onShare = onShare)
    }
}

// 保存为图片 + 分享到 (简洁版/专业版共用): 左"保存到相册", 右"分享到其它 app"
@Composable
private fun SaveShareRow(onSave: () -> Unit, onShare: () -> Unit) {
    Row(
        Modifier
            .fillMaxWidth()
            .padding(horizontal = 24.dp, vertical = Dimens.paddingSm),
        horizontalArrangement = Arrangement.spacedBy(12.dp)
    ) {
        Button(
            onClick = onSave,
            modifier = Modifier.weight(1f).height(40.dp),
            shape = RoundedCornerShape(Dimens.radiusLg),
            colors = ButtonDefaults.buttonColors(
                containerColor = Accent, contentColor = PrimaryDark
            )
        ) {
            Text("保存为图片",
                fontSize = Dimens.fontBody, fontWeight = FontWeight.Medium)
        }
        Button(
            onClick = onShare,
            modifier = Modifier.weight(1f).height(40.dp),
            shape = RoundedCornerShape(Dimens.radiusLg),
            colors = ButtonDefaults.buttonColors(
                containerColor = Primary, contentColor = OnPrimary
            )
        ) {
            Text("分享到",
                fontSize = Dimens.fontBody, fontWeight = FontWeight.Medium)
        }
    }
}

@Composable
private fun SimpleReport(
    arg: BaziResultArg,
    isFemale: Boolean,
    modifier: Modifier = Modifier
) {
    val r = arg.result
    // 由年支决定生肖水印 (第一柱即年柱)
    val yearZhi = r.columns.firstOrNull()?.zhi ?: "子"
    val zodiacRes = zodiacDrawable(yearZhi)

    Box(
        modifier
            .padding(horizontal = Dimens.paddingSm, vertical = Dimens.paddingSm)
            .fillMaxWidth()
            .clip(RoundedCornerShape(10.dp))
            .border(1.dp, GoldLine, RoundedCornerShape(10.dp))
            .background(PaperBg) // 兜底纸色, 防止图片未加载时露白
    ) {
        // 三层结构 (对齐鸿蒙 simpleReport() + Python Bazi.py 原版命书):
        //   ① 底: 宣纸 bz_paper (matchParentSize, 只在 Column 决定 Box 尺寸后按父尺寸铺满, 不撑高 Box)
        //   ② 中: 生肖印 bz_<zodiac> (240dp, 以 Box 中心居中, alpha=1.0)
        //   ③ 顶: 命书正文 (Column, 唯一决定 Box 高度, 文字覆盖生肖印之上)
        // 说明: 生肖 PNG 本身水墨风格 (大量透明像素), 不再二次降透明度,
        //      视觉与 Python img.paste(zodiac, mask=zodiac) 一致.

        // ① 底: 宣纸背景图 (matchParentSize 不参与父测量, Box 高度仍由 Column 决定)
        Image(
            painter = painterResource(R.drawable.bz_paper),
            contentDescription = null,
            contentScale = ContentScale.Crop,
            modifier = Modifier.matchParentSize()
        )

        // ② 中: 生肖印 (以卡片中心居中)
        Image(
            painter = painterResource(zodiacRes),
            contentDescription = null,
            contentScale = ContentScale.Fit,
            alpha = 1.0f,
            modifier = Modifier
                .align(Alignment.Center)
                .size(200.dp)
        )

        // ③ 顶: 命书正文 (紧凑排布, 力求单屏完整显示)
        Column(
            Modifier
                .fillMaxWidth()
                .padding(horizontal = 12.dp, vertical = 10.dp)
        ) {
            // 标题 (单行, 不占额外竖向空间)
            Row(Modifier.fillMaxWidth().padding(bottom = 4.dp)) {
                Spacer(Modifier.weight(1f))
                Text("☯ 八字命书",
                    fontSize = 12.sp, fontWeight = FontWeight.Bold, color = RedColor)
            }

            // 基本信息: 姓名/性别/生肖/年龄 合并为一行
            InfoQuad(
                "姓名", strip(r.userName),
                "性别", if (isFemale) "女" else "男",
                "生肖", strip(r.shengXiao),
                "年龄", strip(r.age)
            )
            // 出生地 / 时间标准 一行
            InfoDuo(
                "出生地", "东经%.2f° 北纬%.2f°".format(arg.longitude, arg.latitude), false,
                "时间标准", if (arg.astEnabled) "真太阳时" else "北京时间(120°E)", false
            )
            // 公历 / 农历各占整行 (内容较长)
            SinglePair("公历", strip(r.solarBirth), false)
            SinglePair("农历", strip(r.lunarBirth), true)

            Spacer(Modifier.height(3.dp))
            Text(jieQiParagraph(arg),
                fontSize = 11.sp, color = InkSoft, lineHeight = 14.sp)
            Spacer(Modifier.height(4.dp))
            ThinDivider(color = GoldLine, alpha = 0.6f)
            Spacer(Modifier.height(4.dp))

            SiZhuBlock(arg.result.columns, isFemale)

            Spacer(Modifier.height(4.dp))
            ThinDivider(color = GoldLine, alpha = 0.6f)
            Spacer(Modifier.height(4.dp))

            Text("${strip(r.qiYun)}　${strip(r.jiaoYun)}",
                fontSize = 11.sp, color = InkSoft, lineHeight = 14.sp)
            Spacer(Modifier.height(3.dp))

            DaYunGrid(arg)
        }
    }
}

@Composable
private fun ThinDivider(color: Color, alpha: Float = 1f) {
    Box(
        Modifier
            .fillMaxWidth()
            .height(1.dp)
            .background(color.copy(alpha = alpha))
    )
}

// 单个『键』值 字段 (紧凑, 单行不换行, 行高收紧)
@Composable
private fun InfoField(k: String, v: String, modifier: Modifier, red: Boolean = false) {
    Row(modifier, verticalAlignment = Alignment.Top) {
        Text("『$k』", fontSize = 12.sp, color = InkSoft, maxLines = 1, lineHeight = 13.sp)
        Text(v, fontSize = 12.sp,
            color = if (red) RedColor else InkColor,
            maxLines = 1, lineHeight = 13.sp,
            modifier = Modifier.padding(start = 1.dp))
    }
}

// 四字段一行 (姓名/性别/生肖/年龄), 对齐命盘图紧凑排版
@Composable
private fun InfoQuad(
    k1: String, v1: String,
    k2: String, v2: String,
    k3: String, v3: String,
    k4: String, v4: String
) {
    Row(Modifier.fillMaxWidth().padding(bottom = 1.dp)) {
        InfoField(k1, v1, Modifier.weight(1.3f))
        InfoField(k2, v2, Modifier.weight(0.9f))
        InfoField(k3, v3, Modifier.weight(0.9f))
        InfoField(k4, v4, Modifier.weight(1.1f))
    }
}

// 两字段一行 (出生地/时间标准)
@Composable
private fun InfoDuo(
    k1: String, v1: String, red1: Boolean,
    k2: String, v2: String, red2: Boolean
) {
    Row(Modifier.fillMaxWidth().padding(bottom = 1.dp)) {
        InfoField(k1, v1, Modifier.weight(1.1f), red1)
        InfoField(k2, v2, Modifier.weight(1f), red2)
    }
}

@Composable
private fun SinglePair(k: String, v: String, red: Boolean) {
    Row(Modifier.fillMaxWidth().padding(bottom = 1.dp), verticalAlignment = Alignment.Top) {
        Text("『$k』", fontSize = 12.sp, color = InkSoft, maxLines = 1, lineHeight = 13.sp)
        Text(v, fontSize = 12.sp,
            color = if (red) RedColor else InkColor,
            maxLines = 1, lineHeight = 13.sp,
            modifier = Modifier.padding(start = 1.dp))
    }
}

@Composable
private fun SiZhuBlock(columns: List<BaziColumn>, isFemale: Boolean) {
    Row(Modifier.fillMaxWidth(), verticalAlignment = Alignment.Top) {
        // 乾造/坤造 竖排标签, 与天干行对齐 (顶部留白跳过主星行)
        Text(if (isFemale) "坤造" else "乾造",
            fontSize = 12.sp, color = InkColor,
            fontFamily = WenYueFont,
            modifier = Modifier.width(24.dp).padding(top = 16.dp))
        columns.forEachIndexed { i, c ->
            Column(
                Modifier.weight(1f),
                horizontalAlignment = Alignment.CenterHorizontally
            ) {
                Text(if (i == 2) "日元" else c.ganShiShen,
                    fontSize = 11.sp, fontWeight = FontWeight.Medium,
                    lineHeight = 12.sp,
                    color = if (i == 2) RedColor else InkSoft)
                Text(c.gan, fontSize = 18.sp,
                    fontFamily = WenYueFont, lineHeight = 19.sp,
                    color = if (i == 2) RedColor else InkColor,
                    modifier = Modifier.padding(top = 1.dp))
                Text(c.zhi, fontSize = 18.sp,
                    fontFamily = WenYueFont, lineHeight = 19.sp,
                    color = InkColor,
                    modifier = Modifier.padding(top = 1.dp))
                Column(
                    Modifier.padding(top = 3.dp),
                    horizontalAlignment = Alignment.CenterHorizontally
                ) {
                    c.cangGan.forEach { cg ->
                        Row {
                            Text(cg.gan, fontSize = 10.sp, lineHeight = 11.sp,
                                fontFamily = WenYueFont, color = InkColor)
                            Text(cg.shiShen, fontSize = 10.sp, lineHeight = 11.sp,
                                color = InkSoft,
                                modifier = Modifier.padding(start = 2.dp))
                        }
                    }
                }
            }
        }
    }
}

private fun jieQiParagraph(arg: BaziResultArg): String {
    var s = ""
    val dq = strip(arg.result.dingQiType)
    if (dq.isNotEmpty()) s += "依据${dq}。"
    val jq = arg.result.jieQi.trim().replace(Regex("\n+"), "；")
    if (jq.isNotEmpty()) s += "${jq}。"
    if (arg.result.siLing.isNotEmpty()) s += "命主月令司令：${arg.result.siLing}"
    return s
}

@Composable
private fun DaYunGrid(arg: BaziResultArg) {
    val cols = arg.result.daYunColumns.take(8)
    Row(Modifier.fillMaxWidth(), verticalAlignment = Alignment.Top) {
        // 大運 竖排标签, 与干支行对齐 (顶部留白跳过起运岁数行)
        Text("大運",
            fontSize = 12.sp, color = InkColor, fontFamily = WenYueFont,
            modifier = Modifier.width(22.dp).padding(top = 18.dp))
        cols.forEach { c ->
            DaYunSimpleCell(arg, c, Modifier.weight(1f))
        }
    }
}

@Composable
private fun DaYunSimpleCell(arg: BaziResultArg, c: BaziColumn, modifier: Modifier) {
    val curDaYun = arg.result.currentDaYun?.startYear ?: -1
    val curLN = arg.result.currentLiuNian?.startYear ?: -1
    val isCur = c.startYear == curDaYun
    Column(
        modifier.padding(horizontal = 1.dp),
        horizontalAlignment = Alignment.CenterHorizontally
    ) {
        Text("${c.startYear - arg.birthYear + 1}",
            fontSize = 10.sp, fontWeight = FontWeight.Medium,
            color = if (isCur) RedColor else InkSoft)
        Column(
            Modifier
                .padding(horizontal = 2.dp, vertical = 1.dp)
                .clip(RoundedCornerShape(3.dp))
                .background(if (isCur) RedBg else Color.Transparent),
            horizontalAlignment = Alignment.CenterHorizontally
        ) {
            Text(c.gan, fontSize = 16.sp,
                fontFamily = WenYueFont,
                color = if (isCur) RedColor else InkColor)
            Text(c.zhi, fontSize = 16.sp,
                fontFamily = WenYueFont,
                color = if (isCur) RedColor else InkColor)
        }
        Text("${c.startYear}", fontSize = 8.sp, color = RedColor,
            modifier = Modifier.padding(top = 1.dp))
        ThinDivider(color = GoldLine, alpha = 0.5f)
        Column(horizontalAlignment = Alignment.CenterHorizontally) {
            (0 until 10).forEach { j ->
                val y = c.startYear + j
                val isCurLN = y == curLN
                Text(
                    BaziCalc.GAN[BaziCalc.lnGan(y)]
                        + BaziCalc.ZHI[BaziCalc.lnZhi(y)],
                    fontSize = 12.sp, maxLines = 1, lineHeight = 13.sp,
                    fontFamily = WenYueFont,
                    color = if (isCurLN) RedColor else InkColor,
                    modifier = Modifier
                        .padding(horizontal = 1.dp)
                        .clip(RoundedCornerShape(3.dp))
                        .background(if (isCurLN) RedBg else Color.Transparent)
                )
            }
        }
        Text("${c.startYear + 9}", fontSize = 8.sp, color = RedColor,
            modifier = Modifier.padding(top = 1.dp))
    }
}

// ════════════════════ 专业版 ════════════════════

@Composable
private fun InfoBar(arg: BaziResultArg, isFemale: Boolean) {
    Column(
        Modifier
            .padding(horizontal = Dimens.paddingMd, vertical = Dimens.paddingMd)
            .fillMaxWidth()
            .clip(RoundedCornerShape(Dimens.radiusMd))
            .background(Surface)
            .padding(Dimens.paddingMd)
    ) {
        Text("${if (isFemale) "坤" else "乾"}  ${strip(arg.result.userName)}",
            fontSize = Dimens.fontSubtitle, fontWeight = FontWeight.Bold,
            color = Primary)
        Text("${strip(arg.result.shengXiao)}  ${strip(arg.result.age)}",
            fontSize = Dimens.fontCaption, color = TextSecondary,
            modifier = Modifier.padding(top = 2.dp))
        Text("公历 ${strip(arg.result.solarBirth)}",
            fontSize = Dimens.fontCaption, color = TextSecondary,
            modifier = Modifier.padding(top = 2.dp))
        Text("农历 ${strip(arg.result.lunarBirth)}",
            fontSize = Dimens.fontCaption, color = TextSecondary,
            modifier = Modifier.padding(top = 2.dp))
    }
}

@Composable
private fun PaiPanTable(arg: BaziResultArg) {
    val columns = buildAllColumns(arg)
    Column(
        Modifier
            .padding(horizontal = Dimens.paddingMd, vertical = Dimens.paddingMd)
            .fillMaxWidth()
            .clip(RoundedCornerShape(Dimens.radiusMd))
            .background(Surface)
            .padding(Dimens.paddingSm)
            .horizontalScroll(rememberScrollState())
    ) {
        // header
        Row(
            Modifier.height(28.dp).fillMaxWidth(),
            verticalAlignment = Alignment.CenterVertically
        ) {
            Box(Modifier.width(LabelW.dp))
            ColTitles.forEachIndexed { i, t ->
                Text(t, fontSize = Dimens.fontCaption,
                    fontWeight = FontWeight.Bold,
                    color = if (i >= 4) Accent else Primary,
                    textAlign = TextAlign.Center,
                    modifier = Modifier.width(ColW.dp))
            }
        }
        ThinDivider(color = DividerColor)

        TableRow("主星", columns.mapIndexed { i, c ->
            ColorBox(c.ganShiShen,
                if (i == 2) Primary else TextSecondary,
                bold = i == 2)
        })
        TableRow("天干", columns.map { ColorBox(it.gan, charColor(it.gan), 26.sp, true) })
        TableRow("地支", columns.map { ColorBox(it.zhi, charColor(it.zhi), 26.sp, true) })

        Row(
            Modifier.padding(vertical = 4.dp).fillMaxWidth(),
            verticalAlignment = Alignment.Top
        ) {
            LabelCell("藏干")
            columns.forEach { c ->
                Column(
                    Modifier.width(ColW.dp),
                    horizontalAlignment = Alignment.CenterHorizontally
                ) {
                    c.cangGan.forEach { cg ->
                        Row {
                            Text(cg.gan,
                                fontSize = Dimens.fontSmall, color = charColor(cg.gan))
                            Text(cg.shiShen,
                                fontSize = Dimens.fontSmall, color = TextSecondary,
                                modifier = Modifier.padding(start = 2.dp))
                        }
                    }
                }
            }
        }

        TableRow("纳音", columns.map { ColorBox(it.nayin, OnSurface) })
        TableRow("星运", columns.map { ColorBox(it.xingYun, OnSurface) })
        TableRow("自坐", columns.map { ColorBox(it.ziZuo, OnSurface) })
        TableRow("空亡", columns.map { ColorBox(it.kongWang, OnSurface) })

        ThinDivider(color = DividerColor)
        Row(
            Modifier.padding(vertical = 4.dp).fillMaxWidth(),
            verticalAlignment = Alignment.Top
        ) {
            LabelCell("神煞")
            columns.forEach { c ->
                Column(
                    Modifier.width(ColW.dp),
                    horizontalAlignment = Alignment.CenterHorizontally
                ) {
                    c.shenSha.forEach { s ->
                        Text(s, fontSize = 9.sp, color = TextSecondary,
                            modifier = Modifier.padding(bottom = 1.dp))
                    }
                }
            }
        }
    }
}

private data class ColorBox(
    val text: String,
    val color: Color,
    val fontSize: TextUnit = TextUnit.Unspecified,
    val bold: Boolean = false
)

@Composable
private fun TableRow(label: String, cells: List<ColorBox>) {
    val firstFs = cells.firstOrNull()?.fontSize ?: TextUnit.Unspecified
    val rowHeight = if (firstFs != TextUnit.Unspecified && firstFs.value > 20f) 44.dp else 26.dp
    Row(
        Modifier.height(rowHeight).fillMaxWidth(),
        verticalAlignment = Alignment.CenterVertically
    ) {
        LabelCell(label)
        cells.forEach { box ->
            Text(box.text,
                fontSize = if (box.fontSize != TextUnit.Unspecified) box.fontSize else Dimens.fontSmall,
                fontWeight = if (box.bold) FontWeight.Bold else FontWeight.Normal,
                color = box.color,
                textAlign = TextAlign.Center,
                modifier = Modifier.width(ColW.dp))
        }
    }
}

@Composable
private fun LabelCell(label: String) {
    Text(label, fontSize = Dimens.fontSmall, color = TextSecondary,
        textAlign = TextAlign.Center,
        modifier = Modifier.width(LabelW.dp))
}

@Composable
private fun SummaryBar(arg: BaziResultArg) {
    Column(
        Modifier
            .padding(horizontal = Dimens.paddingMd, vertical = Dimens.paddingMd)
            .fillMaxWidth()
            .clip(RoundedCornerShape(Dimens.radiusMd))
            .background(Surface)
            .padding(Dimens.paddingMd)
    ) {
        Row(
            Modifier.fillMaxWidth().padding(bottom = Dimens.paddingSm),
            verticalAlignment = Alignment.CenterVertically
        ) {
            Text("司令", fontSize = Dimens.fontCaption, color = TextSecondary)
            Text(arg.result.siLing,
                fontSize = Dimens.fontBody, fontWeight = FontWeight.Bold,
                color = Primary, modifier = Modifier.padding(start = 4.dp))
            Spacer(Modifier.weight(1f))
            Text(wuXingSummary(arg),
                fontSize = Dimens.fontCaption, color = OnSurface)
        }
        Row(Modifier.fillMaxWidth()) {
            statusChips(arg).forEachIndexed { i, label ->
                Box(
                    Modifier
                        .weight(1f)
                        .padding(end = if (i < 4) 4.dp else 0.dp)
                        .height(28.dp)
                        .clip(RoundedCornerShape(Dimens.radiusSm))
                        .background(ElemColors.getOrElse(i) { Primary }),
                    contentAlignment = Alignment.Center
                ) {
                    Text(label, fontSize = Dimens.fontCaption, color = OnPrimary)
                }
            }
        }
    }
}

@Composable
private fun LiuNianCard(arg: BaziResultArg) {
    val curDaYun = arg.result.currentDaYun?.startYear ?: -1
    val curLN = arg.result.currentLiuNian?.startYear ?: -1
    Column(
        Modifier
            .padding(horizontal = Dimens.paddingMd, vertical = Dimens.paddingMd)
            .fillMaxWidth()
            .clip(RoundedCornerShape(Dimens.radiusLg))
            .background(Surface)
            .padding(Dimens.paddingLg)
    ) {
        Text("大运", fontSize = Dimens.fontSubtitle, fontWeight = FontWeight.Bold,
            color = Primary,
            modifier = Modifier.padding(bottom = Dimens.paddingSm))
        Row(Modifier.horizontalScroll(rememberScrollState())) {
            arg.result.daYunColumns.forEach { c ->
                val isCur = c.startYear == curDaYun
                Column(
                    Modifier
                        .width(48.dp)
                        .padding(end = 4.dp)
                        .clip(RoundedCornerShape(Dimens.radiusSm))
                        .background(if (isCur) RedBg else Background)
                        .padding(vertical = 4.dp),
                    horizontalAlignment = Alignment.CenterHorizontally
                ) {
                    Text(c.ganShiShen, fontSize = 9.sp, color = TextSecondary)
                    Text(c.gan, fontSize = 18.sp, fontWeight = FontWeight.Bold,
                        color = charColor(c.gan))
                    Text(c.zhi, fontSize = 18.sp, fontWeight = FontWeight.Bold,
                        color = charColor(c.zhi))
                    Text("${c.startYear}", fontSize = 9.sp, color = TextSecondary,
                        modifier = Modifier.padding(top = 2.dp))
                }
            }
        }

        Text("流年 / 小运（当前大运）",
            fontSize = Dimens.fontSubtitle, fontWeight = FontWeight.Bold,
            color = Primary,
            modifier = Modifier
                .padding(top = Dimens.paddingMd, bottom = Dimens.paddingSm))
        Row(Modifier.horizontalScroll(rememberScrollState())) {
            arg.result.liuNian.forEach { it ->
                LiuNianCell(it, curLN)
            }
        }
    }
}

@Composable
private fun LiuNianCell(item: LiuNianItem, curLN: Int) {
    val isCur = item.year == curLN
    Column(
        Modifier
            .width(54.dp)
            .padding(end = 4.dp)
            .clip(RoundedCornerShape(Dimens.radiusSm))
            .background(if (isCur) RedBg else Background)
            .padding(vertical = 4.dp),
        horizontalAlignment = Alignment.CenterHorizontally
    ) {
        Text("${item.year}", fontSize = 9.sp, color = TextSecondary)
        Text("${item.age}岁", fontSize = 9.sp, color = LunarText)
        Row(modifier = Modifier.padding(top = 2.dp)) {
            Text(item.ganZhi,
                fontSize = 13.sp, fontWeight = FontWeight.Bold, color = OnSurface)
            Text(item.ganShiShen,
                fontSize = 9.sp, color = Accent,
                modifier = Modifier.padding(start = 1.dp))
        }
        Row {
            Text(item.xiaoYun, fontSize = 11.sp, color = TextSecondary)
            Text(item.xiaoYunShiShen,
                fontSize = 9.sp, color = LunarText,
                modifier = Modifier.padding(start = 1.dp))
        }
    }
}

// ─── 工具 ────────────────────────────────────────────────

private fun buildAllColumns(arg: BaziResultArg): List<BaziColumn> {
    val list = mutableListOf<BaziColumn>()
    list += arg.result.columns
    arg.result.currentDaYun?.let { list += it }
    arg.result.currentLiuNian?.let { list += it }
    return list
}

private fun strip(s: String?): String {
    if (s.isNullOrEmpty()) return ""
    val i = s.indexOf('』')
    return if (i >= 0) s.substring(i + 1).trim() else s.trim()
}

private fun charColor(ch: String): Color {
    val ganRange = "甲乙丙丁戊己庚辛壬癸"
    val gi = ganRange.indexOf(ch)
    if (gi >= 0) return ElemColors[gi / 2]
    val zhiRange = "子丑寅卯辰巳午未申酉戌亥"
    val zi = zhiRange.indexOf(ch)
    if (zi >= 0) {
        val zw = intArrayOf(4, 2, 0, 0, 2, 1, 1, 2, 3, 3, 2, 4)
        return ElemColors[zw[zi]]
    }
    return OnSurface
}

private fun wuXingSummary(arg: BaziResultArg): String {
    val names = listOf("木", "火", "土", "金", "水")
    if (arg.result.wuXingCount.isEmpty()) return ""
    var s = "含藏气: "
    for (i in 0..4) {
        if (i < arg.result.wuXingCount.size) {
            s += "${arg.result.wuXingCount[i]}${names[i]} "
        }
    }
    return s.trim()
}

private fun statusChips(arg: BaziResultArg): List<String> {
    val names = listOf("木", "火", "土", "金", "水")
    return (0..4).map { i ->
        val st = arg.result.wuXingStatus.getOrNull(i) ?: ""
        "${names[i]}$st"
    }
}

// ── 录制 Composable 到 Picture (Compose 1.6 兼容, 不依赖 1.7 rememberGraphicsLayer) ──
// 挂在待截图区域最外层; onDrawWithContent 内先录制到 picture 再回放到真实画布,
// 因此背景/子内容都会被捕获, 且滚动容器中仍录制完整高度。
private fun Modifier.recordToPicture(picture: Picture): Modifier = this.drawWithCache {
    val w = size.width.toInt().coerceAtLeast(1)
    val h = size.height.toInt().coerceAtLeast(1)
    onDrawWithContent {
        val pictureCanvas = Canvas(picture.beginRecording(w, h))
        draw(this, this.layoutDirection, pictureCanvas, this.size) {
            this@onDrawWithContent.drawContent()
        }
        picture.endRecording()
        drawIntoCanvas { it.nativeCanvas.drawPicture(picture) }
    }
}

// ── 保存到相册 ────────────────────────────────────────────
private fun saveBaziImage(ctx: Context, scope: CoroutineScope, picture: Picture) {
    scope.launch {
        val bmp = pictureToBitmap(picture)
        if (bmp == null) {
            Toast.makeText(ctx, "截图失败,请稍后再试", Toast.LENGTH_SHORT).show()
            return@launch
        }
        val ok = withContext(Dispatchers.IO) {
            try {
                saveBitmapToGallery(ctx, bmp)
            } finally {
                if (!bmp.isRecycled) bmp.recycle()
            }
        }
        Toast.makeText(ctx, if (ok) "已保存到相册" else "保存失败",
            Toast.LENGTH_SHORT).show()
    }
}

// ── 分享到其它 app (系统分享面板) ─────────────────────────
private fun shareBaziImage(ctx: Context, scope: CoroutineScope, picture: Picture) {
    scope.launch {
        val bmp = pictureToBitmap(picture)
        if (bmp == null) {
            Toast.makeText(ctx, "截图失败,请稍后再试", Toast.LENGTH_SHORT).show()
            return@launch
        }
        // saveBitmap 落盘后会 recycle bmp
        val file = withContext(Dispatchers.IO) {
            EclipseShareUtil.saveBitmap(ctx, bmp, "bazi_${System.currentTimeMillis()}.png")
        }
        EclipseShareUtil.shareImage(ctx, file, "分享八字命盘")
    }
}

// ── 保存图片 (简洁版命书) ─────────────────────────────────
// Picture → Bitmap. Android 28+ 用 Bitmap.createBitmap(Picture) 一步到位;
// 老版本回落到 Canvas.drawPicture() (需要 software canvas)。
private fun pictureToBitmap(picture: Picture): Bitmap? {
    val w = picture.width
    val h = picture.height
    if (w <= 0 || h <= 0) return null
    return try {
        if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.P) {
            Bitmap.createBitmap(picture)
        } else {
            val bmp = Bitmap.createBitmap(w, h, Bitmap.Config.ARGB_8888)
            android.graphics.Canvas(bmp).drawPicture(picture)
            bmp
        }
    } catch (_: Throwable) {
        null
    }
}

// 写入 PNG 到系统相册。
// Android 10+ 走 MediaStore.IS_PENDING 流程; 旧版本写到外部 Pictures 目录。
private fun saveBitmapToGallery(ctx: Context, bitmap: Bitmap): Boolean {
    val name = "bazi_${System.currentTimeMillis()}.png"
    return try {
        if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.Q) {
            val values = ContentValues().apply {
                put(MediaStore.Images.Media.DISPLAY_NAME, name)
                put(MediaStore.Images.Media.MIME_TYPE, "image/png")
                put(MediaStore.Images.Media.RELATIVE_PATH,
                    Environment.DIRECTORY_PICTURES + "/SxwnlCalendar")
                put(MediaStore.Images.Media.IS_PENDING, 1)
            }
            val uri = ctx.contentResolver.insert(
                MediaStore.Images.Media.EXTERNAL_CONTENT_URI, values
            ) ?: return false
            ctx.contentResolver.openOutputStream(uri)?.use {
                bitmap.compress(Bitmap.CompressFormat.PNG, 100, it)
            } ?: return false
            values.clear()
            values.put(MediaStore.Images.Media.IS_PENDING, 0)
            ctx.contentResolver.update(uri, values, null, null)
            true
        } else {
            val dir = File(
                Environment.getExternalStoragePublicDirectory(Environment.DIRECTORY_PICTURES),
                "SxwnlCalendar"
            ).apply { if (!exists()) mkdirs() }
            FileOutputStream(File(dir, name)).use {
                bitmap.compress(Bitmap.CompressFormat.PNG, 100, it)
            }
            true
        }
    } catch (_: Throwable) {
        false
    }
}
