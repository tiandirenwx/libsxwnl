package com.sxwnl.calendar.bridge

import com.sxwnl.calendar.data.AlmanacTopic
import com.sxwnl.calendar.data.BaziColumn
import com.sxwnl.calendar.data.DayAlmanac
import com.sxwnl.calendar.data.DayInfo
import com.sxwnl.calendar.data.DayRTS
import com.sxwnl.calendar.data.GeoCity
import com.sxwnl.calendar.data.GeoProvince
import com.sxwnl.calendar.data.JieQiItem
import com.sxwnl.calendar.data.LiuNianItem
import com.sxwnl.calendar.data.LocalSolarEclipse
import com.sxwnl.calendar.data.LunarEclipseDetail
import com.sxwnl.calendar.data.LunarEclipseItem
import com.sxwnl.calendar.data.LunarMonth
import com.sxwnl.calendar.data.ReverseItem
import com.sxwnl.calendar.data.SolarEclipseItem
import com.sxwnl.calendar.data.SolarEclipsePath
import com.sxwnl.calendar.data.TimezoneGroup
import com.sxwnl.calendar.data.YearCalMonth

/**
 * 与 capi/sxwnl_capi.h 对齐的 JNI 桥接面.
 * 全部 native 方法都是无状态的, 八字使用 handle (Long) 持有 C++ 对象。
 */
object SxwnlNative {

    init {
        System.loadLibrary("sxwnl_jni")
    }

    // ═══ Calendar ════════════════════════════════════════════

    @JvmStatic external fun getDayInfo(year: Int, month: Int, day: Int): DayInfo?
    @JvmStatic external fun getMonthData(year: Int, month: Int): Array<DayInfo>?

    @JvmStatic external fun lunarToSolar(
        year: Int, month: Int, day: Int, isLeap: Boolean
    ): IntArray?

    @JvmStatic external fun solarToLunar(year: Int, month: Int, day: Int): DayInfo?

    @JvmStatic external fun getJieQiList(year: Int): Array<JieQiItem>?
    @JvmStatic external fun getYearLeapMonth(year: Int): Int
    @JvmStatic external fun getLunarMonths(year: Int): Array<LunarMonth>?
    @JvmStatic external fun getLunarMonthDays(
        year: Int, month: Int, isLeap: Boolean, isSpec: Boolean
    ): Int
    @JvmStatic external fun getSolarMonthValidDays(year: Int, month: Int): IntArray?

    // ═══ Year Calendar ═══════════════════════════════════════

    @JvmStatic external fun getYearCalendar(year: Int): Array<YearCalMonth>?

    // ═══ RTS (日月升降/中天) ══════════════════════════════════

    @JvmStatic external fun calcDayRTS(
        year: Int, month: Int, day: Int,
        longitude: Double, latitude: Double, tzHours: Double
    ): DayRTS?

    // ═══ Bazi - 基础 (handle 模式) ════════════════════════════

    @JvmStatic external fun baziCreate(
        name: String,
        gender: Boolean,
        isAst: Boolean,
        isLunar: Boolean,
        isLeap: Boolean,
        isSpec: Boolean,
        year: Int,
        month: Int,
        day: Int,
        hour: Int,
        minute: Int,
        second: Double,
        longitude: Double,
        latitude: Double,
        lifa: Int
    ): Long

    @JvmStatic external fun baziFree(handle: Long)

    @JvmStatic external fun baziGetUserName(handle: Long): String
    @JvmStatic external fun baziGetGender(handle: Long): String
    @JvmStatic external fun baziGetSolarBirth(handle: Long): String
    @JvmStatic external fun baziGetLunarBirth(handle: Long): String
    @JvmStatic external fun baziGetDateOfBirth(handle: Long): String
    @JvmStatic external fun baziGetShengXiao(handle: Long): String
    @JvmStatic external fun baziGetAge(handle: Long): String
    @JvmStatic external fun baziGetLifa(handle: Long): String
    @JvmStatic external fun baziGetDingQiType(handle: Long): String
    @JvmStatic external fun baziGetAst(handle: Long): String
    @JvmStatic external fun baziGetJieQi(handle: Long): String
    @JvmStatic external fun baziGetQiYun(handle: Long): String
    @JvmStatic external fun baziGetJiaoYun(handle: Long): String
    @JvmStatic external fun baziGetSiLing(handle: Long): String

    // ═══ Bazi - 完整盘面 ═════════════════════════════════════

    @JvmStatic external fun baziGetColumns(handle: Long): Array<BaziColumn>?
    @JvmStatic external fun baziGetCurrentDaYun(handle: Long): BaziColumn?
    @JvmStatic external fun baziGetCurrentLiuNian(handle: Long): BaziColumn?
    @JvmStatic external fun baziGetDaYunColumns(handle: Long): Array<BaziColumn>?
    @JvmStatic external fun baziGetWuXingCount(handle: Long, includeCangGan: Boolean): IntArray?
    @JvmStatic external fun baziGetWuXingStatus(handle: Long): Array<String>?
    @JvmStatic external fun baziGetLiuNian(handle: Long, startYear: Int, count: Int): Array<LiuNianItem>?

    // ═══ Bazi - 工具 ═════════════════════════════════════════

    @JvmStatic external fun baziReverse(
        yg: Int, yz: Int, mg: Int, mz: Int,
        dg: Int, dz: Int, hg: Int, hz: Int,
        startYear: Int, endYear: Int
    ): Array<ReverseItem>?

    // ═══ Eclipse Map (结构化, 用于可视化) ══════════════════════

    @JvmStatic external fun searchSolarEclipses(
        year: Int, month: Int, count: Int
    ): Array<SolarEclipseItem>?

    @JvmStatic external fun getSolarEclipsePath(
        year: Int, month: Int, day: Int
    ): SolarEclipsePath?

    @JvmStatic external fun getLocalSolarEclipse(
        year: Int, month: Int, day: Int,
        longitude: Double, latitude: Double, frameInterval: Int
    ): LocalSolarEclipse?

    @JvmStatic external fun searchLunarEclipses(
        year: Int, month: Int, count: Int
    ): Array<LunarEclipseItem>?

    @JvmStatic external fun getLunarEclipseDetail(
        year: Int, month: Int, day: Int, frameInterval: Int
    ): LunarEclipseDetail?

    // ═══ World Map (\u6d77\u5cb8\u7ebf\u8f6e\u5ed3) ═══════════════════════
    //
    // \u8fd4\u56de double[], \u7ecf\u7eac\u5ea6(\u5f27\u5ea6) \u4ea4\u66ff: [lon0, lat0, lon1, lat1, ...]
    // \u6bb5\u95f4\u5206\u9694\u70b9\u7528 1e7 \u6807\u8bb0 (Move-To)\u3002
    @JvmStatic external fun getWorldMapDitu0(): DoubleArray?

    /** idx=1 → ditu1 大图海岸; idx=2 → ditu2 国界 */
    @JvmStatic external fun getWorldMapData(idx: Int): DoubleArray?

    // ═══ 老黄历 (Almanac) ═════════════════════════════════════

    /** 取公历某日老黄历 (二十八宿/黄道黑道/冲煞/方位/彭祖/神煞/宜忌/吉时/用事). */
    @JvmStatic external fun getAlmanac(year: Int, month: Int, day: Int): DayAlmanac?

    /** 老黄历静态知识 (董公总论/口诀/方位 等, 全局只取一次). */
    @JvmStatic external fun getAlmanacTopics(): Array<AlmanacTopic>?

    // ═══ 地理目录 (GeoPostion + JWv/SQv) ══════════════════════
    //
    //   数据完全来自 libsxwnl C++ 内部 (src/geo.cpp), 上层无需维护城市表.

    @JvmStatic external fun geoListProvinces(): Array<GeoProvince>?
    @JvmStatic external fun geoListCities(province: String): Array<GeoCity>?
    @JvmStatic external fun geoSearch(keyword: String, limit: Int): Array<GeoCity>?
    @JvmStatic external fun geoListTimezones(): Array<TimezoneGroup>?
    @JvmStatic external fun geoDefault(): GeoCity?

    // ═══ 纪年字符串工具 ═══════════════════════════════════════

    @JvmStatic external fun yearStrToAstro(s: String): Int
    @JvmStatic external fun astroYearToStr(year: Int, fullStyle: Boolean): String
}
