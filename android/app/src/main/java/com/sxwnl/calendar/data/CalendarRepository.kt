package com.sxwnl.calendar.data

import com.sxwnl.calendar.bridge.SxwnlNative
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.withContext

/**
 * 数据层 — 对 SxwnlNative 的协程封装, 与鸿蒙 NativeBridge.ets 同义。
 *
 *  所有方法均在 IO 调度器执行, UI 线程安全。
 */
object CalendarRepository {

    // ═══ Calendar ════════════════════════════════════════════

    suspend fun getDayInfo(year: Int, month: Int, day: Int): DayInfo? =
        withContext(Dispatchers.IO) { SxwnlNative.getDayInfo(year, month, day) }

    suspend fun getMonthData(year: Int, month: Int): List<DayInfo> =
        withContext(Dispatchers.IO) {
            SxwnlNative.getMonthData(year, month)?.toList() ?: emptyList()
        }

    suspend fun lunarToSolar(
        year: Int, month: Int, day: Int, isLeap: Boolean
    ): SolarDate? = withContext(Dispatchers.IO) {
        SxwnlNative.lunarToSolar(year, month, day, isLeap)?.let {
            SolarDate(it[0], it[1], it[2])
        }
    }

    suspend fun solarToLunar(year: Int, month: Int, day: Int): DayInfo? =
        withContext(Dispatchers.IO) { SxwnlNative.solarToLunar(year, month, day) }

    suspend fun getJieQiList(year: Int): List<JieQiItem> =
        withContext(Dispatchers.IO) {
            SxwnlNative.getJieQiList(year)?.toList() ?: emptyList()
        }

    suspend fun getYearLeapMonth(year: Int): Int =
        withContext(Dispatchers.IO) { SxwnlNative.getYearLeapMonth(year) }

    suspend fun getLunarMonths(year: Int): List<LunarMonth> =
        withContext(Dispatchers.IO) {
            SxwnlNative.getLunarMonths(year)?.toList() ?: emptyList()
        }

    suspend fun getLunarMonthDays(
        year: Int, month: Int, isLeap: Boolean, isSpec: Boolean
    ): Int = withContext(Dispatchers.IO) {
        SxwnlNative.getLunarMonthDays(year, month, isLeap, isSpec)
    }

    suspend fun getSolarMonthValidDays(year: Int, month: Int): List<Int> =
        withContext(Dispatchers.IO) {
            SxwnlNative.getSolarMonthValidDays(year, month)?.toList() ?: emptyList()
        }

    // ═══ Year Calendar ═══════════════════════════════════════

    suspend fun getYearCalendar(year: Int): List<YearCalMonth> =
        withContext(Dispatchers.IO) {
            SxwnlNative.getYearCalendar(year)?.toList() ?: emptyList()
        }

    // ═══ RTS ═════════════════════════════════════════════════

    suspend fun calcDayRTS(
        year: Int, month: Int, day: Int,
        longitude: Double, latitude: Double, tzHours: Double
    ): DayRTS? = withContext(Dispatchers.IO) {
        SxwnlNative.calcDayRTS(year, month, day, longitude, latitude, tzHours)
    }

    // ═══ Bazi (完整封装) ═════════════════════════════════════

    /**
     * 与鸿蒙端 NativeBridge.calcBazi 对齐. 内部 handle 在函数返回前释放,
     * BaziResult 中的所有 getter 已直接拷贝为 Kotlin 字段, 调用方无需关心生命周期。
     */
    suspend fun calcBazi(params: BaziParams): BaziResult? = withContext(Dispatchers.IO) {
        val handle = SxwnlNative.baziCreate(
            params.name,
            params.gender,
            params.astEnabled,
            params.inputMode == 1,         // isLunar
            params.isLeap,
            params.isSpec,
            params.year, params.month, params.day,
            params.hour, params.minute, 0.0,
            params.longitude, params.latitude,
            params.lifa
        )
        if (handle == 0L) return@withContext null

        try {
            buildResult(handle)
        } finally {
            SxwnlNative.baziFree(handle)
        }
    }

    private fun buildResult(handle: Long): BaziResult {
        val solarBirth = SxwnlNative.baziGetSolarBirth(handle)
        // 起始年: 从 solarBirth 取数, 失败用 1900
        val startYear = parseYearFrom(solarBirth) ?: 1900
        return BaziResult().apply {
            userName    = SxwnlNative.baziGetUserName(handle)
            gender      = SxwnlNative.baziGetGender(handle)
            this.solarBirth = solarBirth
            lunarBirth  = SxwnlNative.baziGetLunarBirth(handle)
            dateOfBirth = SxwnlNative.baziGetDateOfBirth(handle)
            shengXiao   = SxwnlNative.baziGetShengXiao(handle)
            age         = SxwnlNative.baziGetAge(handle)
            lifa        = SxwnlNative.baziGetLifa(handle)
            dingQiType  = SxwnlNative.baziGetDingQiType(handle)
            ast         = SxwnlNative.baziGetAst(handle)
            jieQi       = SxwnlNative.baziGetJieQi(handle)
            qiYun       = SxwnlNative.baziGetQiYun(handle)
            jiaoYun     = SxwnlNative.baziGetJiaoYun(handle)
            siLing      = SxwnlNative.baziGetSiLing(handle)
            columns     = SxwnlNative.baziGetColumns(handle)?.toList() ?: emptyList()
            currentDaYun   = SxwnlNative.baziGetCurrentDaYun(handle)
            currentLiuNian = SxwnlNative.baziGetCurrentLiuNian(handle)
            daYunColumns   = SxwnlNative.baziGetDaYunColumns(handle)?.toList() ?: emptyList()
            wuXingCount    = SxwnlNative.baziGetWuXingCount(handle, true) ?: IntArray(5)
            wuXingStatus   = SxwnlNative.baziGetWuXingStatus(handle)
                ?: arrayOf("", "", "", "", "")
            liuNian        = SxwnlNative.baziGetLiuNian(handle, startYear, 120)?.toList() ?: emptyList()
        }
    }

    /** "公历 2024年01月01日 ..." → 2024 */
    private fun parseYearFrom(birth: String): Int? {
        val m = Regex("""(-?\d+)年""").find(birth) ?: return null
        return m.groupValues[1].toIntOrNull()
    }

    // ═══ Bazi - 工具 ═════════════════════════════════════════

    suspend fun baziReverse(
        yg: Int, yz: Int, mg: Int, mz: Int,
        dg: Int, dz: Int, hg: Int, hz: Int,
        startYear: Int, endYear: Int
    ): List<ReverseItem> = withContext(Dispatchers.IO) {
        SxwnlNative.baziReverse(yg, yz, mg, mz, dg, dz, hg, hz, startYear, endYear)
            ?.toList() ?: emptyList()
    }

    // ═══ Eclipse Map (结构化) ════════════════════════════════

    suspend fun searchSolarEclipses(
        year: Int, month: Int, count: Int
    ): List<SolarEclipseItem> = withContext(Dispatchers.IO) {
        SxwnlNative.searchSolarEclipses(year, month, count)?.toList() ?: emptyList()
    }

    suspend fun getSolarEclipsePath(
        year: Int, month: Int, day: Int
    ): SolarEclipsePath? = withContext(Dispatchers.IO) {
        SxwnlNative.getSolarEclipsePath(year, month, day)
    }

    suspend fun getLocalSolarEclipse(
        year: Int, month: Int, day: Int,
        longitude: Double, latitude: Double, frameInterval: Int = 60
    ): LocalSolarEclipse? = withContext(Dispatchers.IO) {
        SxwnlNative.getLocalSolarEclipse(year, month, day, longitude, latitude, frameInterval)
    }

    suspend fun searchLunarEclipses(
        year: Int, month: Int, count: Int
    ): List<LunarEclipseItem> = withContext(Dispatchers.IO) {
        SxwnlNative.searchLunarEclipses(year, month, count)?.toList() ?: emptyList()
    }

    suspend fun getLunarEclipseDetail(
        year: Int, month: Int, day: Int, frameInterval: Int = 60
    ): LunarEclipseDetail? = withContext(Dispatchers.IO) {
        SxwnlNative.getLunarEclipseDetail(year, month, day, frameInterval)
    }

    /** 内置世界地图海岸线 (lon, lat) 弧度, 段间 1e7 分隔; 用于绘制食带地图. */
    suspend fun getWorldMapDitu0(): DoubleArray = withContext(Dispatchers.IO) {
        SxwnlNative.getWorldMapDitu0() ?: DoubleArray(0)
    }

    /** ditu1 大图海岸 (4200×2100 像素, 弧度); 段间 1e7 分隔 */
    suspend fun getWorldMapDitu1(): DoubleArray = withContext(Dispatchers.IO) {
        SxwnlNative.getWorldMapData(1) ?: DoubleArray(0)
    }

    /** ditu2 大图国界 */
    suspend fun getWorldMapDitu2(): DoubleArray = withContext(Dispatchers.IO) {
        SxwnlNative.getWorldMapData(2) ?: DoubleArray(0)
    }
}
