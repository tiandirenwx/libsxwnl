package com.sxwnl.calendar.data

/**
 * 出生地经纬度数据 (源自寿星天文历 JWv 压缩字符串).
 *
 * 数据规模: 32 个省级行政区, 每个 30-200 个地级市/区县, 合计 ~3500 条
 * 第一级 = 省 / 直辖市 / 特别行政区, 第二级 = 城市 / 区县
 * 八字计算只取经纬度, 时区一律按北京时间 (UT+8) 处理.
 *
 * 与 iOS `Cities.swift` / 鸿蒙 `Cities.ets` 字段结构对齐.
 */
// tz 字段供日食/RTS 等需要观测点时区的模块使用; 八字仅使用经纬度
data class City(val name: String, val lon: Double, val lat: Double, val tz: Double = 8.0)

data class Region(val name: String, val cities: List<City>)

object Cities {

    val REGIONS: List<Region> = buildRegions()

    val DEFAULT: City = REGIONS.first().cities.first()

    // 供日食/RTS 等模块使用的简单扁平预设(取首批 ~30 个 JWv 内置项)
    val PRESET: List<City> = REGIONS.first().cities.take(30)

    private fun buildRegions(): List<Region> {
        val out = ArrayList<Region>(JWV_RAW.size)
        for (raw in JWV_RAW) {
            val parts = raw.split(' ')
            if (parts.size < 2) continue
            val regionName = parts[0]
            val cities = ArrayList<City>(parts.size - 1)
            for (i in 1 until parts.size) {
                val entry = parts[i]
                if (entry.length < 5) continue
                val (lon, lat) = decodeJW(entry.substring(0, 4))
                cities.add(City(entry.substring(4), lon, lat))
            }
            if (cities.isNotEmpty()) out.add(Region(regionName, cities))
        }
        return out
    }

    // 4 字符压缩: [纬度度, 纬度分, 经度度-73, 经度分]
    // '0'-'9' → 0-9, 'A'-'Z' → 10-35, 'a'-'z' → 36-61
    private fun decodeJW(code: String): Pair<Double, Double> {
        val v = IntArray(4)
        for (i in 0..3) {
            var c = code[i].code
            c = when {
                c > 96 -> c - (97 - 36)
                c > 64 -> c - (65 - 10)
                else   -> c - 48
            }
            v[i] = c
        }
        val lon = v[2] + v[3] / 60.0 + 73.0
        val lat = v[0] + v[1] / 60.0
        return lon to lat
    }
}
