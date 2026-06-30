package com.sxwnl.calendar.data

import com.sxwnl.calendar.bridge.SxwnlNative

/**
 *  出生地经纬度目录 — 薄包装 SxwnlNative.geoXxx().
 *
 *  数据来源:  libsxwnl C++ 层的 `GeoPosition` 单例 (src/geo.cpp 中的 `JWv`).
 *
 *  设计:
 *  ─────────────────────────────────────────────────────────────
 *  · 本文件**不再持有任何硬编码城市表**, 跨端唯一数据源 = C++ 核心.
 *  · 公共 API (REGIONS / DEFAULT / PRESET 与旧版兼容) 通过 `by lazy`
 *    在首次访问时一次性加载 — 整个 App 共享同一份结果.
 *  · libsxwnl 是 in-process native, 一次性拉 32 省 × ~100 城市约 50ms,
 *    之后所有访问都是普通对象索引, 不再触发 JNI.
 *  · native 库故障兜底: 返回单条 FALLBACK 城市, 防止下游空数组解引用崩溃.
 *
 *  与 iOS `Cities.swift` / 鸿蒙 `Cities.ets` 字段结构对齐.
 */

// tz 字段供日食/RTS 等需要观测点时区的模块使用; 八字仅使用经纬度.
data class City(val name: String, val lon: Double, val lat: Double, val tz: Double = 8.0)

data class Region(val name: String, val cities: List<City>)

object Cities {

    /** 兜底默认: native 完全不可用时使用. */
    private val FALLBACK = City(name = "天安门", lon = 116.3833, lat = 39.9, tz = 8.0)

    val REGIONS: List<Region> by lazy { loadRegions() }

    val DEFAULT: City by lazy {
        REGIONS.firstOrNull()?.cities?.firstOrNull() ?: FALLBACK
    }

    /** 旧版兼容: 扁平预设 (供日食观测点快速选择条等模块使用), 取首批 30 个. */
    val PRESET: List<City> by lazy {
        REGIONS.firstOrNull()?.cities?.take(30) ?: listOf(FALLBACK)
    }

    private fun loadRegions(): List<Region> {
        // 任何 native 异常 (UnsatisfiedLinkError / 数据损坏等) 都不能让 App 启动崩溃.
        val provinces = try {
            SxwnlNative.geoListProvinces() ?: emptyArray()
        } catch (t: Throwable) {
            return listOf(Region(name = FALLBACK.name, cities = listOf(FALLBACK)))
        }
        val out = ArrayList<Region>(provinces.size)
        for (p in provinces) {
            val raw = try {
                SxwnlNative.geoListCities(p.province) ?: continue
            } catch (t: Throwable) {
                continue
            }
            if (raw.isEmpty()) continue
            val cities = raw.map {
                City(name = it.district, lon = it.longitude,
                     lat = it.latitude, tz = it.timezone)
            }
            out += Region(name = p.province, cities = cities)
        }
        return if (out.isEmpty())
            listOf(Region(name = FALLBACK.name, cities = listOf(FALLBACK)))
        else out
    }
}
