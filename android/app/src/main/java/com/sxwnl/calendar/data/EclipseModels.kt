package com.sxwnl.calendar.data

data class EclipseGeoPoint(
    val longitude: Double,
    val latitude: Double,
    val jd: Double
)

data class SolarEclipseItem(
    val jd: Double,
    val year: Int,
    val month: Int,
    val day: Int,
    val hour: Int,
    val minute: Int,
    val type: String,
    val typeName: String,
    val gamma: Double,
    val magnitude: Double,
    val centerLon: Double,
    val centerLat: Double,
    val pathWidth: Double,
    val duration: Double,
    val saros: Int = 0,
    val sarosMember: Int = 0,
    val season: String = ""
) {
    val hasCenter: Boolean get() = centerLon < 99 && centerLat < 99
}

data class SolarEclipsePath(
    val type: String,
    val typeName: String,
    val gamma: Double,
    val magnitude: Double,
    val pathWidthKm: Double,
    val durationSeconds: Double,
    val maxEclipseJd: Double,
    val centerLine: Array<EclipseGeoPoint>,
    val umbraNorth: Array<EclipseGeoPoint>,
    val umbraSouth: Array<EclipseGeoPoint>,
    val penumbraNorth: Array<EclipseGeoPoint>,
    val penumbraSouth: Array<EclipseGeoPoint>,
    val partialStart: EclipseGeoPoint,
    val centralStart: EclipseGeoPoint,
    val greatestEclipse: EclipseGeoPoint,
    val centralEnd: EclipseGeoPoint,
    val partialEnd: EclipseGeoPoint
)

data class EclipseFrame(
    val sunRadius: Double,
    val moonX: Double,
    val moonY: Double,
    val moonRadius: Double,
    val jd: Double,
    val magnitude: Double
)

data class LocalSolarEclipse(
    val type: String,
    val maxMagnitude: Double,
    val moonSunRatio: Double,
    val frameCount: Double,
    val frames: Array<EclipseFrame>,
    val tC1: Double,
    val tC2: Double,
    val tMax: Double,
    val tC3: Double,
    val tC4: Double,
    val tSunrise: Double,
    val tSunset: Double
)

data class LunarEclipseItem(
    val jd: Double,
    val year: Int,
    val month: Int,
    val day: Int,
    val hour: Int,
    val minute: Int,
    val type: String,
    val typeName: String,
    val magnitude: Double,
    val saros: Int = 0,
    val sarosMember: Int = 0,
    val season: String = ""
)

data class LunarEclipseFrame(
    val moonX: Double,
    val moonY: Double,
    val moonRadius: Double,
    val umbraRadius: Double,
    val penumbraRadius: Double,
    val jd: Double,
    val coverage: Double
)

data class LunarEclipseDetail(
    val type: String,
    val typeName: String,
    val maxMagnitude: Double,
    val frames: Array<LunarEclipseFrame>,
    val tP1: Double,
    val tU1: Double,
    val tU2: Double,
    val tMax: Double,
    val tU3: Double,
    val tU4: Double,
    val tP4: Double
)
