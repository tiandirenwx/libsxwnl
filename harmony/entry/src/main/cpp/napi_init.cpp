#include <napi/native_api.h>
#include "napi_utils.h"
#include "napi_calendar.h"
#include "napi_bazi.h"
#include "napi_eclipse_map.h"

static napi_value Init(napi_env env, napi_value exports) {
    static const NExport entries[] = {
        // Calendar
        {"getDayInfo",        NapiGetDayInfo},
        {"getMonthData",      NapiGetMonthData},
        {"lunarToSolar",      NapiLunarToSolar},
        {"solarToLunar",      NapiSolarToLunar},
        {"getJieQiList",      NapiGetJieQiList},
        {"getYearLeapMonth",  NapiGetYearLeapMonth},
        {"getLunarMonths",    NapiGetLunarMonths},
        {"getLunarMonthDays", NapiGetLunarMonthDays},
        {"getLunarDayName",   NapiGetLunarDayName},
        {"getSolarMonthValidDays", NapiGetSolarMonthValidDays},
        {"getYearCalendar",   NapiGetYearCalendar},
        {"calcDayRTS",        NapiCalcDayRTS},
        {"getAlmanac",        NapiGetAlmanac},
        {"getAlmanacTopics",  NapiGetAlmanacTopics},
        // Geo (城市目录)
        {"geoListProvinces",  NapiGeoListProvinces},
        {"geoListCities",     NapiGeoListCities},
        {"geoSearch",         NapiGeoSearch},
        {"geoListTimezones",  NapiGeoListTimezones},
        {"geoDefault",        NapiGeoDefault},
        // Bazi
        {"calcBazi",          NapiCalcBazi},
        {"baziReverse",       NapiBaziReverse},
        // Eclipse Map (structured data for visualization)
        {"searchSolarEclipses",   NapiSearchSolarEclipses},
        {"getSolarEclipsePath",   NapiGetSolarEclipsePath},
        {"getLocalSolarEclipse",  NapiGetLocalSolarEclipse},
        {"searchLunarEclipses",   NapiSearchLunarEclipses},
        {"getLunarEclipseDetail",  NapiGetLunarEclipseDetail},
        // World Map (海岸线轮廓, 用于绘制食带地图)
        {"getWorldMapDitu0",       NapiGetWorldMapDitu0},
        {"getWorldMapData",        NapiGetWorldMapData},
    };
    napi_export_all(env, exports, entries, sizeof(entries) / sizeof(entries[0]));
    return exports;
}

NAPI_MODULE(sxwnl_napi, Init)
