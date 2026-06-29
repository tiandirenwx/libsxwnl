#include "napi_calendar.h"
#include "napi_utils.h"
#include "sxwnl_capi.h"

static napi_value dayInfoToNapi(napi_env env, const SxwnlDayInfo &d) {
    return NObj(env)
        .i("solarYear",   d.solar_year)   .i("solarMonth",  d.solar_month)  .i("solarDay", d.solar_day)
        .i("lunarYear",   d.lunar_year)    .i("lunarMonth",  d.lunar_month)  .i("lunarDay", d.lunar_day)
        .b("isLeapMonth", d.is_leap_month) .i("weekDay",     d.week_day)
        .i("yearGan",     d.year_gan)      .i("yearZhi",     d.year_zhi)
        .i("monthGan",    d.month_gan)     .i("monthZhi",    d.month_zhi)
        .i("dayGan",      d.day_gan)       .i("dayZhi",      d.day_zhi)
        .i("jieQi",       d.jie_qi)        .i("yueXiang",    d.yue_xiang)
        .i("constellation", d.constellation).i("jian12",      d.jian12)
        .s("yearGZ",         d.year_gz)        .s("monthGZ",     d.month_gz)
        .s("dayGZ",          d.day_gz)         .s("lunarMonthName", d.lunar_month_name)
        .s("lunarDayName",   d.lunar_day_name) .s("jieQiName",   d.jie_qi_name)
        .s("jieQiTime",      d.jie_qi_time)    .s("shengXiao",   d.sheng_xiao)
        .s("constellationName", d.constellation_name)
        .s("weekName",       d.week_name)      .s("yueXiangName", d.yue_xiang_name)
        .s("yueXiangTime",   d.yue_xiang_time)
        .s("jian12Name",     d.jian12_name)
        // 节日/纪年/回历/纳音
        .s("holiday", d.holiday) .s("major", d.major)
        .s("minor",   d.minor)   .s("misc",  d.misc)
        .b("isOffDay", d.is_off_day)
        .i("lunarJunYear",    d.lunar_jun_year)
        .i("lunarLichunYear", d.lunar_lichun_year)
        .i("huangdiYear",     d.huangdi_year)
        .i("moslemYear",  d.moslem_year)
        .i("moslemMonth", d.moslem_month)
        .i("moslemDay",   d.moslem_day)
        .s("yearNaYin",  d.year_nayin)
        .s("monthNaYin", d.month_nayin)
        .s("dayNaYin",   d.day_nayin)
        .d("julianDay",  d.julian_day);
}

napi_value NapiGetDayInfo(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 3);
    SxwnlDayInfo d;
    if (sxwnl_get_day_info(a.intAt(0), a.intAt(1), a.intAt(2), &d) != 0)
        return js_null(env);
    return dayInfoToNapi(env, d);
}

napi_value NapiGetMonthData(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 2);
    SxwnlDayInfo days[31];
    int n = sxwnl_get_month_days(a.intAt(0), a.intAt(1), days, 31);
    if (n < 0) return js_null(env);

    NArr arr(env, n);
    for (int i = 0; i < n; i++) arr.push(dayInfoToNapi(env, days[i]));
    return arr;
}

napi_value NapiLunarToSolar(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 4);
    int oy, om, od;
    if (sxwnl_lunar_to_solar(a.intAt(0), a.intAt(1), a.intAt(2), a.boolAt(3), &oy, &om, &od) != 0)
        return js_null(env);
    return NObj(env).i("year", oy).i("month", om).i("day", od);
}

napi_value NapiSolarToLunar(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 3);
    SxwnlDayInfo d;
    if (sxwnl_solar_to_lunar(a.intAt(0), a.intAt(1), a.intAt(2), &d) != 0)
        return js_null(env);
    return dayInfoToNapi(env, d);
}

napi_value NapiGetJieQiList(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 1);
    SxwnlJieQiItem items[30];
    int n = sxwnl_get_jieqi_list(a.intAt(0), items, 30);
    if (n < 0) return js_null(env);

    NArr arr(env, n);
    for (int i = 0; i < n; i++) {
        arr.push(NObj(env)
            .i("idx", items[i].idx)
            .i("solarMonth", items[i].solar_month)
            .i("solarDay", items[i].solar_day)
            .s("name", items[i].name)
            .s("time", items[i].time));
    }
    return arr;
}

napi_value NapiGetYearLeapMonth(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 1);
    return js_int(env, sxwnl_get_year_leap_month(a.intAt(0)));
}

napi_value NapiGetLunarMonthDays(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 4);
    return js_int(env, sxwnl_get_lunar_month_days(a.intAt(0), a.intAt(1), a.boolAt(2), a.boolAt(3)));
}

napi_value NapiGetSolarMonthValidDays(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 2);
    int days[31];
    int n = sxwnl_get_solar_month_valid_days(a.intAt(0), a.intAt(1), days, 31);
    NArr arr(env, n > 0 ? n : 0);
    for (int i = 0; i < n; i++) arr.push(days[i]);
    return arr;
}

napi_value NapiGetLunarMonths(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 1);
    SxwnlLunarMonth ms[20];
    int n = sxwnl_get_lunar_months(a.intAt(0), ms, 20);
    NArr arr(env, n > 0 ? n : 0);
    for (int i = 0; i < n; i++) {
        arr.push(NObj(env)
            .i("month", ms[i].month)
            .b("isLeap", ms[i].is_leap != 0)
            .b("isSpec", ms[i].is_spec != 0)
            .s("name", ms[i].name));
    }
    return arr;
}

napi_value NapiGetYearCalendar(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 1);
    SxwnlYearCalMonth months[16];
    int n = sxwnl_get_year_calendar(a.intAt(0), months, 16);
    if (n < 0) return js_null(env);

    NArr arr(env, n);
    for (int i = 0; i < n; i++) {
        const SxwnlYearCalMonth &m = months[i];
        NArr jqArr(env, m.jq_count);
        for (int k = 0; k < m.jq_count; k++) {
            const SxwnlYearCalJieQi &jq = m.jieqi[k];
            jqArr.push(NObj(env)
                .i("idx",        jq.idx)
                .s("name",       jq.name)
                .s("gz",         jq.gz)
                .i("solarMonth", jq.solar_month)
                .i("solarDay",   jq.solar_day)
                .s("time",       jq.time)
                .i("dayOffset",  jq.day_offset)
                .s("dayName",    jq.day_name));
        }
        arr.push(NObj(env)
            .i("monthIdx",   m.month_idx)
            .s("monthName",  m.month_name)
            .b("isLeap",     m.is_leap != 0)
            .b("isSpec",     m.is_spec != 0)
            .i("dayCount",   m.day_count)
            .i("solarYear",  m.solar_year)
            .i("solarMonth", m.solar_month)
            .i("solarDay",   m.solar_day)
            .s("shuoGZ",     m.shuo_gz)
            .v("jieQi",      jqArr));
    }
    return arr;
}

// 日月升降/中天 (单天版本, 对应 JS RTS1 月历底栏)
//   参数: (year, month, day, longitude, latitude, tzHours)
napi_value NapiCalcDayRTS(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 6);
    SxwnlDayRTS r{};
    if (sxwnl_calc_day_rts(a.intAt(0), a.intAt(1), a.intAt(2),
                            a.dblAt(3), a.dblAt(4), a.dblAt(5), &r) != 0)
        return js_null(env);

    return NObj(env)
        .s("sunRise",      r.sun_rise)
        .s("sunSet",       r.sun_set)
        .s("sunMeridian",  r.sun_meridian)
        .s("moonRise",     r.moon_rise)
        .s("moonSet",      r.moon_set)
        .s("moonMeridian", r.moon_meridian)
        .s("civilDawn",    r.civil_dawn)
        .s("civilDusk",    r.civil_dusk)
        .s("dayLength",    r.day_length)
        .s("lightLength",  r.light_length);
}
