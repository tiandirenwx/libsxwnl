#include "napi_eclipse_map.h"
#include "napi_utils.h"
#include "sxwnl_eclipse_map.h"
#include "sxwnl_capi.h"
#include <cstring>

// ─── Helper: GeoPoint → napi object ───

static napi_value geoPointToNapi(napi_env env, const SxwnlGeoPoint &p) {
    return NObj(env)
        .d("longitude", p.longitude)
        .d("latitude", p.latitude)
        .d("jd", p.jd);
}

static napi_value geoPointArrayToNapi(napi_env env, const SxwnlGeoPoint *pts, int count) {
    NArr arr(env, count);
    for (int i = 0; i < count; i++)
        arr.push(geoPointToNapi(env, pts[i]));
    return arr;
}

// ═══ Solar Eclipse Search ═══

napi_value NapiSearchSolarEclipses(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 3);
    SxwnlSolarEclipseList *list = sxwnl_solar_eclipse_search(
        a.intAt(0), a.intAt(1), a.intAt(2));
    if (!list) return js_null(env);

    NArr arr(env, list->count);
    for (int i = 0; i < list->count; i++) {
        auto &item = list->items[i];
        arr.push(NObj(env)
            .d("jd", item.jd)
            .i("year", item.year).i("month", item.month).i("day", item.day)
            .i("hour", item.hour).i("minute", item.minute)
            .s("type", item.type).s("typeName", item.type_name)
            .d("gamma", item.gamma).d("magnitude", item.magnitude)
            .d("centerLon", item.center_lon).d("centerLat", item.center_lat)
            .d("pathWidth", item.path_width).d("duration", item.duration)
            .i("saros", item.saros).i("sarosMember", item.saros_member)
            .s("season", item.season));
    }

    sxwnl_solar_eclipse_list_free(list);
    return arr;
}

// ═══ Solar Eclipse Path ═══

napi_value NapiGetSolarEclipsePath(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 3);
    SxwnlSolarEclipsePath *path = sxwnl_solar_eclipse_path(
        a.intAt(0), a.intAt(1), a.intAt(2));
    if (!path) return js_null(env);

    napi_value result = NObj(env)
        .s("type", path->type)
        .s("typeName", path->type_name)
        .d("gamma", path->gamma)
        .d("magnitude", path->magnitude)
        .d("pathWidthKm", path->path_width_km)
        .d("durationSeconds", path->duration_seconds)
        .d("maxEclipseJd", path->max_eclipse_jd)
        .v("centerLine", geoPointArrayToNapi(env, path->center_line, path->center_count))
        .v("umbraNorth", geoPointArrayToNapi(env, path->umbra_north, path->umbra_north_count))
        .v("umbraSouth", geoPointArrayToNapi(env, path->umbra_south, path->umbra_south_count))
        .v("penumbraNorth", geoPointArrayToNapi(env, path->penumbra_north, path->penumbra_north_count))
        .v("penumbraSouth", geoPointArrayToNapi(env, path->penumbra_south, path->penumbra_south_count))
        .v("partialStart", geoPointToNapi(env, path->partial_start))
        .v("centralStart", geoPointToNapi(env, path->central_start))
        .v("greatestEclipse", geoPointToNapi(env, path->greatest_eclipse))
        .v("centralEnd", geoPointToNapi(env, path->central_end))
        .v("partialEnd", geoPointToNapi(env, path->partial_end));

    sxwnl_solar_eclipse_path_free(path);
    return result;
}

// ═══ Solar Eclipse Local ═══

napi_value NapiGetLocalSolarEclipse(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 6);
    SxwnlLocalSolarEclipse *data = sxwnl_solar_eclipse_local(
        a.intAt(0), a.intAt(1), a.intAt(2),
        a.dblAt(3), a.dblAt(4), a.intAt(5));
    if (!data) return js_null(env);

    NArr frames(env, data->count);
    for (int i = 0; i < data->count; i++) {
        auto &f = data->frames[i];
        frames.push(NObj(env)
            .d("sunRadius", f.sun_radius)
            .d("moonX", f.moon_x).d("moonY", f.moon_y)
            .d("moonRadius", f.moon_radius)
            .d("jd", f.jd).d("magnitude", f.magnitude));
    }

    napi_value result = NObj(env)
        .s("type", data->type)
        .d("maxMagnitude", data->max_magnitude)
        .d("moonSunRatio", data->moon_sun_ratio)
        .v("frames", frames)
        .d("tC1", data->t_c1).d("tC2", data->t_c2)
        .d("tMax", data->t_max).d("tC3", data->t_c3)
        .d("tC4", data->t_c4)
        .d("tSunrise", data->t_sunrise).d("tSunset", data->t_sunset);

    sxwnl_solar_eclipse_local_free(data);
    return result;
}

// ═══ Lunar Eclipse Search ═══

napi_value NapiSearchLunarEclipses(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 3);
    SxwnlLunarEclipseList *list = sxwnl_lunar_eclipse_search(
        a.intAt(0), a.intAt(1), a.intAt(2));
    if (!list) return js_null(env);

    NArr arr(env, list->count);
    for (int i = 0; i < list->count; i++) {
        auto &item = list->items[i];
        arr.push(NObj(env)
            .d("jd", item.jd)
            .i("year", item.year).i("month", item.month).i("day", item.day)
            .i("hour", item.hour).i("minute", item.minute)
            .s("type", item.type).s("typeName", item.type_name)
            .d("magnitude", item.magnitude)
            .i("saros", item.saros).i("sarosMember", item.saros_member)
            .s("season", item.season));
    }

    sxwnl_lunar_eclipse_list_free(list);
    return arr;
}

// ═══ Lunar Eclipse Detail ═══

napi_value NapiGetLunarEclipseDetail(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 4);
    SxwnlLunarEclipseDetail *data = sxwnl_lunar_eclipse_detail(
        a.intAt(0), a.intAt(1), a.intAt(2), a.intAt(3));
    if (!data) return js_null(env);

    NArr frames(env, data->count);
    for (int i = 0; i < data->count; i++) {
        auto &f = data->frames[i];
        frames.push(NObj(env)
            .d("moonX", f.moon_x).d("moonY", f.moon_y)
            .d("moonRadius", f.moon_radius)
            .d("umbraRadius", f.umbra_radius)
            .d("penumbraRadius", f.penumbra_radius)
            .d("jd", f.jd).d("coverage", f.coverage));
    }

    napi_value result = NObj(env)
        .s("type", data->type)
        .s("typeName", data->type_name)
        .d("maxMagnitude", data->max_magnitude)
        .v("frames", frames)
        .d("tP1", data->t_p1).d("tU1", data->t_u1)
        .d("tU2", data->t_u2).d("tMax", data->t_max)
        .d("tU3", data->t_u3).d("tU4", data->t_u4)
        .d("tP4", data->t_p4);

    sxwnl_lunar_eclipse_detail_free(data);
    return result;
}

// ═══ World Map (海岸线轮廓) ═══
//
// 返回 number[] 普通数组, 经纬度(弧度) 交替: [lon0, lat0, lon1, lat1, ...]
// 段间分隔点用 1e7 标记 (Move-To)。
//
// 重要: HarmonyOS ArkTS 用的是 Ark 引擎不是 V8, NAPI 用 napi_create_typedarray
// 返回的 Float64Array 在 ArkTS 严格模式下访问 .length 会抛
// "This value does not have a [[TypedArrayName]] internal slot." (在国行真机
// 6.1.1/API 24 上 100% 复现). 改用 napi_create_array + napi_set_element 返回
// 普通 JS number[] 数组, ArkTS 完全兼容.
//
static napi_value make_number_array(napi_env env, const double *data, int n) {
    napi_value arr;
    napi_create_array_with_length(env, (size_t)n, &arr);
    for (int i = 0; i < n; ++i) {
        napi_value v;
        napi_create_double(env, data[i], &v);
        napi_set_element(env, arr, (uint32_t)i, v);
    }
    return arr;
}

napi_value NapiGetWorldMapDitu0(napi_env env, napi_callback_info /*info*/) {
    // 协议: 缓冲不足时返回 -need; sz==0 时返回 0.
    int probe = sxwnl_world_map_get_ditu0(nullptr, 0);
    int need = probe < 0 ? -probe : probe;
    if (need <= 0) return js_null(env);

    std::vector<double> buf(need);
    int got = sxwnl_world_map_get_ditu0(buf.data(), need);
    if (got <= 0) return js_null(env);
    if (got > need) got = need;

    return make_number_array(env, buf.data(), got);
}

// idx=1 → ditu1 大图海岸; idx=2 → ditu2 国界. 首次调用自动加载内建数据.
napi_value NapiGetWorldMapData(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 1);
    int idx = a.intAt(0);
    int probe = sxwnl_world_map_get_data(idx, nullptr, 0);
    int need = probe < 0 ? -probe : probe;
    if (need <= 0) return js_null(env);

    std::vector<double> buf(need);
    int got = sxwnl_world_map_get_data(idx, buf.data(), need);
    if (got <= 0) return js_null(env);
    if (got > need) got = need;

    return make_number_array(env, buf.data(), got);
}
