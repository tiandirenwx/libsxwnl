#include <jni.h>
#include <string>
#include <cstring>
#include <vector>
#include <mutex>
#include "sxwnl_eclipse_map.h"
#include "sxwnl_capi.h"

#define JNI_PKG "com/sxwnl/calendar/bridge/SxwnlNative"

// ─── Helper: build GeoPoint jobject ───
//
// 全局缓存 EclipseGeoPoint 的 class/ctor. JNI 的 GlobalRef 一旦创建在 VM
// 生命周期内有效, 但 ensureGeoPointClass 在并发首次调用时需要互斥保护,
// 否则两次 NewGlobalRef 会泄漏一个 GlobalRef。

static jclass gGeoPointCls = nullptr;
static jmethodID gGeoPointCtor = nullptr;
static std::once_flag gGeoPointInitFlag;

static void ensureGeoPointClass(JNIEnv *env) {
    std::call_once(gGeoPointInitFlag, [env]() {
        jclass local = env->FindClass("com/sxwnl/calendar/data/EclipseGeoPoint");
        if (local == nullptr) {
            env->ExceptionClear();
            return;
        }
        gGeoPointCls = (jclass)env->NewGlobalRef(local);
        env->DeleteLocalRef(local);
        if (gGeoPointCls != nullptr) {
            gGeoPointCtor = env->GetMethodID(gGeoPointCls, "<init>", "(DDD)V");
        }
    });
}

static jobject newGeoPoint(JNIEnv *env, const SxwnlGeoPoint &p) {
    ensureGeoPointClass(env);
    if (!gGeoPointCls || !gGeoPointCtor) return nullptr;
    return env->NewObject(gGeoPointCls, gGeoPointCtor,
        (jdouble)p.longitude, (jdouble)p.latitude, (jdouble)p.jd);
}

// 大数组拷贝时确保有足够的 local ref 容量, 创建完每个元素后立即释放局部引用
static jobjectArray newGeoPointArray(JNIEnv *env, const SxwnlGeoPoint *pts, int count) {
    ensureGeoPointClass(env);
    if (!gGeoPointCls || !gGeoPointCtor) return nullptr;
    if (count < 0) count = 0;
    // 预留容量: 至少 16, 否则按需 +8
    env->EnsureLocalCapacity(count + 16);
    jobjectArray arr = env->NewObjectArray(count, gGeoPointCls, nullptr);
    if (!arr) return nullptr;
    for (int i = 0; i < count; i++) {
        jobject o = env->NewObject(gGeoPointCls, gGeoPointCtor,
            (jdouble)pts[i].longitude,
            (jdouble)pts[i].latitude,
            (jdouble)pts[i].jd);
        env->SetObjectArrayElement(arr, i, o);
        if (o) env->DeleteLocalRef(o);
    }
    return arr;
}

// ═══ Solar Eclipse Search ═══

extern "C" JNIEXPORT jobjectArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_searchSolarEclipses(
    JNIEnv *env, jclass, jint year, jint month, jint count) {

    SxwnlSolarEclipseList *list = sxwnl_solar_eclipse_search(year, month, count);
    if (!list || list->count == 0) {
        sxwnl_solar_eclipse_list_free(list);
        return nullptr;
    }

    jclass cls = env->FindClass("com/sxwnl/calendar/data/SolarEclipseItem");
    if (!cls) { env->ExceptionClear(); sxwnl_solar_eclipse_list_free(list); return nullptr; }
    jmethodID ctor = env->GetMethodID(cls, "<init>",
        "(DIIIIILjava/lang/String;Ljava/lang/String;DDDDDDIILjava/lang/String;)V");
    if (!ctor) { env->ExceptionClear(); env->DeleteLocalRef(cls); sxwnl_solar_eclipse_list_free(list); return nullptr; }

    // 每条记录: 1 obj + 3 string + cls = 5 ref
    env->EnsureLocalCapacity(list->count * 4 + 16);

    jobjectArray arr = env->NewObjectArray(list->count, cls, nullptr);

    for (int i = 0; i < list->count; i++) {
        auto &item = list->items[i];
        jstring jType  = env->NewStringUTF(item.type);
        jstring jTName = env->NewStringUTF(item.type_name);
        jstring jSeason= env->NewStringUTF(item.season);
        jobject obj = env->NewObject(cls, ctor,
            (jdouble)item.jd,
            (jint)item.year, (jint)item.month, (jint)item.day,
            (jint)item.hour, (jint)item.minute,
            jType, jTName,
            (jdouble)item.gamma, (jdouble)item.magnitude,
            (jdouble)item.center_lon, (jdouble)item.center_lat,
            (jdouble)item.path_width, (jdouble)item.duration,
            (jint)item.saros, (jint)item.saros_member,
            jSeason);
        env->SetObjectArrayElement(arr, i, obj);
        if (obj)    env->DeleteLocalRef(obj);
        if (jType)  env->DeleteLocalRef(jType);
        if (jTName) env->DeleteLocalRef(jTName);
        if (jSeason)env->DeleteLocalRef(jSeason);
    }
    env->DeleteLocalRef(cls);
    sxwnl_solar_eclipse_list_free(list);
    return arr;
}

// ═══ Solar Eclipse Path ═══

extern "C" JNIEXPORT jobject JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_getSolarEclipsePath(
    JNIEnv *env, jclass, jint year, jint month, jint day) {

    SxwnlSolarEclipsePath *path = sxwnl_solar_eclipse_path(year, month, day);
    if (!path) return nullptr;

    jclass cls = env->FindClass("com/sxwnl/calendar/data/SolarEclipsePath");
    if (!cls) { env->ExceptionClear(); sxwnl_solar_eclipse_path_free(path); return nullptr; }
    jmethodID ctor = env->GetMethodID(cls, "<init>",
        "(Ljava/lang/String;Ljava/lang/String;DDDDD"
        "[Lcom/sxwnl/calendar/data/EclipseGeoPoint;"
        "[Lcom/sxwnl/calendar/data/EclipseGeoPoint;"
        "[Lcom/sxwnl/calendar/data/EclipseGeoPoint;"
        "[Lcom/sxwnl/calendar/data/EclipseGeoPoint;"
        "[Lcom/sxwnl/calendar/data/EclipseGeoPoint;"
        "Lcom/sxwnl/calendar/data/EclipseGeoPoint;"
        "Lcom/sxwnl/calendar/data/EclipseGeoPoint;"
        "Lcom/sxwnl/calendar/data/EclipseGeoPoint;"
        "Lcom/sxwnl/calendar/data/EclipseGeoPoint;"
        "Lcom/sxwnl/calendar/data/EclipseGeoPoint;)V");
    if (!ctor) { env->ExceptionClear(); env->DeleteLocalRef(cls); sxwnl_solar_eclipse_path_free(path); return nullptr; }

    // 5 个数组(每个可能 360 帧) + 5 个点 + 2 个 string + cls + obj
    // 数组内引用已在 newGeoPointArray 内即时释放, 这里只需 12 个外层引用
    env->EnsureLocalCapacity(32);

    jstring jType  = env->NewStringUTF(path->type);
    jstring jTName = env->NewStringUTF(path->type_name);
    jobjectArray aCenter = newGeoPointArray(env, path->center_line, path->center_count);
    jobjectArray aUN     = newGeoPointArray(env, path->umbra_north, path->umbra_north_count);
    jobjectArray aUS     = newGeoPointArray(env, path->umbra_south, path->umbra_south_count);
    jobjectArray aPN     = newGeoPointArray(env, path->penumbra_north, path->penumbra_north_count);
    jobjectArray aPS     = newGeoPointArray(env, path->penumbra_south, path->penumbra_south_count);
    jobject pPS = newGeoPoint(env, path->partial_start);
    jobject pCS = newGeoPoint(env, path->central_start);
    jobject pGE = newGeoPoint(env, path->greatest_eclipse);
    jobject pCE = newGeoPoint(env, path->central_end);
    jobject pPE = newGeoPoint(env, path->partial_end);

    jobject obj = env->NewObject(cls, ctor,
        jType, jTName,
        (jdouble)path->gamma, (jdouble)path->magnitude,
        (jdouble)path->path_width_km, (jdouble)path->duration_seconds,
        (jdouble)path->max_eclipse_jd,
        aCenter, aUN, aUS, aPN, aPS,
        pPS, pCS, pGE, pCE, pPE);

    if (jType)  env->DeleteLocalRef(jType);
    if (jTName) env->DeleteLocalRef(jTName);
    if (aCenter) env->DeleteLocalRef(aCenter);
    if (aUN) env->DeleteLocalRef(aUN);
    if (aUS) env->DeleteLocalRef(aUS);
    if (aPN) env->DeleteLocalRef(aPN);
    if (aPS) env->DeleteLocalRef(aPS);
    if (pPS) env->DeleteLocalRef(pPS);
    if (pCS) env->DeleteLocalRef(pCS);
    if (pGE) env->DeleteLocalRef(pGE);
    if (pCE) env->DeleteLocalRef(pCE);
    if (pPE) env->DeleteLocalRef(pPE);
    env->DeleteLocalRef(cls);

    sxwnl_solar_eclipse_path_free(path);
    return obj;
}

// ═══ Solar Eclipse Local ═══

extern "C" JNIEXPORT jobject JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_getLocalSolarEclipse(
    JNIEnv *env, jclass, jint year, jint month, jint day,
    jdouble longitude, jdouble latitude, jint frameInterval) {

    SxwnlLocalSolarEclipse *data = sxwnl_solar_eclipse_local(
        year, month, day, longitude, latitude, frameInterval);
    if (!data) return nullptr;

    // Build frames array
    jclass frameCls = env->FindClass("com/sxwnl/calendar/data/EclipseFrame");
    if (!frameCls) { env->ExceptionClear(); sxwnl_solar_eclipse_local_free(data); return nullptr; }
    jmethodID frameCtor = env->GetMethodID(frameCls, "<init>", "(DDDDDD)V");
    if (!frameCtor) { env->ExceptionClear(); env->DeleteLocalRef(frameCls); sxwnl_solar_eclipse_local_free(data); return nullptr; }

    // 帧数可能上千: +16 容量预留, 每帧及时释放即可
    env->EnsureLocalCapacity(64);

    jobjectArray frames = env->NewObjectArray(data->count, frameCls, nullptr);
    for (int i = 0; i < data->count; i++) {
        auto &f = data->frames[i];
        jobject fo = env->NewObject(frameCls, frameCtor,
            (jdouble)f.sun_radius, (jdouble)f.moon_x, (jdouble)f.moon_y,
            (jdouble)f.moon_radius, (jdouble)f.jd, (jdouble)f.magnitude);
        env->SetObjectArrayElement(frames, i, fo);
        if (fo) env->DeleteLocalRef(fo);
    }

    jclass cls = env->FindClass("com/sxwnl/calendar/data/LocalSolarEclipse");
    if (!cls) { env->ExceptionClear(); env->DeleteLocalRef(frameCls); env->DeleteLocalRef(frames); sxwnl_solar_eclipse_local_free(data); return nullptr; }
    jmethodID ctor = env->GetMethodID(cls, "<init>",
        "(Ljava/lang/String;DDD[Lcom/sxwnl/calendar/data/EclipseFrame;DDDDDDD)V");
    if (!ctor) { env->ExceptionClear(); env->DeleteLocalRef(cls); env->DeleteLocalRef(frameCls); env->DeleteLocalRef(frames); sxwnl_solar_eclipse_local_free(data); return nullptr; }

    jstring jType = env->NewStringUTF(data->type);
    jobject obj = env->NewObject(cls, ctor,
        jType,
        (jdouble)data->max_magnitude, (jdouble)data->moon_sun_ratio,
        (jdouble)data->count,
        frames,
        (jdouble)data->t_c1, (jdouble)data->t_c2,
        (jdouble)data->t_max, (jdouble)data->t_c3,
        (jdouble)data->t_c4, (jdouble)data->t_sunrise,
        (jdouble)data->t_sunset);

    if (jType) env->DeleteLocalRef(jType);
    env->DeleteLocalRef(frames);
    env->DeleteLocalRef(cls);
    env->DeleteLocalRef(frameCls);

    sxwnl_solar_eclipse_local_free(data);
    return obj;
}

// ═══ Lunar Eclipse Search ═══

extern "C" JNIEXPORT jobjectArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_searchLunarEclipses(
    JNIEnv *env, jclass, jint year, jint month, jint count) {

    SxwnlLunarEclipseList *list = sxwnl_lunar_eclipse_search(year, month, count);
    if (!list || list->count == 0) {
        sxwnl_lunar_eclipse_list_free(list);
        return nullptr;
    }

    jclass cls = env->FindClass("com/sxwnl/calendar/data/LunarEclipseItem");
    if (!cls) { env->ExceptionClear(); sxwnl_lunar_eclipse_list_free(list); return nullptr; }
    jmethodID ctor = env->GetMethodID(cls, "<init>",
        "(DIIIIILjava/lang/String;Ljava/lang/String;DIILjava/lang/String;)V");
    if (!ctor) { env->ExceptionClear(); env->DeleteLocalRef(cls); sxwnl_lunar_eclipse_list_free(list); return nullptr; }

    env->EnsureLocalCapacity(list->count * 4 + 16);
    jobjectArray arr = env->NewObjectArray(list->count, cls, nullptr);

    for (int i = 0; i < list->count; i++) {
        auto &item = list->items[i];
        jstring jType  = env->NewStringUTF(item.type);
        jstring jTName = env->NewStringUTF(item.type_name);
        jstring jSeason= env->NewStringUTF(item.season);
        jobject obj = env->NewObject(cls, ctor,
            (jdouble)item.jd,
            (jint)item.year, (jint)item.month, (jint)item.day,
            (jint)item.hour, (jint)item.minute,
            jType, jTName,
            (jdouble)item.magnitude,
            (jint)item.saros, (jint)item.saros_member,
            jSeason);
        env->SetObjectArrayElement(arr, i, obj);
        if (obj)    env->DeleteLocalRef(obj);
        if (jType)  env->DeleteLocalRef(jType);
        if (jTName) env->DeleteLocalRef(jTName);
        if (jSeason)env->DeleteLocalRef(jSeason);
    }
    env->DeleteLocalRef(cls);
    sxwnl_lunar_eclipse_list_free(list);
    return arr;
}

// ═══ Lunar Eclipse Detail ═══

extern "C" JNIEXPORT jobject JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_getLunarEclipseDetail(
    JNIEnv *env, jclass, jint year, jint month, jint day, jint frameInterval) {

    SxwnlLunarEclipseDetail *data = sxwnl_lunar_eclipse_detail(year, month, day, frameInterval);
    if (!data) return nullptr;

    // Build frames array
    jclass frameCls = env->FindClass("com/sxwnl/calendar/data/LunarEclipseFrame");
    if (!frameCls) { env->ExceptionClear(); sxwnl_lunar_eclipse_detail_free(data); return nullptr; }
    jmethodID frameCtor = env->GetMethodID(frameCls, "<init>", "(DDDDDDD)V");
    if (!frameCtor) { env->ExceptionClear(); env->DeleteLocalRef(frameCls); sxwnl_lunar_eclipse_detail_free(data); return nullptr; }

    env->EnsureLocalCapacity(64);
    jobjectArray frames = env->NewObjectArray(data->count, frameCls, nullptr);
    for (int i = 0; i < data->count; i++) {
        auto &f = data->frames[i];
        jobject fo = env->NewObject(frameCls, frameCtor,
            (jdouble)f.moon_x, (jdouble)f.moon_y, (jdouble)f.moon_radius,
            (jdouble)f.umbra_radius, (jdouble)f.penumbra_radius,
            (jdouble)f.jd, (jdouble)f.coverage);
        env->SetObjectArrayElement(frames, i, fo);
        if (fo) env->DeleteLocalRef(fo);
    }

    jclass cls = env->FindClass("com/sxwnl/calendar/data/LunarEclipseDetail");
    if (!cls) { env->ExceptionClear(); env->DeleteLocalRef(frameCls); env->DeleteLocalRef(frames); sxwnl_lunar_eclipse_detail_free(data); return nullptr; }
    jmethodID ctor = env->GetMethodID(cls, "<init>",
        "(Ljava/lang/String;Ljava/lang/String;D"
        "[Lcom/sxwnl/calendar/data/LunarEclipseFrame;"
        "DDDDDDD)V");
    if (!ctor) { env->ExceptionClear(); env->DeleteLocalRef(cls); env->DeleteLocalRef(frameCls); env->DeleteLocalRef(frames); sxwnl_lunar_eclipse_detail_free(data); return nullptr; }

    jstring jType  = env->NewStringUTF(data->type);
    jstring jTName = env->NewStringUTF(data->type_name);
    jobject obj = env->NewObject(cls, ctor,
        jType, jTName,
        (jdouble)data->max_magnitude,
        frames,
        (jdouble)data->t_p1, (jdouble)data->t_u1,
        (jdouble)data->t_u2, (jdouble)data->t_max,
        (jdouble)data->t_u3, (jdouble)data->t_u4,
        (jdouble)data->t_p4);

    if (jType)  env->DeleteLocalRef(jType);
    if (jTName) env->DeleteLocalRef(jTName);
    env->DeleteLocalRef(frames);
    env->DeleteLocalRef(cls);
    env->DeleteLocalRef(frameCls);

    sxwnl_lunar_eclipse_detail_free(data);
    return obj;
}

// ═══ World Map (海岸线轮廓) ═══
//
// 返回 double[]，经纬度(弧度) 交替: [lon0, lat0, lon1, lat1, ...]
// 段间分隔点用 1e7 标记 (Move-To)。
//
// 注: capi 的 size 查询约定: 传 (nullptr, 0) → out_max(0) < sz, 返回 -sz (负数);
// 因此调用方必须取负号. 之前实现遗漏了 ditu0 的负号, 导致 ditu0 永远拿不到.
//
// 防御性兜底:
//   1. `got > need` 不应发生(C 协议保证), 但万一发生会让 SetDoubleArrayRegion
//      越界读 buf -> SIGSEGV. 主动 clamp.
//   2. NewDoubleArray 在内存压力极大时可能返回 nullptr; 不检查会让下一行崩.
extern "C" JNIEXPORT jdoubleArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_getWorldMapDitu0(JNIEnv *env, jclass) {
    int probe = sxwnl_world_map_get_ditu0(nullptr, 0);
    int need = probe < 0 ? -probe : probe;
    if (need <= 0) return nullptr;

    std::vector<double> buf(need);
    int got = sxwnl_world_map_get_ditu0(buf.data(), need);
    if (got <= 0) return nullptr;
    if (got > need) got = need;

    jdoubleArray arr = env->NewDoubleArray(got);
    if (!arr) { env->ExceptionClear(); return nullptr; }
    env->SetDoubleArrayRegion(arr, 0, got, buf.data());
    return arr;
}

extern "C" JNIEXPORT jdoubleArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_getWorldMapData(
    JNIEnv *env, jclass, jint idx) {
    int probe = sxwnl_world_map_get_data(idx, nullptr, 0);
    int need = probe < 0 ? -probe : probe;
    if (need <= 0) return nullptr;

    std::vector<double> buf(need);
    int got = sxwnl_world_map_get_data(idx, buf.data(), need);
    if (got <= 0) return nullptr;
    if (got > need) got = need;

    jdoubleArray arr = env->NewDoubleArray(got);
    if (!arr) { env->ExceptionClear(); return nullptr; }
    env->SetDoubleArrayRegion(arr, 0, got, buf.data());
    return arr;
}
