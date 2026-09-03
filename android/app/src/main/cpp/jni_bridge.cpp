// ═══════════════════════════════════════════════════════════════════════
//  jni_bridge.cpp — Kotlin/Java ↔ libsxwnl C API 桥接 (对齐鸿蒙 NAPI)
//
//  设计原则:
//  - DayInfo / BaziColumn / YearCalMonth 等字段众多, 使用 GetFieldID +
//    SetXxxField 方式填充, 避免长签名构造器导致后续扩展困难。
//  - 八字仍使用 handle (Long) 长持有 C++ 对象, 多次 getter 直接读取。
//  - 所有动态字符串调用方负责 sxwnl_string_free。
// ═══════════════════════════════════════════════════════════════════════

#include <jni.h>
#include <string>
#include <cstring>
#include <memory>
#include <vector>
#include "sxwnl_capi.h"

#define PKG "com/sxwnl/calendar/data/"

// ── helpers ─────────────────────────────────────────────────────────────

static jstring newUTF(JNIEnv *env, const char *s) {
    return env->NewStringUTF(s ? s : "");
}

// 设置 Kotlin data class 字段(默认值的 var 字段) - String
static void setStr(JNIEnv *env, jclass cls, jobject obj,
                   const char *field, const char *value) {
    jfieldID f = env->GetFieldID(cls, field, "Ljava/lang/String;");
    if (f == nullptr) { env->ExceptionClear(); return; }
    jstring s = env->NewStringUTF(value ? value : "");
    env->SetObjectField(obj, f, s);
    env->DeleteLocalRef(s);
}
static void setInt(JNIEnv *env, jclass cls, jobject obj,
                   const char *field, jint value) {
    jfieldID f = env->GetFieldID(cls, field, "I");
    if (f == nullptr) { env->ExceptionClear(); return; }
    env->SetIntField(obj, f, value);
}
static void setBool(JNIEnv *env, jclass cls, jobject obj,
                    const char *field, jboolean value) {
    jfieldID f = env->GetFieldID(cls, field, "Z");
    if (f == nullptr) { env->ExceptionClear(); return; }
    env->SetBooleanField(obj, f, value);
}
static void setDouble(JNIEnv *env, jclass cls, jobject obj,
                      const char *field, jdouble value) {
    jfieldID f = env->GetFieldID(cls, field, "D");
    if (f == nullptr) { env->ExceptionClear(); return; }
    env->SetDoubleField(obj, f, value);
}
static void setObj(JNIEnv *env, jclass cls, jobject obj,
                   const char *field, const char *sig, jobject value) {
    jfieldID f = env->GetFieldID(cls, field, sig);
    if (f == nullptr) { env->ExceptionClear(); return; }
    env->SetObjectField(obj, f, value);
}

// 构造空对象 (调用 data class 无参构造器 - 所有字段都有默认值)
static jobject newObj(JNIEnv *env, const char *className, jclass *outCls = nullptr) {
    jclass cls = env->FindClass(className);
    if (cls == nullptr) return nullptr;
    jmethodID ctor = env->GetMethodID(cls, "<init>", "()V");
    if (ctor == nullptr) {
        // data class 可能没有 ()V 默认构造, 改用 newObject 失败后由调用方自行处理
        env->ExceptionClear();
        if (outCls) *outCls = cls;
        return nullptr;
    }
    jobject o = env->NewObject(cls, ctor);
    if (outCls) *outCls = cls;
    return o;
}

// ─── DayInfo 填充 ────────────────────────────────────────────────────────

static jobject newDayInfoObj(JNIEnv *env, const SxwnlDayInfo &d) {
    jclass cls;
    jobject obj = newObj(env, PKG "DayInfo", &cls);
    if (obj == nullptr) return nullptr;

    setInt(env, cls, obj, "solarYear",   d.solar_year);
    setInt(env, cls, obj, "solarMonth",  d.solar_month);
    setInt(env, cls, obj, "solarDay",    d.solar_day);
    setInt(env, cls, obj, "lunarYear",   d.lunar_year);
    setInt(env, cls, obj, "lunarMonth",  d.lunar_month);
    setInt(env, cls, obj, "lunarDay",    d.lunar_day);
    setBool(env, cls, obj, "isLeapMonth", d.is_leap_month);
    setInt(env, cls, obj, "weekDay",     d.week_day);

    setInt(env, cls, obj, "yearGan",  d.year_gan);
    setInt(env, cls, obj, "yearZhi",  d.year_zhi);
    setInt(env, cls, obj, "monthGan", d.month_gan);
    setInt(env, cls, obj, "monthZhi", d.month_zhi);
    setInt(env, cls, obj, "dayGan",   d.day_gan);
    setInt(env, cls, obj, "dayZhi",   d.day_zhi);

    setInt(env, cls, obj, "jieQi",        d.jie_qi);
    setInt(env, cls, obj, "lipuJieQi",    d.lipu_jie_qi);
    setInt(env, cls, obj, "yueXiang",     d.yue_xiang);
    setInt(env, cls, obj, "constellation", d.constellation);
    setInt(env, cls, obj, "jian12",        d.jian12);

    setStr(env, cls, obj, "yearGZ",            d.year_gz);
    setStr(env, cls, obj, "monthGZ",           d.month_gz);
    setStr(env, cls, obj, "dayGZ",             d.day_gz);
    setStr(env, cls, obj, "lunarMonthName",    d.lunar_month_name);
    setStr(env, cls, obj, "lunarDayName",      d.lunar_day_name);
    setStr(env, cls, obj, "jieQiName",         d.jie_qi_name);
    setStr(env, cls, obj, "jieQiTime",         d.jie_qi_time);
    setStr(env, cls, obj, "lipuJieQiName",     d.lipu_jie_qi_name);
    setStr(env, cls, obj, "shengXiao",         d.sheng_xiao);
    setStr(env, cls, obj, "constellationName", d.constellation_name);
    setStr(env, cls, obj, "weekName",          d.week_name);
    setStr(env, cls, obj, "yueXiangName",      d.yue_xiang_name);
    setStr(env, cls, obj, "jian12Name",        d.jian12_name);

    setStr(env, cls, obj, "holiday", d.holiday);
    setStr(env, cls, obj, "major",   d.major);
    setStr(env, cls, obj, "minor",   d.minor);
    setStr(env, cls, obj, "misc",    d.misc);
    setBool(env, cls, obj, "isOffDay", d.is_off_day);

    setInt(env, cls, obj, "lunarJunYear",    d.lunar_jun_year);
    setInt(env, cls, obj, "lunarLichunYear", d.lunar_lichun_year);
    setInt(env, cls, obj, "huangdiYear",     d.huangdi_year);

    setInt(env, cls, obj, "moslemYear",  d.moslem_year);
    setInt(env, cls, obj, "moslemMonth", d.moslem_month);
    setInt(env, cls, obj, "moslemDay",   d.moslem_day);

    setStr(env, cls, obj, "yearNaYin",  d.year_nayin);
    setStr(env, cls, obj, "monthNaYin", d.month_nayin);
    setStr(env, cls, obj, "dayNaYin",   d.day_nayin);

    setStr(env, cls, obj, "yueXiangTime", d.yue_xiang_time);
    setDouble(env, cls, obj, "julianDay", d.julian_day);

    env->DeleteLocalRef(cls);
    return obj;
}

// ═══ Calendar ════════════════════════════════════════════════════════════

extern "C" JNIEXPORT jobject JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_getDayInfo(
    JNIEnv *env, jclass, jint year, jint month, jint day) {
    SxwnlDayInfo d;
    if (sxwnl_get_day_info(year, month, day, &d) != 0) return nullptr;
    return newDayInfoObj(env, d);
}

extern "C" JNIEXPORT jobjectArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_getMonthData(
    JNIEnv *env, jclass, jint year, jint month) {
    SxwnlDayInfo days[31];
    int n = sxwnl_get_month_days(year, month, days, 31);
    if (n <= 0) return nullptr;

    jclass cls = env->FindClass(PKG "DayInfo");
    jobjectArray arr = env->NewObjectArray(n, cls, nullptr);
    for (int i = 0; i < n; i++) {
        jobject o = newDayInfoObj(env, days[i]);
        env->SetObjectArrayElement(arr, i, o);
        env->DeleteLocalRef(o);
    }
    env->DeleteLocalRef(cls);
    return arr;
}

extern "C" JNIEXPORT jintArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_lunarToSolar(
    JNIEnv *env, jclass, jint year, jint month, jint day, jboolean isLeap) {
    int oy, om, od;
    if (sxwnl_lunar_to_solar(year, month, day, isLeap, &oy, &om, &od) != 0)
        return nullptr;
    jintArray arr = env->NewIntArray(3);
    jint buf[] = {oy, om, od};
    env->SetIntArrayRegion(arr, 0, 3, buf);
    return arr;
}

extern "C" JNIEXPORT jobject JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_solarToLunar(
    JNIEnv *env, jclass, jint year, jint month, jint day) {
    SxwnlDayInfo d;
    if (sxwnl_solar_to_lunar(year, month, day, &d) != 0) return nullptr;
    return newDayInfoObj(env, d);
}

extern "C" JNIEXPORT jobjectArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_getJieQiList(
    JNIEnv *env, jclass, jint year) {
    SxwnlJieQiItem items[30];
    int n = sxwnl_get_jieqi_list(year, items, 30);
    if (n <= 0) return nullptr;

    jclass cls = env->FindClass(PKG "JieQiItem");
    jobjectArray arr = env->NewObjectArray(n, cls, nullptr);
    for (int i = 0; i < n; i++) {
        jobject obj = newObj(env, PKG "JieQiItem");
        setInt(env, cls, obj, "idx", items[i].idx);
        setInt(env, cls, obj, "solarMonth", items[i].solar_month);
        setInt(env, cls, obj, "solarDay", items[i].solar_day);
        setStr(env, cls, obj, "name", items[i].name);
        setStr(env, cls, obj, "time", items[i].time);
        env->SetObjectArrayElement(arr, i, obj);
        env->DeleteLocalRef(obj);
    }
    env->DeleteLocalRef(cls);
    return arr;
}

extern "C" JNIEXPORT jint JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_getYearLeapMonth(
    JNIEnv *, jclass, jint year) {
    return sxwnl_get_year_leap_month(year);
}

extern "C" JNIEXPORT jobjectArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_getLunarMonths(
    JNIEnv *env, jclass, jint year) {
    SxwnlLunarMonth lm[16];
    int n = sxwnl_get_lunar_months(year, lm, 16);
    if (n <= 0) return nullptr;

    jclass cls = env->FindClass(PKG "LunarMonth");
    jobjectArray arr = env->NewObjectArray(n, cls, nullptr);
    for (int i = 0; i < n; i++) {
        jobject obj = newObj(env, PKG "LunarMonth");
        setInt(env, cls, obj, "month", lm[i].month);
        setBool(env, cls, obj, "isLeap", lm[i].is_leap != 0);
        setBool(env, cls, obj, "isSpec", lm[i].is_spec != 0);
        setStr(env, cls, obj, "name", lm[i].name);
        env->SetObjectArrayElement(arr, i, obj);
        env->DeleteLocalRef(obj);
    }
    env->DeleteLocalRef(cls);
    return arr;
}

extern "C" JNIEXPORT jint JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_getLunarMonthDays(
    JNIEnv *, jclass, jint year, jint month, jboolean isLeap, jboolean isSpec) {
    return sxwnl_get_lunar_month_days(year, month, isLeap, isSpec);
}

extern "C" JNIEXPORT jstring JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_getLunarDayName(
    JNIEnv *env, jclass, jint day) {
    char buf[16] = {0};
    sxwnl_get_lunar_day_name(day, buf, sizeof(buf));
    return env->NewStringUTF(buf);
}

extern "C" JNIEXPORT jintArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_getSolarMonthValidDays(
    JNIEnv *env, jclass, jint year, jint month) {
    int buf[31];
    int n = sxwnl_get_solar_month_valid_days(year, month, buf, 31);
    if (n <= 0) return nullptr;
    jintArray arr = env->NewIntArray(n);
    env->SetIntArrayRegion(arr, 0, n, buf);
    return arr;
}

// ═══ Year Calendar ═══════════════════════════════════════════════════════

extern "C" JNIEXPORT jobjectArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_getYearCalendar(
    JNIEnv *env, jclass, jint year) {
    SxwnlYearCalMonth months[16];
    int n = sxwnl_get_year_calendar(year, months, 16);
    if (n <= 0) return nullptr;

    jclass monthCls = env->FindClass(PKG "YearCalMonth");
    jclass jqCls = env->FindClass(PKG "YearCalJieQi");
    jobjectArray arr = env->NewObjectArray(n, monthCls, nullptr);

    for (int i = 0; i < n; i++) {
        const auto &m = months[i];
        jobject obj = newObj(env, PKG "YearCalMonth");
        setInt(env, monthCls, obj, "monthIdx", m.month_idx);
        setStr(env, monthCls, obj, "monthName", m.month_name);
        setBool(env, monthCls, obj, "isLeap", m.is_leap != 0);
        setBool(env, monthCls, obj, "isSpec", m.is_spec != 0);
        setInt(env, monthCls, obj, "dayCount", m.day_count);
        setInt(env, monthCls, obj, "solarYear", m.solar_year);
        setInt(env, monthCls, obj, "solarMonth", m.solar_month);
        setInt(env, monthCls, obj, "solarDay", m.solar_day);
        setStr(env, monthCls, obj, "shuoGZ", m.shuo_gz);

        // jieQi 数组
        int jq_n = m.jq_count;
        jobjectArray jqArr = env->NewObjectArray(jq_n, jqCls, nullptr);
        for (int k = 0; k < jq_n; k++) {
            const auto &j = m.jieqi[k];
            jobject jo = newObj(env, PKG "YearCalJieQi");
            setInt(env, jqCls, jo, "idx", j.idx);
            setStr(env, jqCls, jo, "name", j.name);
            setStr(env, jqCls, jo, "gz", j.gz);
            setInt(env, jqCls, jo, "solarMonth", j.solar_month);
            setInt(env, jqCls, jo, "solarDay", j.solar_day);
            setStr(env, jqCls, jo, "time", j.time);
            setInt(env, jqCls, jo, "dayOffset", j.day_offset);
            setStr(env, jqCls, jo, "dayName", j.day_name);
            setInt(env, jqCls, jo, "accMonth", j.acc_month);
            setInt(env, jqCls, jo, "accDay", j.acc_day);
            env->SetObjectArrayElement(jqArr, k, jo);
            env->DeleteLocalRef(jo);
        }
        setObj(env, monthCls, obj, "jieQi",
               "[L" PKG "YearCalJieQi;", jqArr);
        env->DeleteLocalRef(jqArr);

        env->SetObjectArrayElement(arr, i, obj);
        env->DeleteLocalRef(obj);
    }
    env->DeleteLocalRef(monthCls);
    env->DeleteLocalRef(jqCls);
    return arr;
}

// ═══ RTS ═════════════════════════════════════════════════════════════════

extern "C" JNIEXPORT jobject JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_calcDayRTS(
    JNIEnv *env, jclass, jint year, jint month, jint day,
    jdouble longitude, jdouble latitude, jdouble tzHours) {
    SxwnlDayRTS r;
    if (sxwnl_calc_day_rts(year, month, day, longitude, latitude, tzHours, &r) != 0)
        return nullptr;

    jclass cls;
    jobject obj = newObj(env, PKG "DayRTS", &cls);
    if (obj == nullptr) return nullptr;
    setStr(env, cls, obj, "sunRise",     r.sun_rise);
    setStr(env, cls, obj, "sunSet",      r.sun_set);
    setStr(env, cls, obj, "sunMeridian", r.sun_meridian);
    setStr(env, cls, obj, "moonRise",    r.moon_rise);
    setStr(env, cls, obj, "moonSet",     r.moon_set);
    setStr(env, cls, obj, "moonMeridian", r.moon_meridian);
    setStr(env, cls, obj, "civilDawn",   r.civil_dawn);
    setStr(env, cls, obj, "civilDusk",   r.civil_dusk);
    setStr(env, cls, obj, "dayLength",   r.day_length);
    setStr(env, cls, obj, "lightLength", r.light_length);
    env->DeleteLocalRef(cls);
    return obj;
}

// ═══ Bazi ════════════════════════════════════════════════════════════════

extern "C" JNIEXPORT jlong JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_baziCreate(
    JNIEnv *env, jclass,
    jstring name, jboolean gender, jboolean isAst,
    jboolean isLunar, jboolean isLeap, jboolean isSpec,
    jint year, jint month, jint day, jint hour, jint minute, jdouble second,
    jdouble longitude, jdouble latitude, jint lifa) {

    const char *nameStr = env->GetStringUTFChars(name, nullptr);
    SxwnlBaziParams p{};
    strncpy(p.name, nameStr ? nameStr : "", sizeof(p.name) - 1);
    env->ReleaseStringUTFChars(name, nameStr);

    p.gender = gender;
    p.is_ast = isAst;
    p.is_lunar = isLunar;
    p.is_leap = isLeap;
    p.is_spec = isSpec;
    p.year = year; p.month = month; p.day = day;
    p.hour = hour; p.minute = minute; p.second = second;
    p.longitude = longitude; p.latitude = latitude;
    p.lifa = lifa;

    return reinterpret_cast<jlong>(sxwnl_bazi_create(&p));
}

extern "C" JNIEXPORT void JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_baziFree(JNIEnv *, jclass, jlong handle) {
    sxwnl_bazi_free(reinterpret_cast<SxwnlBazi>(handle));
}

#define BAZI_STR_GETTER(JNAME, CFUNC) \
extern "C" JNIEXPORT jstring JNICALL \
Java_com_sxwnl_calendar_bridge_SxwnlNative_##JNAME(JNIEnv *env, jclass, jlong h) { \
    return env->NewStringUTF(CFUNC(reinterpret_cast<SxwnlBazi>(h))); \
}

BAZI_STR_GETTER(baziGetUserName,    sxwnl_bazi_get_user_name)
BAZI_STR_GETTER(baziGetGender,      sxwnl_bazi_get_gender)
BAZI_STR_GETTER(baziGetSolarBirth,  sxwnl_bazi_get_solar_birth)
BAZI_STR_GETTER(baziGetLunarBirth,  sxwnl_bazi_get_lunar_birth)
BAZI_STR_GETTER(baziGetDateOfBirth, sxwnl_bazi_get_date_of_birth)
BAZI_STR_GETTER(baziGetShengXiao,   sxwnl_bazi_get_sheng_xiao)
BAZI_STR_GETTER(baziGetAge,         sxwnl_bazi_get_age)
BAZI_STR_GETTER(baziGetLifa,        sxwnl_bazi_get_lifa)
BAZI_STR_GETTER(baziGetDingQiType,  sxwnl_bazi_get_ding_qi_type)
BAZI_STR_GETTER(baziGetAst,         sxwnl_bazi_get_ast)
BAZI_STR_GETTER(baziGetJieQi,       sxwnl_bazi_get_jie_qi)
BAZI_STR_GETTER(baziGetQiYun,       sxwnl_bazi_get_qi_yun)
BAZI_STR_GETTER(baziGetJiaoYun,     sxwnl_bazi_get_jiao_yun)
BAZI_STR_GETTER(baziGetSiLing,      sxwnl_bazi_get_si_ling)

#undef BAZI_STR_GETTER

// ─── BaziColumn 填充 ─────────────────────────────────────────────────────

static jobject newBaziColumnObj(JNIEnv *env, const SxwnlBaziColumn &c) {
    jclass cls;
    jobject obj = newObj(env, PKG "BaziColumn", &cls);
    if (obj == nullptr) return nullptr;

    setStr(env, cls, obj, "gan", c.gan);
    setStr(env, cls, obj, "zhi", c.zhi);
    setStr(env, cls, obj, "ganShiShen", c.gan_shi_shen);
    setStr(env, cls, obj, "nayin", c.nayin);
    setStr(env, cls, obj, "xingYun", c.xing_yun);
    setStr(env, cls, obj, "ziZuo", c.zi_zuo);
    setStr(env, cls, obj, "kongWang", c.kong_wang);
    setInt(env, cls, obj, "startYear", c.start_year);

    // 藏干
    jclass cgCls = env->FindClass(PKG "CangGanItem");
    jobjectArray cgArr = env->NewObjectArray(c.cang_gan_count, cgCls, nullptr);
    for (int i = 0; i < c.cang_gan_count; i++) {
        jobject cgo = newObj(env, PKG "CangGanItem");
        setStr(env, cgCls, cgo, "gan", c.cang_gan[i]);
        setStr(env, cgCls, cgo, "shiShen", c.cang_gan_shi_shen[i]);
        env->SetObjectArrayElement(cgArr, i, cgo);
        env->DeleteLocalRef(cgo);
    }
    setObj(env, cls, obj, "cangGan",
           "[L" PKG "CangGanItem;", cgArr);
    env->DeleteLocalRef(cgArr);
    env->DeleteLocalRef(cgCls);

    // 神煞
    jclass strCls = env->FindClass("java/lang/String");
    jobjectArray ssArr = env->NewObjectArray(c.shen_sha_count, strCls, nullptr);
    for (int i = 0; i < c.shen_sha_count; i++) {
        jstring s = env->NewStringUTF(c.shen_sha[i]);
        env->SetObjectArrayElement(ssArr, i, s);
        env->DeleteLocalRef(s);
    }
    setObj(env, cls, obj, "shenSha", "[Ljava/lang/String;", ssArr);
    env->DeleteLocalRef(ssArr);
    env->DeleteLocalRef(strCls);

    env->DeleteLocalRef(cls);
    return obj;
}

extern "C" JNIEXPORT jobjectArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_baziGetColumns(JNIEnv *env, jclass, jlong h) {
    SxwnlBaziColumn cols[4];
    int n = sxwnl_bazi_get_columns(reinterpret_cast<SxwnlBazi>(h), cols);
    if (n <= 0) return nullptr;
    jclass cls = env->FindClass(PKG "BaziColumn");
    jobjectArray arr = env->NewObjectArray(n, cls, nullptr);
    for (int i = 0; i < n; i++) {
        jobject o = newBaziColumnObj(env, cols[i]);
        env->SetObjectArrayElement(arr, i, o);
        env->DeleteLocalRef(o);
    }
    env->DeleteLocalRef(cls);
    return arr;
}

extern "C" JNIEXPORT jobject JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_baziGetCurrentDaYun(JNIEnv *env, jclass, jlong h) {
    SxwnlBaziColumn c{};
    if (sxwnl_bazi_get_current_da_yun(reinterpret_cast<SxwnlBazi>(h), &c) != 0) return nullptr;
    return newBaziColumnObj(env, c);
}

extern "C" JNIEXPORT jobject JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_baziGetCurrentLiuNian(JNIEnv *env, jclass, jlong h) {
    SxwnlBaziColumn c{};
    if (sxwnl_bazi_get_current_liu_nian(reinterpret_cast<SxwnlBazi>(h), &c) != 0) return nullptr;
    return newBaziColumnObj(env, c);
}

extern "C" JNIEXPORT jobjectArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_baziGetDaYunColumns(JNIEnv *env, jclass, jlong h) {
    int count = sxwnl_bazi_get_da_yun_count(reinterpret_cast<SxwnlBazi>(h));
    if (count <= 0) return nullptr;
    std::vector<SxwnlBaziColumn> buf(count);
    int n = sxwnl_bazi_get_da_yun_columns(reinterpret_cast<SxwnlBazi>(h), buf.data(), count);
    if (n <= 0) return nullptr;

    jclass cls = env->FindClass(PKG "BaziColumn");
    jobjectArray arr = env->NewObjectArray(n, cls, nullptr);
    for (int i = 0; i < n; i++) {
        jobject o = newBaziColumnObj(env, buf[i]);
        env->SetObjectArrayElement(arr, i, o);
        env->DeleteLocalRef(o);
    }
    env->DeleteLocalRef(cls);
    return arr;
}

extern "C" JNIEXPORT jintArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_baziGetWuXingCount(
    JNIEnv *env, jclass, jlong h, jboolean includeCangGan) {
    int buf[5] = {0};
    sxwnl_bazi_get_wuxing_count(reinterpret_cast<SxwnlBazi>(h), buf, includeCangGan);
    jintArray arr = env->NewIntArray(5);
    env->SetIntArrayRegion(arr, 0, 5, buf);
    return arr;
}

extern "C" JNIEXPORT jobjectArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_baziGetWuXingStatus(JNIEnv *env, jclass, jlong h) {
    char buf[5][8] = {{0}};
    sxwnl_bazi_get_wuxing_status(reinterpret_cast<SxwnlBazi>(h), buf);

    jclass strCls = env->FindClass("java/lang/String");
    jobjectArray arr = env->NewObjectArray(5, strCls, nullptr);
    for (int i = 0; i < 5; i++) {
        jstring s = env->NewStringUTF(buf[i]);
        env->SetObjectArrayElement(arr, i, s);
        env->DeleteLocalRef(s);
    }
    env->DeleteLocalRef(strCls);
    return arr;
}

extern "C" JNIEXPORT jobjectArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_baziGetLiuNian(
    JNIEnv *env, jclass, jlong h, jint startYear, jint count) {
    if (count <= 0) count = 100;
    std::vector<SxwnlLiuNianItem> buf(count);
    int n = sxwnl_bazi_get_liu_nian(reinterpret_cast<SxwnlBazi>(h), startYear,
                                    buf.data(), count);
    if (n <= 0) return nullptr;

    jclass cls = env->FindClass(PKG "LiuNianItem");
    jobjectArray arr = env->NewObjectArray(n, cls, nullptr);
    for (int i = 0; i < n; i++) {
        jobject obj = newObj(env, PKG "LiuNianItem");
        setInt(env, cls, obj, "year", buf[i].year);
        setInt(env, cls, obj, "age", buf[i].age);
        setStr(env, cls, obj, "ganZhi", buf[i].gan_zhi);
        setStr(env, cls, obj, "ganShiShen", buf[i].gan_shi_shen);
        setStr(env, cls, obj, "zhiShiShen", buf[i].zhi_shi_shen);
        setStr(env, cls, obj, "xiaoYun", buf[i].xiao_yun);
        setStr(env, cls, obj, "xiaoYunShiShen", buf[i].xiao_yun_shi_shen);
        env->SetObjectArrayElement(arr, i, obj);
        env->DeleteLocalRef(obj);
    }
    env->DeleteLocalRef(cls);
    return arr;
}

extern "C" JNIEXPORT jobjectArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_baziReverse(
    JNIEnv *env, jclass,
    jint yg, jint yz, jint mg, jint mz, jint dg, jint dz, jint hg, jint hz,
    jint startYear, jint endYear) {
    int sz[8] = {yg, yz, mg, mz, dg, dz, hg, hz};
    SxwnlReverseItem buf[64];
    int n = sxwnl_bazi_reverse(sz, startYear, endYear, buf, 64);
    if (n <= 0) return nullptr;

    jclass cls = env->FindClass(PKG "ReverseItem");
    jobjectArray arr = env->NewObjectArray(n, cls, nullptr);
    for (int i = 0; i < n; i++) {
        jobject obj = newObj(env, PKG "ReverseItem");
        setInt(env, cls, obj, "year",  buf[i].year);
        setInt(env, cls, obj, "month", buf[i].month);
        setInt(env, cls, obj, "day",   buf[i].day);
        setInt(env, cls, obj, "hour",  buf[i].hour);
        setStr(env, cls, obj, "yearGZ",  buf[i].ganzhi[0]);
        setStr(env, cls, obj, "monthGZ", buf[i].ganzhi[1]);
        setStr(env, cls, obj, "dayGZ",   buf[i].ganzhi[2]);
        setStr(env, cls, obj, "hourGZ",  buf[i].ganzhi[3]);
        env->SetObjectArrayElement(arr, i, obj);
        env->DeleteLocalRef(obj);
    }
    env->DeleteLocalRef(cls);
    return arr;
}

// ═══ 纪年字符串工具 ══════════════════════════════════════════════════════

extern "C" JNIEXPORT jint JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_yearStrToAstro(
    JNIEnv *env, jclass, jstring s) {
    const char *str = env->GetStringUTFChars(s, nullptr);
    int32_t r = sxwnl_year_str_to_astro(str);
    env->ReleaseStringUTFChars(s, str);
    return r;
}

extern "C" JNIEXPORT jstring JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_astroYearToStr(
    JNIEnv *env, jclass, jint year, jboolean fullStyle) {
    char buf[64] = {0};
    if (sxwnl_astro_year_to_str(year, fullStyle, buf, 64) == 0) {
        return env->NewStringUTF(buf);
    }
    return env->NewStringUTF("");
}

// ═══ Almanac (老黄历) ════════════════════════════════════════════════════
//
//  注意: SxwnlAlmanac 是大型 POD struct (~13KB), 用栈分配可能触发 JNI 栈
//  限制告警, 改用 std::vector<uint8_t> 堆内存承载. 所有数组按 *_count 严格
//  截断, 防止上层接收到未初始化数据.

static jobjectArray newStringArray(JNIEnv *env, const char *items, int item_count,
                                    size_t item_stride) {
    jclass strCls = env->FindClass("java/lang/String");
    jobjectArray arr = env->NewObjectArray(item_count, strCls, nullptr);
    for (int i = 0; i < item_count; i++) {
        jstring s = env->NewStringUTF(items + i * item_stride);
        env->SetObjectArrayElement(arr, i, s);
        env->DeleteLocalRef(s);
    }
    env->DeleteLocalRef(strCls);
    return arr;
}

extern "C" JNIEXPORT jobject JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_getAlmanac(
    JNIEnv *env, jclass, jint year, jint month, jint day) {
    // 堆上分配避免 ~13KB 在 JNI 栈; 用 unique_ptr 保证: (1) 正确对齐
    //   (alignof(SxwnlAlmanac) 由 new 满足, 而 std::vector<uint8_t> 只保证
    //   alignof(uint8_t)=1, 是潜在 UB)  (2) 零初始化 (避免 *_count 字段
    //   未定义)  (3) 退出时自动释放, 无显式 delete 出错风险.
    auto a = std::make_unique<SxwnlAlmanac>();
    if (sxwnl_get_almanac(year, month, day, a.get()) != 0) return nullptr;

    jclass cls;
    jobject obj = newObj(env, PKG "DayAlmanac", &cls);
    if (obj == nullptr) return nullptr;

    setStr(env, cls, obj, "xiu",       a->xiu);
    setStr(env, cls, obj, "xiuZheng",  a->xiu_zheng);
    setStr(env, cls, obj, "xiuAnimal", a->xiu_animal);
    setStr(env, cls, obj, "xiuLuck",   a->xiu_luck);
    setStr(env, cls, obj, "xiuGong",   a->xiu_gong);
    setStr(env, cls, obj, "twelveGod", a->twelve_god);
    setStr(env, cls, obj, "huangHei",  a->huang_hei);
    setBool(env, cls, obj, "isHuangDao", a->is_huang_dao);
    setStr(env, cls, obj, "chongShengXiao", a->chong_sheng_xiao);
    setStr(env, cls, obj, "chongGanZhi",    a->chong_gan_zhi);
    setStr(env, cls, obj, "sha",            a->sha);
    setStr(env, cls, obj, "xiShenFang",  a->xi_shen_fang);
    setStr(env, cls, obj, "yangGuiFang", a->yang_gui_fang);
    setStr(env, cls, obj, "yinGuiFang",  a->yin_gui_fang);
    setStr(env, cls, obj, "fuShenFang",  a->fu_shen_fang);
    setStr(env, cls, obj, "caiShenFang", a->cai_shen_fang);
    setStr(env, cls, obj, "pengZuGan", a->peng_zu_gan);
    setStr(env, cls, obj, "pengZuZhi", a->peng_zu_zhi);

    // ── quotes ──
    {
        int n = a->quote_count;
        if (n < 0) n = 0;
        if (n > SXWNL_ALMANAC_QUOTE_MAX) n = SXWNL_ALMANAC_QUOTE_MAX;
        jclass qCls = env->FindClass(PKG "AlmanacQuote");
        jobjectArray arr = env->NewObjectArray(n, qCls, nullptr);
        for (int i = 0; i < n; i++) {
            jobject qo = newObj(env, PKG "AlmanacQuote");
            setStr(env, qCls, qo, "source", a->quotes[i].source);
            setStr(env, qCls, qo, "title",  a->quotes[i].title);
            setStr(env, qCls, qo, "luck",   a->quotes[i].luck);
            setStr(env, qCls, qo, "body",   a->quotes[i].body);
            env->SetObjectArrayElement(arr, i, qo);
            env->DeleteLocalRef(qo);
        }
        setObj(env, cls, obj, "quotes", "[L" PKG "AlmanacQuote;", arr);
        env->DeleteLocalRef(arr);
        env->DeleteLocalRef(qCls);
    }

    // ── shenSha ──
    {
        int n = a->shen_sha_count;
        if (n < 0) n = 0;
        if (n > SXWNL_ALMANAC_SHENSHA_MAX) n = SXWNL_ALMANAC_SHENSHA_MAX;
        jclass sCls = env->FindClass(PKG "ShenSha");
        jobjectArray arr = env->NewObjectArray(n, sCls, nullptr);
        for (int i = 0; i < n; i++) {
            jobject so = newObj(env, PKG "ShenSha");
            setStr(env, sCls, so, "name",    a->shen_sha[i].name);
            setBool(env, sCls, so, "isLucky", a->shen_sha[i].is_lucky);
            setInt(env, sCls, so, "weight",  a->shen_sha[i].weight);
            env->SetObjectArrayElement(arr, i, so);
            env->DeleteLocalRef(so);
        }
        setObj(env, cls, obj, "shenSha", "[L" PKG "ShenSha;", arr);
        env->DeleteLocalRef(arr);
        env->DeleteLocalRef(sCls);
    }

    // ── yi / ji (定长字符串数组) ──
    {
        int yn = a->yi_count;
        if (yn < 0) yn = 0;
        if (yn > SXWNL_ALMANAC_TEXT_LIST_ITEM_MAX) yn = SXWNL_ALMANAC_TEXT_LIST_ITEM_MAX;
        jobjectArray yArr = newStringArray(env,
            reinterpret_cast<const char *>(a->yi), yn, SXWNL_ALMANAC_TEXT_LIST_LEN);
        setObj(env, cls, obj, "yi", "[Ljava/lang/String;", yArr);
        env->DeleteLocalRef(yArr);

        int jn = a->ji_count;
        if (jn < 0) jn = 0;
        if (jn > SXWNL_ALMANAC_TEXT_LIST_ITEM_MAX) jn = SXWNL_ALMANAC_TEXT_LIST_ITEM_MAX;
        jobjectArray jArr = newStringArray(env,
            reinterpret_cast<const char *>(a->ji), jn, SXWNL_ALMANAC_TEXT_LIST_LEN);
        setObj(env, cls, obj, "ji", "[Ljava/lang/String;", jArr);
        env->DeleteLocalRef(jArr);
    }

    // ── luckyHours ──
    {
        int n = a->lucky_hour_count;
        if (n < 0) n = 0;
        if (n > SXWNL_ALMANAC_LUCKY_HOUR_MAX) n = SXWNL_ALMANAC_LUCKY_HOUR_MAX;
        jclass lCls = env->FindClass(PKG "LuckyHour");
        jobjectArray arr = env->NewObjectArray(n, lCls, nullptr);
        for (int i = 0; i < n; i++) {
            jobject lo = newObj(env, PKG "LuckyHour");
            setStr(env, lCls, lo, "name", a->lucky_hours[i].name);
            setInt(env, lCls, lo, "zhi", a->lucky_hours[i].zhi);
            env->SetObjectArrayElement(arr, i, lo);
            env->DeleteLocalRef(lo);
        }
        setObj(env, cls, obj, "luckyHours", "[L" PKG "LuckyHour;", arr);
        env->DeleteLocalRef(arr);
        env->DeleteLocalRef(lCls);
    }

    // ── events (用事择吉) ──
    {
        int n = a->event_count;
        if (n < 0) n = 0;
        if (n > SXWNL_ALMANAC_EVENT_MAX) n = SXWNL_ALMANAC_EVENT_MAX;
        jclass eCls = env->FindClass(PKG "EventAdvice");
        jobjectArray arr = env->NewObjectArray(n, eCls, nullptr);
        for (int i = 0; i < n; i++) {
            jobject eo = newObj(env, PKG "EventAdvice");
            setStr(env, eCls, eo, "event",    a->events[i].event);
            setBool(env, eCls, eo, "suitable", a->events[i].suitable);
            setStr(env, eCls, eo, "reason",   a->events[i].reason);
            env->SetObjectArrayElement(arr, i, eo);
            env->DeleteLocalRef(eo);
        }
        setObj(env, cls, obj, "events", "[L" PKG "EventAdvice;", arr);
        env->DeleteLocalRef(arr);
        env->DeleteLocalRef(eCls);
    }

    // ── notes ──
    {
        int n = a->note_count;
        if (n < 0) n = 0;
        if (n > SXWNL_ALMANAC_NOTE_MAX) n = SXWNL_ALMANAC_NOTE_MAX;
        jobjectArray arr = newStringArray(env,
            reinterpret_cast<const char *>(a->notes), n, sizeof(a->notes[0]));
        setObj(env, cls, obj, "notes", "[Ljava/lang/String;", arr);
        env->DeleteLocalRef(arr);
    }

    env->DeleteLocalRef(cls);
    return obj;
}

extern "C" JNIEXPORT jobjectArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_getAlmanacTopics(JNIEnv *env, jclass) {
    // 静态知识 ≤ 64 条; 堆分配 (单条 ~1KB)
    constexpr int kMax = 64;
    std::vector<SxwnlAlmanacTopic> buf(kMax);
    int n = sxwnl_get_almanac_topics(buf.data(), kMax);
    if (n <= 0) return nullptr;
    if (n > kMax) n = kMax;

    jclass cls = env->FindClass(PKG "AlmanacTopic");
    jobjectArray arr = env->NewObjectArray(n, cls, nullptr);
    for (int i = 0; i < n; i++) {
        jobject o = newObj(env, PKG "AlmanacTopic");
        setStr(env, cls, o, "category", buf[i].category);
        setStr(env, cls, o, "title",    buf[i].title);
        setStr(env, cls, o, "body",     buf[i].body);
        env->SetObjectArrayElement(arr, i, o);
        env->DeleteLocalRef(o);
    }
    env->DeleteLocalRef(cls);
    return arr;
}

// ═══ Geo (地理目录) ══════════════════════════════════════════════════════
//
//  与鸿蒙 NAPI 完全对齐: 国内 32 省 + 4000+ 市县 全部驻留在 libsxwnl C++ 层,
//  上层只查目录, 不维护副本.

static jobject newGeoCityObj(JNIEnv *env, const SxwnlGeoCity &c) {
    jclass cls;
    jobject o = newObj(env, PKG "GeoCity", &cls);
    if (o == nullptr) return nullptr;
    setStr(env, cls, o, "province",  c.province);
    setStr(env, cls, o, "district",  c.district);
    setDouble(env, cls, o, "longitude", c.longitude);
    setDouble(env, cls, o, "latitude",  c.latitude);
    setDouble(env, cls, o, "timezone",  c.timezone);
    env->DeleteLocalRef(cls);
    return o;
}

extern "C" JNIEXPORT jobjectArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_geoListProvinces(JNIEnv *env, jclass) {
    constexpr int kMax = 64;
    SxwnlGeoProvince buf[kMax];
    int n = sxwnl_geo_list_provinces(buf, kMax);
    if (n <= 0) return nullptr;
    if (n > kMax) n = kMax;

    jclass cls = env->FindClass(PKG "GeoProvince");
    jobjectArray arr = env->NewObjectArray(n, cls, nullptr);
    for (int i = 0; i < n; i++) {
        jobject o = newObj(env, PKG "GeoProvince");
        setStr(env, cls, o, "province",  buf[i].province);
        setInt(env, cls, o, "cityCount", buf[i].city_count);
        env->SetObjectArrayElement(arr, i, o);
        env->DeleteLocalRef(o);
    }
    env->DeleteLocalRef(cls);
    return arr;
}

extern "C" JNIEXPORT jobjectArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_geoListCities(
    JNIEnv *env, jclass, jstring province) {
    const char *p = env->GetStringUTFChars(province, nullptr);
    // 单省最多约 300 城; 给 512 余量
    constexpr int kMax = 512;
    std::vector<SxwnlGeoCity> buf(kMax);
    int n = sxwnl_geo_list_cities(p ? p : "", buf.data(), kMax);
    env->ReleaseStringUTFChars(province, p);
    if (n <= 0) return nullptr;
    if (n > kMax) n = kMax;

    jclass cls = env->FindClass(PKG "GeoCity");
    jobjectArray arr = env->NewObjectArray(n, cls, nullptr);
    for (int i = 0; i < n; i++) {
        jobject o = newGeoCityObj(env, buf[i]);
        env->SetObjectArrayElement(arr, i, o);
        env->DeleteLocalRef(o);
    }
    env->DeleteLocalRef(cls);
    return arr;
}

extern "C" JNIEXPORT jobjectArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_geoSearch(
    JNIEnv *env, jclass, jstring keyword, jint limit) {
    if (limit <= 0) limit = 64;
    if (limit > 1024) limit = 1024;
    const char *kw = env->GetStringUTFChars(keyword, nullptr);
    std::vector<SxwnlGeoCity> buf(limit);
    int n = sxwnl_geo_search(kw ? kw : "", buf.data(), limit);
    env->ReleaseStringUTFChars(keyword, kw);
    if (n <= 0) return nullptr;
    if (n > limit) n = limit;

    jclass cls = env->FindClass(PKG "GeoCity");
    jobjectArray arr = env->NewObjectArray(n, cls, nullptr);
    for (int i = 0; i < n; i++) {
        jobject o = newGeoCityObj(env, buf[i]);
        env->SetObjectArrayElement(arr, i, o);
        env->DeleteLocalRef(o);
    }
    env->DeleteLocalRef(cls);
    return arr;
}

extern "C" JNIEXPORT jobjectArray JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_geoListTimezones(JNIEnv *env, jclass) {
    constexpr int kMax = 128;
    std::vector<SxwnlTimezoneGroup> buf(kMax);
    int n = sxwnl_geo_list_timezone_groups(buf.data(), kMax);
    if (n <= 0) return nullptr;
    if (n > kMax) n = kMax;

    jclass cls = env->FindClass(PKG "TimezoneGroup");
    jobjectArray arr = env->NewObjectArray(n, cls, nullptr);
    for (int i = 0; i < n; i++) {
        jobject o = newObj(env, PKG "TimezoneGroup");
        setStr(env, cls, o, "continent", buf[i].continent);
        setStr(env, cls, o, "country",   buf[i].country);
        setDouble(env, cls, o, "timezone", buf[i].timezone);

        int cc = buf[i].city_count;
        if (cc < 0) cc = 0;
        if (cc > 8) cc = 8;
        jobjectArray cArr = newStringArray(env,
            reinterpret_cast<const char *>(buf[i].cities), cc, sizeof(buf[i].cities[0]));
        setObj(env, cls, o, "cities", "[Ljava/lang/String;", cArr);
        env->DeleteLocalRef(cArr);

        env->SetObjectArrayElement(arr, i, o);
        env->DeleteLocalRef(o);
    }
    env->DeleteLocalRef(cls);
    return arr;
}

extern "C" JNIEXPORT jobject JNICALL
Java_com_sxwnl_calendar_bridge_SxwnlNative_geoDefault(JNIEnv *env, jclass) {
    SxwnlGeoCity c{};
    if (sxwnl_geo_default(&c) != 0) return nullptr;
    return newGeoCityObj(env, c);
}
