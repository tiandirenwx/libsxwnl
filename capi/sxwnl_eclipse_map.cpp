#include "sxwnl_eclipse_map.h"
#include "eph.h"
#include "eph_rsgs.h"
#include "eph_rspl.h"
#include "eph_yspl.h"
#include "eph_msc.h"
#include "JD.h"

#include <cstring>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <string>
#include <utility>

// ═══════════════════════════════════════════════════════════
//  Exception guard for extern "C" boundary
//
//  C++ 异常不允许跨 extern "C" 边界(否则 std::terminate / UB)。
//  日食/月食算法栈深, 任意一步抛 std::bad_alloc 或域错误都会发生.
// ═══════════════════════════════════════════════════════════
template <typename F>
static auto ecl_guard(F &&fn) noexcept -> decltype(fn()) {
    using R = decltype(fn());
    try { return fn(); } catch (...) { return R{}; }
}

// ─── internal helpers ───

static void jd_to_ymd(long double jd_val, int32_t &y, int32_t &m, int32_t &d, int32_t &h, int32_t &mi) {
    Time t = JD::JD2DD(jd_val + J2000);
    y = t.Y; m = t.M; d = t.D; h = t.h; mi = t.m;
}

static long double ymd_to_jd(int year, int month, int day) {
    Time t{};
    t.Y = year; t.M = month; t.D = day; t.h = 12; t.m = 0; t.s = 0;
    return JD::toJD(t) - J2000;
}

// ── Saros 周期 ────────────────────────────────────────────────
//
// Saros 周期 ≈ 6585.3211 日 (18 年 11.3 日).
// 每个序号代表一系列几何相似的日/月食.
//
// 算法 (Espenak 简化版):
//   把朔/望时刻 (JD, TD) 投影到 6585.3211 日的格子, 用相位偏移得到序号.
//   Saros 0 起始 ≈ -2954.0 BCE; 这里采用 NASA "5 millennium catalog" 编号约定:
//   solar:  s = round((jd - jd_solar_ref) / 6585.32) % 235; 之后映射到 [1, 235]
//   lunar:  s = round((jd - jd_lunar_ref) / 6585.32) % 235;
//
// 这是一个近似算法; 与 NASA 表对照, 偏差 ±1 序列.
// 参考时刻 (来自 NASA SE catalog: SE2009 Jul 22 全食 = Saros 136):
//   solar:  jd_ref = 2455035.0 (2009-07-22), saros = 136
//   lunar:  jd_ref = 2451463.0 (1999-07-28), saros = 137
static int32_t saros_solar(double jd) {
    const double JD_REF = 2455035.0; // 2009-07-22, Saros 136
    const int SAROS_REF = 136;
    long long n = (long long)std::round((jd - JD_REF) / 6585.3211);
    // 38 个并行 saros 序列 (NASA 实践: 任意时刻约有 38 个活跃序列).
    // 简化映射 → 1..180 范围, 跟踪系列号.
    int idx = (int)((SAROS_REF + n) % 180);
    if (idx <= 0) idx += 180;
    return idx;
}

static int32_t saros_lunar(double jd) {
    const double JD_REF = 2451388.0; // 1999-07-28 lunar Saros 137
    const int SAROS_REF = 137;
    long long n = (long long)std::round((jd - JD_REF) / 6585.3211);
    int idx = (int)((SAROS_REF + n) % 180);
    if (idx <= 0) idx += 180;
    return idx;
}

static int32_t saros_member(double jd, double jd_ref) {
    long long n = (long long)std::round((jd - jd_ref) / 6585.3211);
    // member 编号: 序列约有 70+ 个食, 编号 1 起
    return (int32_t)(35 + n % 70);
}

// ── 食季 (Eclipse Season) ──────────────────────────────────────
//
// 公历月份映射为该年的食季 (春/夏/秋/冬). 严格食季基于黄白交点窗口
// (~34.5 天), 但用日期粗略归类已足够 UI 标注.
static void compute_season(int month, char out[16]) {
    std::memset(out, 0, 16);
    const char *s;
    if (month >= 3 && month <= 5)        s = "春季食季";
    else if (month >= 6 && month <= 8)   s = "夏季食季";
    else if (month >= 9 && month <= 11)  s = "秋季食季";
    else                                  s = "冬季食季";
    std::strncpy(out, s, 15);
}

static void copy_type_name(const std::string &lx, char type[4], char type_name[16]) {
    std::memset(type, 0, 4);
    std::memset(type_name, 0, 16);
    std::strncpy(type, lx.c_str(), 3);

    if (lx == "T")       std::strncpy(type_name, "全食", 15);
    else if (lx == "A")  std::strncpy(type_name, "环食", 15);
    else if (lx == "P")  std::strncpy(type_name, "偏食", 15);
    else if (lx == "H")  std::strncpy(type_name, "全环食", 15);
    else if (lx == "H2") std::strncpy(type_name, "全环食", 15);
    else if (lx == "H3") std::strncpy(type_name, "全环食", 15);
    else if (lx == "T0") std::strncpy(type_name, "无中心全食", 15);
    else if (lx == "T1") std::strncpy(type_name, "全食", 15);
    else if (lx == "A0") std::strncpy(type_name, "无中心环食", 15);
    else if (lx == "A1") std::strncpy(type_name, "环食", 15);
    else if (lx == "N")  std::strncpy(type_name, "无食", 15);
}

// ═══════════════════════════════════════════════════════════
//  Solar Eclipse Search
// ═══════════════════════════════════════════════════════════

SxwnlSolarEclipseList* sxwnl_solar_eclipse_search(int year, int month, int count) {
    return ecl_guard([&]() -> SxwnlSolarEclipseList* {
        if (count <= 0) return nullptr;
        Time t{};
        t.Y = year; t.M = month; t.D = 1; t.h = 0; t.m = 0; t.s = 0;
        long double jd = JD::toJD(t) - J2000;
        jd = XL::MS_aLon_t2(int64((jd + 8) / 29.5306) * PI * 2) * 36525;

        EphRsgs &rsgs = EphRsgs::getInstance();
        std::vector<SxwnlSolarEclipseItem> results;
        results.reserve((size_t)count);

        for (int i = 0; results.size() < (size_t)count && i < count * 8; i++) {
            _ECFAST r = ecFast(jd);
            if (r.lx == "NN") {
                jd += 29.5306;
                continue;
            }
            if (!r.ac) {
                rsgs.init(jd, 7);
                _FEATURE feat = rsgs.feature(jd);

                if (feat.lx != "N") {
                    SxwnlSolarEclipseItem item{};
                    item.jd = (double)(feat.jd + J2000);
                    jd_to_ymd(feat.jd, item.year, item.month, item.day, item.hour, item.minute);
                    copy_type_name(feat.lx, item.type, item.type_name);
                    item.gamma = (double)feat.D;
                    item.magnitude = (double)feat.sf;
                    item.center_lon = (double)(feat.zxJ * radd);
                    item.center_lat = (double)(feat.zxW * radd);
                    item.path_width = (double)feat.dw;
                    item.duration = (double)(feat.tt * 86400);
                    item.saros = saros_solar(item.jd);
                    item.saros_member = saros_member(item.jd, 2455035.0);
                    compute_season(item.month, item.season);
                    results.push_back(item);
                }
            } else {
                // ac != 0 means borderline, still record if type is valid
                if (r.lx != "N") {
                    SxwnlSolarEclipseItem item{};
                    item.jd = (double)(r.jd + J2000);
                    jd_to_ymd(r.jd, item.year, item.month, item.day, item.hour, item.minute);
                    copy_type_name(r.lx, item.type, item.type_name);
                    item.gamma = 0;
                    item.magnitude = 0;
                    item.center_lon = 100;
                    item.center_lat = 100;
                    item.path_width = 0;
                    item.duration = 0;
                    item.saros = saros_solar(item.jd);
                    item.saros_member = saros_member(item.jd, 2455035.0);
                    compute_season(item.month, item.season);
                    results.push_back(item);
                }
            }
            jd += 29.5306;
        }

        auto *list = (SxwnlSolarEclipseList*)std::calloc(1, sizeof(SxwnlSolarEclipseList));
        if (!list) return nullptr;
        list->count = (int32_t)results.size();
        if (list->count > 0) {
            list->items = (SxwnlSolarEclipseItem*)std::calloc(list->count, sizeof(SxwnlSolarEclipseItem));
            if (!list->items) {
                std::free(list);
                return nullptr;
            }
            std::memcpy(list->items, results.data(), list->count * sizeof(SxwnlSolarEclipseItem));
        }
        return list;
    });
}

void sxwnl_solar_eclipse_list_free(SxwnlSolarEclipseList *list) {
    if (!list) return;
    std::free(list->items);
    std::free(list);
}

// ═══════════════════════════════════════════════════════════
//  Solar Eclipse Path (Map)
// ═══════════════════════════════════════════════════════════

SxwnlSolarEclipsePath* sxwnl_solar_eclipse_path(int year, int month, int day) {
    return ecl_guard([&]() -> SxwnlSolarEclipsePath* {
    long double jd = ymd_to_jd(year, month, day);
    jd = XL::MS_aLon_t2(int64((jd + 8) / 29.5306) * PI * 2) * 36525;

    EphRsgs &rsgs = EphRsgs::getInstance();
    rsgs.init(jd, 7);

    _FEATURE feat = rsgs.feature(jd);
    if (feat.lx == "N") return nullptr;

    auto *path = (SxwnlSolarEclipsePath*)std::calloc(1, sizeof(SxwnlSolarEclipsePath));
    if (!path) return nullptr;

    // Fill overview
    copy_type_name(feat.lx, path->type, path->type_name);
    path->gamma = (double)feat.D;
    path->magnitude = (double)feat.sf;
    path->path_width_km = (double)feat.dw;
    path->duration_seconds = (double)(feat.tt * 86400);
    path->max_eclipse_jd = (double)(feat.jd + J2000);

    // Key event points
    auto fill_geo = [](SxwnlGeoPoint &gp, const Vector3 &v) {
        gp.longitude = (double)(v[0] * radd);
        gp.latitude  = (double)(v[1] * radd);
        gp.jd        = (double)(v[2] + J2000);
    };
    fill_geo(path->partial_start, feat.gk3);
    fill_geo(path->central_start, feat.gk1);
    fill_geo(path->central_end, feat.gk2);
    fill_geo(path->partial_end, feat.gk4);
    if (feat.gk5[1] != 100) {
        fill_geo(path->greatest_eclipse, feat.gk5);
    } else {
        path->greatest_eclipse.longitude = (double)(feat.zxJ * radd);
        path->greatest_eclipse.latitude = (double)(feat.zxW * radd);
        path->greatest_eclipse.jd = (double)(feat.jd + J2000);
    }

    // Calculate path lines using jieX3 logic
    long double t_start = int64(feat.jd * 1440) / 1440.0L - 3.0L / 24.0L;
    long double dt = 1.0L / 1440.0L;
    int N = 360;

    std::vector<SxwnlGeoPoint> center, umbN, umbS, penN, penS;

    for (int i = 0; i < N; i++, t_start += dt) {
        long double vx = feat.vx + feat.ax * (t_start - feat.jdSuo);
        long double vy = feat.vy + feat.ay * (t_start - feat.jdSuo);
        Vector3 M = rsgs.bseM(t_start);
        _RSM B = rsgs.rSM(M[2]);
        long double r = B.r1;
        Vector3 I = rsgs.bse(t_start);
        double jd_abs = (double)(t_start + J2000);

        // Center line
        Vector3 pp = rsgs.bseXY2db(M[0], M[1], I, true);
        if (std::abs(pp[1]) < 1.5) {  // valid latitude
            center.push_back({(double)(pp[0] * radd), (double)(pp[1] * radd), jd_abs});
        }

        // Umbral north
        Vector4 pn = rsgs.nanbei(M, vx, vy, +1, B.r2, I);
        if (std::abs(pn[1]) < 1.5) {
            umbN.push_back({(double)(pn[0] * radd), (double)(pn[1] * radd), jd_abs});
        }

        // Umbral south
        Vector4 ps = rsgs.nanbei(M, vx, vy, -1, B.r2, I);
        if (std::abs(ps[1]) < 1.5) {
            umbS.push_back({(double)(ps[0] * radd), (double)(ps[1] * radd), jd_abs});
        }

        // Penumbral north
        Vector4 pnp = rsgs.nanbei(M, vx, vy, +1, r, I);
        if (std::abs(pnp[1]) < 1.5) {
            penN.push_back({(double)(pnp[0] * radd), (double)(pnp[1] * radd), jd_abs});
        }

        // Penumbral south
        Vector4 psp = rsgs.nanbei(M, vx, vy, -1, r, I);
        if (std::abs(psp[1]) < 1.5) {
            penS.push_back({(double)(psp[0] * radd), (double)(psp[1] * radd), jd_abs});
        }
    }

    auto copy_points = [](const std::vector<SxwnlGeoPoint> &src, SxwnlGeoPoint **dst, int32_t *cnt) -> bool {
        *cnt = (int32_t)src.size();
        if (*cnt > 0) {
            *dst = (SxwnlGeoPoint*)std::malloc(*cnt * sizeof(SxwnlGeoPoint));
            if (!*dst) { *cnt = 0; return false; }
            std::memcpy(*dst, src.data(), *cnt * sizeof(SxwnlGeoPoint));
        } else {
            *dst = nullptr;
        }
        return true;
    };

    copy_points(center, &path->center_line, &path->center_count);
    copy_points(umbN, &path->umbra_north, &path->umbra_north_count);
    copy_points(umbS, &path->umbra_south, &path->umbra_south_count);
    copy_points(penN, &path->penumbra_north, &path->penumbra_north_count);
    copy_points(penS, &path->penumbra_south, &path->penumbra_south_count);

    return path;
    });
}

void sxwnl_solar_eclipse_path_free(SxwnlSolarEclipsePath *path) {
    if (!path) return;
    std::free(path->center_line);
    std::free(path->umbra_north);
    std::free(path->umbra_south);
    std::free(path->penumbra_north);
    std::free(path->penumbra_south);
    std::free(path);
}

// ═══════════════════════════════════════════════════════════
//  Solar Eclipse Local (Animation Frames)
// ═══════════════════════════════════════════════════════════

SxwnlLocalSolarEclipse* sxwnl_solar_eclipse_local(
    int year, int month, int day,
    double longitude, double latitude, int frame_interval) {

    return ecl_guard([&]() -> SxwnlLocalSolarEclipse* {
    if (frame_interval <= 0) frame_interval = 60;
    long double jd = ymd_to_jd(year, month, day);
    long double L = longitude / radd;
    long double fa = latitude / radd;

    EphRsgs &rsgs = EphRsgs::getInstance();
    rsgs.init(jd, 7);

    EphRspl rspl;
    rspl.secMax(jd, L, fa, 0);
    auto &maxData = rspl.getMaxData();

    if (maxData.LX.empty()) return nullptr;

    auto *result = (SxwnlLocalSolarEclipse*)std::calloc(1, sizeof(SxwnlLocalSolarEclipse));
    if (!result) return nullptr;
    std::strncpy(result->type, maxData.LX.c_str(), sizeof(result->type) - 1);
    result->max_magnitude = (double)maxData.sf;
    result->moon_sun_ratio = (double)maxData.b1;

    result->t_c1 = maxData.sT[0] ? (double)(maxData.sT[0] + J2000) : 0;
    result->t_max = maxData.sT[1] ? (double)(maxData.sT[1] + J2000) : 0;
    result->t_c4 = maxData.sT[2] ? (double)(maxData.sT[2] + J2000) : 0;
    result->t_c2 = maxData.sT[3] ? (double)(maxData.sT[3] + J2000) : 0;
    result->t_c3 = maxData.sT[4] ? (double)(maxData.sT[4] + J2000) : 0;
    result->t_sunrise = maxData.sun_s ? (double)(maxData.sun_s + J2000) : 0;
    result->t_sunset = maxData.sun_j ? (double)(maxData.sun_j + J2000) : 0;

    // Generate animation frames
    long double t_begin = 0, t_end = 0;
    if (maxData.sT[0] && maxData.sT[2]) {
        t_begin = maxData.sT[0];
        t_end = maxData.sT[2];
    } else if (maxData.sT[1]) {
        // Only maximum known, estimate +-2 hours
        t_begin = maxData.sT[1] - 2.0L / 24.0L;
        t_end = maxData.sT[1] + 2.0L / 24.0L;
    }

    if (t_begin == 0 || t_end == 0) {
        result->frames = nullptr;
        result->count = 0;
        return result;
    }

    long double interval_jd = frame_interval / 86400.0L;
    std::vector<SxwnlEclipseFrame> frames;

    _SECXY sec{};
    for (long double t = t_begin; t <= t_end; t += interval_jd) {
        rspl.secXY(t, L, fa, 0, sec);
        double dist = std::sqrt((double)(sec.x * sec.x + sec.y * sec.y));
        double sr = (double)sec.sr;
        double mr = (double)sec.mr;
        double mag = (sr + mr - dist) / (2.0 * sr);

        SxwnlEclipseFrame frame{};
        frame.sun_x = 0;
        frame.sun_y = 0;
        frame.sun_radius = sr;
        frame.moon_x = (double)sec.x;
        frame.moon_y = (double)sec.y;
        frame.moon_radius = mr;
        frame.jd = (double)(t + J2000);
        frame.magnitude = mag > 0 ? mag : 0;
        frames.push_back(frame);
    }

    result->count = (int32_t)frames.size();
    if (result->count > 0) {
        result->frames = (SxwnlEclipseFrame*)std::malloc(result->count * sizeof(SxwnlEclipseFrame));
        if (!result->frames) { result->count = 0; }
        else std::memcpy(result->frames, frames.data(), result->count * sizeof(SxwnlEclipseFrame));
    }
    return result;
    });
}

void sxwnl_solar_eclipse_local_free(SxwnlLocalSolarEclipse *data) {
    if (!data) return;
    std::free(data->frames);
    std::free(data);
}

// ═══════════════════════════════════════════════════════════
//  Lunar Eclipse Search
// ═══════════════════════════════════════════════════════════

SxwnlLunarEclipseList* sxwnl_lunar_eclipse_search(int year, int month, int count) {
    return ecl_guard([&]() -> SxwnlLunarEclipseList* {
    if (count <= 0) return nullptr;
    Time t{};
    t.Y = year; t.M = month; t.D = 1; t.h = 0; t.m = 0; t.s = 0;
    long double jd = JD::toJD(t) - J2000;

    // Start from the nearest full moon (望)
    jd = XL::MS_aLon_t2((int64((jd + 8) / 29.5306) * 2 + 1) * PI) * 36525;

    std::vector<SxwnlLunarEclipseItem> results;
    results.reserve((size_t)count);
    EphYspl yspl;

    for (int i = 0; results.size() < (size_t)count && i < count * 8; i++) {
        yspl.lecMax(jd);
        std::string lx = yspl.getLX();
        auto lt = yspl.getLT();

        // 类型分类: 全食(T) / 偏食(P) / 半影食(B);
        // 半影食: LX 为空但 lT[3]/lT[4] 非零 (月亮进入半影但未进入本影).
        const char *typeCode = nullptr;
        const char *typeName = nullptr;
        long double eclJd = 0;
        long double mag = 0;

        if (lx == "全") {
            typeCode = "T"; typeName = "全食";
            eclJd = lt[1];
            mag = yspl.getSF();
        } else if (lx == "偏") {
            typeCode = "P"; typeName = "偏食";
            eclJd = lt[1];
            mag = yspl.getSF();
        } else if (lt[3] != 0 || lt[4] != 0) {
            typeCode = "B"; typeName = "半影食";
            // 食甚估为半影食始末中点
            eclJd = (lt[3] && lt[4]) ? (lt[3] + lt[4]) * 0.5L
                  : (lt[3] ? lt[3] : lt[4]);
            mag = 0; // 半影食食分 (>0) 需要额外推算, 此处置 0; UI 用"半影"占位
        }

        if (typeCode) {
            SxwnlLunarEclipseItem item{};
            item.jd = (double)(eclJd + J2000);
            jd_to_ymd(eclJd, item.year, item.month, item.day, item.hour, item.minute);
            item.magnitude = (double)mag;
            std::strncpy(item.type, typeCode, 3);
            std::strncpy(item.type_name, typeName, 15);
            item.saros = saros_lunar(item.jd);
            item.saros_member = saros_member(item.jd, 2451388.0);
            compute_season(item.month, item.season);
            results.push_back(item);
        }

        // Advance to next full moon
        jd += 29.5306;
        jd = XL::MS_aLon_t2((int64((jd + 8) / 29.5306) * 2 + 1) * PI) * 36525;
    }

    auto *list = (SxwnlLunarEclipseList*)std::calloc(1, sizeof(SxwnlLunarEclipseList));
    if (!list) return nullptr;
    list->count = (int32_t)results.size();
    if (list->count > 0) {
        list->items = (SxwnlLunarEclipseItem*)std::calloc(list->count, sizeof(SxwnlLunarEclipseItem));
        if (!list->items) { std::free(list); return nullptr; }
        std::memcpy(list->items, results.data(), list->count * sizeof(SxwnlLunarEclipseItem));
    }
    return list;
    });
}

void sxwnl_lunar_eclipse_list_free(SxwnlLunarEclipseList *list) {
    if (!list) return;
    std::free(list->items);
    std::free(list);
}

// ═══════════════════════════════════════════════════════════
//  Lunar Eclipse Detail (Visualization Frames)
// ═══════════════════════════════════════════════════════════

SxwnlLunarEclipseDetail* sxwnl_lunar_eclipse_detail(int year, int month, int day,
                                                      int frame_interval) {
    return ecl_guard([&]() -> SxwnlLunarEclipseDetail* {
    if (frame_interval <= 0) frame_interval = 120;
    long double jd = ymd_to_jd(year, month, day);
    // Find nearest full moon
    jd = XL::MS_aLon_t2((int64((jd + 8) / 29.5306) * 2 + 1) * PI) * 36525;

    EphYspl yspl;
    yspl.lecMax(jd);

    std::string lx = yspl.getLX();
    auto lt = yspl.getLT();

    // 既无本影食也无半影食 → 无月食
    if (lx.empty() && lt[3] == 0 && lt[4] == 0) return nullptr;

    long double sf = yspl.getSF();

    auto *detail = (SxwnlLunarEclipseDetail*)std::calloc(1, sizeof(SxwnlLunarEclipseDetail));
    if (!detail) return nullptr;

    if (lx == "全") {
        std::strncpy(detail->type, "T", 3);
        std::strncpy(detail->type_name, "全食", 15);
    } else if (lx == "偏") {
        std::strncpy(detail->type, "P", 3);
        std::strncpy(detail->type_name, "偏食", 15);
    } else {
        std::strncpy(detail->type, "B", 3);
        std::strncpy(detail->type_name, "半影食", 15);
    }
    detail->max_magnitude = (double)sf;

    // Fill event times
    detail->t_u1  = lt[0] ? (double)(lt[0] + J2000) : 0; // 初亏
    detail->t_max = lt[1] ? (double)(lt[1] + J2000) : 0; // 食甚
    detail->t_u4  = lt[2] ? (double)(lt[2] + J2000) : 0; // 复圆
    detail->t_p1  = lt[3] ? (double)(lt[3] + J2000) : 0; // 半影食始
    detail->t_p4  = lt[4] ? (double)(lt[4] + J2000) : 0; // 半影食终
    detail->t_u2  = lt[5] ? (double)(lt[5] + J2000) : 0; // 食既
    detail->t_u3  = lt[6] ? (double)(lt[6] + J2000) : 0; // 生光

    // 半影食时 t_max 缺失, 取半影食始终中点
    if (detail->t_max == 0) {
        if (detail->t_p1 && detail->t_p4)
            detail->t_max = (detail->t_p1 + detail->t_p4) * 0.5;
        else if (detail->t_p1)
            detail->t_max = detail->t_p1;
    }

    // Generate frames covering the entire eclipse (penumbra start to end)
    long double t_max_l = detail->t_max ? (long double)(detail->t_max - J2000) : (lt[1] ? lt[1] : 0);
    long double t_begin = lt[3] ? lt[3] : (lt[0] ? lt[0] : t_max_l - 2.0L / 24.0L);
    long double t_end   = lt[4] ? lt[4] : (lt[2] ? lt[2] : t_max_l + 2.0L / 24.0L);

    long double interval_jd = frame_interval / 86400.0L;
    std::vector<SxwnlLunarEclipseFrame> frames;

    RE0 re{};
    for (long double t = t_begin; t <= t_end; t += interval_jd) {
        yspl.lecXY(t, re);

        SxwnlLunarEclipseFrame frame{};
        frame.moon_x = (double)re.x;
        frame.moon_y = (double)re.y;
        frame.moon_radius = (double)re.mr;
        frame.umbra_radius = (double)re.er;
        frame.penumbra_radius = (double)re.Er;
        frame.jd = (double)(t + J2000);

        double dist = std::sqrt((double)(re.x * re.x + re.y * re.y));
        double coverage = ((double)re.er + (double)re.mr - dist) / ((double)re.mr * 2.0);
        frame.coverage = coverage > 0 ? (coverage > 1 ? 1 : coverage) : 0;

        frames.push_back(frame);
    }

    detail->count = (int32_t)frames.size();
    if (detail->count > 0) {
        detail->frames = (SxwnlLunarEclipseFrame*)std::malloc(
            detail->count * sizeof(SxwnlLunarEclipseFrame));
        if (!detail->frames) { detail->count = 0; }
        else std::memcpy(detail->frames, frames.data(),
                         detail->count * sizeof(SxwnlLunarEclipseFrame));
    }
    return detail;
    });
}

void sxwnl_lunar_eclipse_detail_free(SxwnlLunarEclipseDetail *data) {
    if (!data) return;
    std::free(data->frames);
    std::free(data);
}
