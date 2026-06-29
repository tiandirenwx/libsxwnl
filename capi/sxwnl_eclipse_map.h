#ifndef SXWNL_ECLIPSE_MAP_H
#define SXWNL_ECLIPSE_MAP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdbool.h>

/* ═══════════════════════════════════════════════════════════
 *  Eclipse Map Visualization C API
 *
 *  Provides structured data for rendering eclipse paths on maps
 *  and eclipse animations. Supports both solar and lunar eclipses.
 *
 *  Memory contract:
 *  - sxwnl_eclipse_xxx_free() frees any returned struct pointer
 *  - Arrays inside structs are owned by the parent; freed together
 * ═══════════════════════════════════════════════════════════ */

// ─── Common types ───

typedef struct {
    double longitude;   // degrees, east positive
    double latitude;    // degrees, north positive
    double jd;          // Julian Day (TD)
} SxwnlGeoPoint;

// ─── Solar Eclipse: Search Result (structured) ───

typedef struct {
    double jd;          // Julian Day of maximum eclipse (TD)
    int32_t year;
    int32_t month;
    int32_t day;
    int32_t hour;
    int32_t minute;
    char type[4];       // "T","A","P","H","T0","A0","T1","A1","H2","H3"
    char type_name[16]; // "全食","环食","偏食" etc (UTF-8)
    double gamma;       // shadow axis distance to earth center
    double magnitude;   // eclipse magnitude (食分)
    double center_lon;  // center point longitude (100 = no center)
    double center_lat;  // center point latitude
    double path_width;  // shadow path width in km (0 for partial)
    double duration;    // max duration in seconds (0 for partial)
    int32_t saros;      // Saros series number
    int32_t saros_member; // member within series
    char season[16];    // 食季: "春季食季","秋季食季" etc
} SxwnlSolarEclipseItem;

typedef struct {
    SxwnlSolarEclipseItem *items;
    int32_t count;
} SxwnlSolarEclipseList;

// ─── Solar Eclipse: Map Path Data ───

typedef struct {
    SxwnlGeoPoint *center_line;       // center line path points
    int32_t center_count;

    SxwnlGeoPoint *umbra_north;       // umbral shadow north limit
    int32_t umbra_north_count;

    SxwnlGeoPoint *umbra_south;       // umbral shadow south limit
    int32_t umbra_south_count;

    SxwnlGeoPoint *penumbra_north;    // penumbral shadow north limit
    int32_t penumbra_north_count;

    SxwnlGeoPoint *penumbra_south;    // penumbral shadow south limit
    int32_t penumbra_south_count;

    // Eclipse overview
    char type[4];
    char type_name[16];
    double gamma;
    double magnitude;
    double path_width_km;
    double duration_seconds;
    double max_eclipse_jd;            // JD of greatest eclipse

    // Key event points (lon, lat, jd)
    SxwnlGeoPoint partial_start;      // first external contact (偏食始)
    SxwnlGeoPoint central_start;      // first central contact (中心始)
    SxwnlGeoPoint greatest_eclipse;   // greatest eclipse point
    SxwnlGeoPoint central_end;        // last central contact (中心终)
    SxwnlGeoPoint partial_end;        // last external contact (偏食终)
} SxwnlSolarEclipsePath;

// ─── Solar Eclipse: Local Animation Frame ───

typedef struct {
    double sun_x;       // sun center (reference, always 0)
    double sun_y;
    double sun_radius;  // angular radius in radians
    double moon_x;      // moon center offset from sun (radians)
    double moon_y;
    double moon_radius; // angular radius in radians
    double jd;          // frame time (TD)
    double magnitude;   // instantaneous magnitude at this frame
} SxwnlEclipseFrame;

typedef struct {
    SxwnlEclipseFrame *frames;
    int32_t count;

    // Event times (JD TD; 0 = not applicable)
    double t_c1;        // first contact (初亏)
    double t_c2;        // second contact (食既/环食始; 0 if partial only)
    double t_max;       // maximum eclipse (食甚)
    double t_c3;        // third contact (生光/环食终; 0 if partial only)
    double t_c4;        // fourth contact (复圆)
    double t_sunrise;   // local sunrise
    double t_sunset;    // local sunset

    double max_magnitude;   // max magnitude at this location
    double moon_sun_ratio;  // moon/sun diameter ratio at max
    char type[8];           // "偏","全","环"
} SxwnlLocalSolarEclipse;

// ─── Lunar Eclipse: Search Result (structured) ───

typedef struct {
    double jd;          // Julian Day of maximum eclipse (TD)
    int32_t year;
    int32_t month;
    int32_t day;
    int32_t hour;
    int32_t minute;
    char type[4];       // "N"(none),"B"(penumbral),"P"(partial),"T"(total)
    char type_name[16]; // "全食","偏食","半影食"
    double magnitude;   // eclipse magnitude (食分)
    int32_t saros;      // Saros series number
    int32_t saros_member; // member within series
    char season[16];    // 食季
} SxwnlLunarEclipseItem;

typedef struct {
    SxwnlLunarEclipseItem *items;
    int32_t count;
} SxwnlLunarEclipseList;

// ─── Lunar Eclipse: Visualization Data ───

typedef struct {
    double moon_x;          // moon center x offset (radians, earth shadow at origin)
    double moon_y;          // moon center y offset
    double moon_radius;     // moon angular radius (radians)
    double umbra_radius;    // earth umbral shadow radius (radians)
    double penumbra_radius; // earth penumbral shadow radius (radians)
    double jd;              // frame time (TD)
    double coverage;        // fraction of moon diameter covered by umbra
} SxwnlLunarEclipseFrame;

typedef struct {
    SxwnlLunarEclipseFrame *frames;
    int32_t count;

    // Event times (JD TD; 0 = not applicable)
    double t_p1;        // penumbral eclipse begins (半影食始)
    double t_u1;        // partial eclipse begins (初亏)
    double t_u2;        // total eclipse begins (食既; 0 if partial only)
    double t_max;       // maximum eclipse (食甚)
    double t_u3;        // total eclipse ends (生光; 0 if partial only)
    double t_u4;        // partial eclipse ends (复圆)
    double t_p4;        // penumbral eclipse ends (半影食终)

    double max_magnitude;
    char type[4];       // "P","T"
    char type_name[16]; // "偏食","全食"
} SxwnlLunarEclipseDetail;

// ═══════════════════════════════════════════════════════════
//  API Functions
// ═══════════════════════════════════════════════════════════

// --- Solar Eclipse ---

// Search for solar eclipses starting from (year, month), returning up to `count` results
SxwnlSolarEclipseList* sxwnl_solar_eclipse_search(int year, int month, int count);
void sxwnl_solar_eclipse_list_free(SxwnlSolarEclipseList *list);

// Get the map path for a solar eclipse nearest to the given date
// The date should be close to a new moon (within a few days)
SxwnlSolarEclipsePath* sxwnl_solar_eclipse_path(int year, int month, int day);
void sxwnl_solar_eclipse_path_free(SxwnlSolarEclipsePath *path);

// Get animation frames for a local solar eclipse observation
// Returns null if no eclipse visible at (longitude, latitude) around the given date
// frame_interval: seconds between frames (suggest 60-120 for smooth animation)
SxwnlLocalSolarEclipse* sxwnl_solar_eclipse_local(
    int year, int month, int day,
    double longitude, double latitude, int frame_interval);
void sxwnl_solar_eclipse_local_free(SxwnlLocalSolarEclipse *data);

// --- Lunar Eclipse ---

// Search for lunar eclipses starting from (year, month), returning up to `count` results
SxwnlLunarEclipseList* sxwnl_lunar_eclipse_search(int year, int month, int count);
void sxwnl_lunar_eclipse_list_free(SxwnlLunarEclipseList *list);

// Get visualization frames for a lunar eclipse nearest to the given date
// frame_interval: seconds between frames (suggest 120-300)
SxwnlLunarEclipseDetail* sxwnl_lunar_eclipse_detail(int year, int month, int day,
                                                      int frame_interval);
void sxwnl_lunar_eclipse_detail_free(SxwnlLunarEclipseDetail *data);

#ifdef __cplusplus
}
#endif

#endif // SXWNL_ECLIPSE_MAP_H
