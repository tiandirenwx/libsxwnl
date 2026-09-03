#ifndef SXWNL_CAPI_H
#define SXWNL_CAPI_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdbool.h>

/* ═══════════════════════════════════════════════════════════
 *  Cross-platform C API for libsxwnl
 *
 *  Memory contract:
 *  - sxwnl_xxx_create()  → caller owns, must call sxwnl_xxx_free()
 *  - sxwnl_xxx_get_yyy() → returns pointers valid until parent is freed
 *  - sxwnl_string_free() → free any returned char*
 * ═══════════════════════════════════════════════════════════ */

// ─── Day Info (value type, no manual free needed) ───

typedef struct {
    int32_t solar_year;
    int32_t solar_month;
    int32_t solar_day;

    int32_t lunar_year;
    int32_t lunar_month;
    int32_t lunar_day;
    bool    is_leap_month;

    int32_t week_day;           // 0=Sunday

    int32_t year_gan;           // 天干 index 0-9
    int32_t year_zhi;           // 地支 index 0-11
    int32_t month_gan;
    int32_t month_zhi;
    int32_t day_gan;
    int32_t day_zhi;

    int32_t jie_qi;             // -1 if none, 0-23 (天文口径: qi_accurate 精确交气时刻所在日)
    int32_t yue_xiang;          // -1 if none, 0-3
    int32_t constellation;      // 星座 0-11
    int32_t jian12;             // 建除十二 0-11

    char year_gz[8];            // "甲子" UTF-8, null-terminated
    char month_gz[8];
    char day_gz[8];
    char lunar_month_name[16];  // "正月", "闰二月"
    char lunar_day_name[12];    // "初一", "十五"
    char jie_qi_name[12];       // "" if none (天文口径, 与 jie_qi 对应)
    char jie_qi_time[32];       // "" if none (天文精确交气时刻 HH:MM:SS)
    char sheng_xiao[8];         // "鼠"
    char constellation_name[12];
    char week_name[16];         // "星期一"
    char yue_xiang_name[12];    // "朔","望" or ""
    char jian12_name[8];        // "建","除"...

    // ── 节日 (追加字段, 向后兼容) ────────────────
    char holiday[128];          // A 类: 法定/重要节日, 空串表示无
    char major[256];            // B 类: 中等节日
    char minor[256];            // C 类: 一般节日
    char misc[64];              // 杂节: 数九/三伏/入梅/出梅
    bool is_off_day;            // 是否为放假日(法定 + 周末)

    // ── 纪年扩展 ────────────────────────────────
    int32_t lunar_jun_year;     // 农历纪年(春节为界)
    int32_t lunar_lichun_year;  // 农历纪年(立春为界)
    int32_t huangdi_year;       // 黄帝纪年

    // ── 回历(Hijri) ─────────────────────────────
    int32_t moslem_year;        // 回历年
    int32_t moslem_month;       // 回历月 1-12
    int32_t moslem_day;         // 回历日 1-30

    // ── 纳音(由年/月/日干支推得) ─────────────────
    char year_nayin[16];        // 年柱纳音, 例 "海中金"
    char month_nayin[16];       // 月柱纳音
    char day_nayin[16];         // 日柱纳音

    // ── 月相时刻 ──────────────────────────────
    char yue_xiang_time[32];    // 该日月相极值时刻 "HH:MM:SS", 空串表示无

    // ── 儒略日(对应当日 12:00, 与 sxwnl/lunar.js Bd0+J2000 一致) ─
    double julian_day;

    // ── 历谱口径节气 (追加字段, 向后兼容) ──────────────
    //  与上面 jie_qi/jie_qi_name(天文口径)的区别:
    //  历谱口径采用整日表 + QB 定气修正(对应 sxwnl 网页版 ob.Ljq),
    //  对古代(1645 年以前)节气所在公历日可能与天文口径相差 1 天,
    //  与权威 sxwnl 网页版日历"节气日"一致; 1645 年以后两者相同。
    //  历谱口径是"整日"概念, 无精确时刻(精确时刻仍取 jie_qi_time)。
    //  端上日历格子的节气标签建议用 lipu_jie_qi_name, 精确交节时刻用 jie_qi_time。
    int32_t lipu_jie_qi;        // -1 if none, 0-23 (历谱口径)
    char    lipu_jie_qi_name[12]; // "" if none (历谱口径)
} SxwnlDayInfo;

// ─── Calendar functions ───

int  sxwnl_get_day_info(int year, int month, int day, SxwnlDayInfo *out);
int  sxwnl_get_month_days(int year, int month, SxwnlDayInfo *out, int max_count);
int  sxwnl_lunar_to_solar(int year, int month, int day, bool is_leap,
                           int *out_year, int *out_month, int *out_day);
int  sxwnl_solar_to_lunar(int year, int month, int day, SxwnlDayInfo *out);
int  sxwnl_get_year_leap_month(int year);

// 某农历年的逐月信息(正月/闰四月/后九月/壹月...), 供输入选择
typedef struct {
    int  month;       // 农历月序 1-12 (1=正月), 用于农历->公历
    int  is_leap;     // 闰月
    int  is_spec;     // 特殊月(SYmc)
    char name[16];    // 显示名: 正月/闰四月/后九月/壹月...
} SxwnlLunarMonth;

int  sxwnl_get_lunar_months(int year, SxwnlLunarMonth *out, int max_count);

// 某农历月的天数(29 或 30)
int  sxwnl_get_lunar_month_days(int year, int month, bool is_leap, bool is_spec);

// 农历日名(初一/初二../三十). day 为 1..30; 越界写入数字串兜底。
// 写入 out(以 '\0' 结尾), 返回写入的字节数(不含结尾符); 失败返回 0。
// 底层复用 Rmc[] 表, 供各端(月历/八字选择器)统一取名, 避免各平台重复硬编码。
int  sxwnl_get_lunar_day_name(int day, char *out, int cap);

// 某公历年月中"真实存在"的日号列表(写入 out, 升序), 返回天数。
// 依据本库历法(儒略<->格里), 自动跳过 1582-10-05..14 改历缺失的 10 天,
// 并自然给出各月正确长度(平/闰年二月、30/31 天等)。out 建议容量 >= 31。
int  sxwnl_get_solar_month_valid_days(int year, int month, int *out, int max_count);

// ─── 回历 ────────────────────────────────────────────────
// 公历 → 回历
int  sxwnl_solar_to_moslem(int year, int month, int day,
                           int *out_h_year, int *out_h_month, int *out_h_day);
// 回历 → 公历
int  sxwnl_moslem_to_solar(int h_year, int h_month, int h_day,
                           int *out_year, int *out_month, int *out_day);
// 由回历日期获取完整日历信息
int  sxwnl_get_day_info_by_moslem(int h_year, int h_month, int h_day,
                                  SxwnlDayInfo *out);

// ─── 综合月历 (公/农/回 三历 + 干支 + 节气 + 月相 + 节日) ──
typedef struct {
    int32_t year;                // 公历年
    int32_t month;               // 公历月 1-12
    int32_t first_week_day;      // 月首星期 0-6 (0=日)
    int32_t day_count;           // 本月公历天数
    int32_t total_weeks;         // 本月覆盖的周数
    int32_t first_julian_day;    // 月首儒略日 (J2000 起算)

    int32_t year_gan;            // 干支纪年(春节为界): 天干 0-9
    int32_t year_zhi;            // 干支纪年(春节为界): 地支 0-11
    char    year_gz[8];          // "甲辰" 等
    char    sheng_xiao[8];       // "龙"
    char    nianhao[256];        // 年号(可为多个朝代拼接)
} SxwnlCalendarMonthHeader;

// 一次性算出整月数据.
//   year, month: 公历年月
//   out_header (可空): 整月级别元数据
//   out_days   (可空): 每日详情数组. 容量需 >= 31.
//   max_days   : out_days 容量
// 返回写入的天数(本月实际天数). 错误返回 -1.
int  sxwnl_get_calendar_month(int year, int month,
                              SxwnlCalendarMonthHeader *out_header,
                              SxwnlDayInfo *out_days, int max_days);

// ─── JieQi ───

typedef struct {
    int32_t idx;                // 节气索引 0-23
    int32_t solar_month;
    int32_t solar_day;
    char    name[12];
    char    time[32];
} SxwnlJieQiItem;

int  sxwnl_get_jieqi_list(int year, SxwnlJieQiItem *out, int max_count);

// ─── 年历 (按农历月聚合, 参考 sxwnl/lunar.js 中 nianLi2HTML) ───
//
// 给定公历年, 返回该年覆盖到的农历月列表(通常 12-14 个);
// 每个月含: 朔日干支+公历日期、月名/大小/闰、本月内出现的节气详情。
// 适合用来生成"年历"页面.
// ─────────────────────────────────────────────────────────

typedef struct {
    int32_t idx;                // 节气索引 0..23 (0=冬至)
    char    name[12];           // 节气名称, 如 "小寒"
    char    gz[8];              // 节气日的干支(以日柱算)
    int32_t solar_month;        // 历谱口径公历月 1-12 (整日气表所在日, 古代平气)
    int32_t solar_day;          // 历谱口径公历日
    char    time[32];           // 节气精确交气时刻 "HH:MM:SS" (天文定气)
    int32_t day_offset;         // 距本农历月首(初一)的天数 0..29
    char    day_name[12];       // 本月内日名称, 如 "初一"/"十五"
    // 精确交气(天文定气)所在公历日期. 1645 年前平气历谱日可能与其差 1 天;
    // 现代(定气历)则与 solar_month/solar_day 相同. 用于完整还原
    // sxwnl nianLiHTML 的 "历谱MM-DD(精确日 时:分:秒)" 展示.
    int32_t acc_month;          // 精确交气公历月 1-12
    int32_t acc_day;            // 精确交气公历日
} SxwnlYearCalJieQi;

typedef struct {
    int32_t month_idx;          // 农历月名索引 0..11 (0=正)
    char    month_name[16];     // "正月"/"闰二月"/"后九月"等
    int32_t is_leap;            // 是否闰月
    int32_t is_spec;            // 是否特殊月(SYmc)
    int32_t day_count;          // 本农历月天数 29/30
    int32_t solar_year;         // 朔日所属公历年
    int32_t solar_month;        // 朔日公历月
    int32_t solar_day;          // 朔日公历日
    char    shuo_gz[8];         // 朔日干支
    int32_t jq_count;           // 本月内节气个数
    SxwnlYearCalJieQi jieqi[4]; // 通常 1-2 个, 预留 4
} SxwnlYearCalMonth;

// 返回写入的农历月数(本年内合朔月). 错误返回 -1.
// out 容量建议 >= 14
int  sxwnl_get_year_calendar(int year, SxwnlYearCalMonth *out, int max_count);

// ─── Bazi ───

typedef struct {
    char shi_shen[8];
    char tian_gan[8];
    char di_zhi[8];
    char cang_gan1[8];
    char cang_gan2[8];
    char cang_gan3[8];
} SxwnlPillar;

typedef struct {
    char    gan_zhi[8];
    int32_t start_year;
    int32_t end_year;
} SxwnlDaYun;

typedef void* SxwnlBazi;

typedef struct {
    char     name[64];
    bool     gender;        // false=男, true=女
    bool     is_ast;        // 真太阳时
    double   longitude;
    double   latitude;
    int32_t  lifa;          // LiFaType enum (0-12)
    bool     is_lunar;
    bool     is_leap;       // 闰月
    bool     is_spec;       // 特殊月(古代历法, SYmc)
    int32_t  year, month, day;
    int32_t  hour, minute;
    double   second;
} SxwnlBaziParams;

SxwnlBazi sxwnl_bazi_create(const SxwnlBaziParams *params);
void      sxwnl_bazi_free(SxwnlBazi bazi);

const char*    sxwnl_bazi_get_user_name(SxwnlBazi bazi);
const char*    sxwnl_bazi_get_gender(SxwnlBazi bazi);
const char*    sxwnl_bazi_get_solar_birth(SxwnlBazi bazi);
const char*    sxwnl_bazi_get_lunar_birth(SxwnlBazi bazi);
const char*    sxwnl_bazi_get_date_of_birth(SxwnlBazi bazi);
const char*    sxwnl_bazi_get_sheng_xiao(SxwnlBazi bazi);
const char*    sxwnl_bazi_get_age(SxwnlBazi bazi);
const char*    sxwnl_bazi_get_lifa(SxwnlBazi bazi);
const char*    sxwnl_bazi_get_ding_qi_type(SxwnlBazi bazi);
const char*    sxwnl_bazi_get_ast(SxwnlBazi bazi);
const char*    sxwnl_bazi_get_jie_qi(SxwnlBazi bazi);
const char*    sxwnl_bazi_get_qi_yun(SxwnlBazi bazi);
const char*    sxwnl_bazi_get_jiao_yun(SxwnlBazi bazi);

int sxwnl_bazi_get_si_zhu(SxwnlBazi bazi, SxwnlPillar out[4]);
int sxwnl_bazi_get_da_yun_count(SxwnlBazi bazi);
int sxwnl_bazi_get_da_yun(SxwnlBazi bazi, SxwnlDaYun *out, int max_count);
int sxwnl_bazi_get_fleeting_years(SxwnlBazi bazi, char (*out)[8], int max_count);

// 时辰干支 (standalone, no bazi context needed)
void sxwnl_get_hour_gz(int day_gan, int hour, char out[8]);

// ─── Bazi 完整排盘列信息(命理派生数据) ───

// 一柱(年/月/日/时/大运/流年)的完整盘面信息
typedef struct {
    char gan[8];                 // 天干
    char zhi[8];                 // 地支
    char gan_shi_shen[12];       // 主星(天干十神); 日柱填性别"男"/"女"
    char nayin[16];              // 纳音
    char xing_yun[8];            // 星运: 日主在该地支的十二长生
    char zi_zuo[8];              // 自坐: 本柱天干在本柱地支的十二长生
    char kong_wang[12];          // 空亡(两地支)
    char cang_gan[3][8];         // 藏干天干
    char cang_gan_shi_shen[3][8];// 藏干十神(简称)
    int  cang_gan_count;
    char shen_sha[12][20];       // 神煞名
    int  shen_sha_count;
    int  start_year;             // 大运/流年: 起始公历年; 四柱填0
} SxwnlBaziColumn;

// 四柱(年月日时)完整列信息, 返回填充的列数(4)
int sxwnl_bazi_get_columns(SxwnlBazi bazi, SxwnlBaziColumn out[4]);
// 当前此刻所行大运列(0成功 -1未起运/失败)
int sxwnl_bazi_get_current_da_yun(SxwnlBazi bazi, SxwnlBaziColumn *out);
// 当前流年列
int sxwnl_bazi_get_current_liu_nian(SxwnlBazi bazi, SxwnlBaziColumn *out);
// 大运列表(完整列)
int sxwnl_bazi_get_da_yun_columns(SxwnlBazi bazi, SxwnlBaziColumn *out, int max_count);

// 司令天干名(人元司令分野)
const char* sxwnl_bazi_get_si_ling(SxwnlBazi bazi);
// 五行统计(木火土金水); include_cang_gan=true 含藏气
void sxwnl_bazi_get_wuxing_count(SxwnlBazi bazi, int out[5], bool include_cang_gan);
// 五行旺衰(木火土金水) → "旺""相""休""囚""死"
void sxwnl_bazi_get_wuxing_status(SxwnlBazi bazi, char out[5][8]);

// 流年/小运项
typedef struct {
    int  year;                   // 公历年
    int  age;                    // 虚岁
    char gan_zhi[8];             // 流年干支
    char gan_shi_shen[8];        // 流年天干十神(简称)
    char zhi_shi_shen[8];        // 流年地支本气十神(简称)
    char xiao_yun[8];            // 小运干支
    char xiao_yun_shi_shen[8];   // 小运天干十神(简称)
} SxwnlLiuNianItem;

// 指定起始公历年起 count 年的流年+小运
int sxwnl_bazi_get_liu_nian(SxwnlBazi bazi, int start_year,
                            SxwnlLiuNianItem *out, int max_count);

// ─── 四柱反推(由四柱干支反查公历日期) ───

// 一个反推结果
typedef struct {
    int  year;          // 公历年(可为负, 表示天文年, year<=0 即公元前)
    int  month;         // 公历月 1-12
    int  day;           // 公历日
    int  hour;          // 时辰起始小时 0-23; 时辰未知时为 -1
    char ganzhi[4][8];  // 年月日时四柱干支
} SxwnlReverseItem;

// 四柱反推: sz[8]=年干年支月干月支日干日支时干时支(索引);
// 若时柱未知, 传 sz[6]=sz[7]=-1; 在 [start_year,end_year] 公历年内搜索。
// 返回匹配数量(写入 out, 最多 max_count)。
int sxwnl_bazi_reverse(const int sz[8], int start_year, int end_year,
                       SxwnlReverseItem *out, int max_count);

// ─── 节日 (独立查询, 与 sxwnl_get_day_info 共享数据) ───

typedef struct {
    char holiday[128];   // A 类
    char major[256];     // B 类
    char minor[256];     // C 类
    char misc[64];       // 杂节
    bool is_off_day;     // 放假
} SxwnlFestival;

// 公历(year-month-day) 一日的节日信息
int sxwnl_get_festivals(int year, int month, int day, SxwnlFestival *out);

// ─── 纪年 / 时间字符串工具 ───

// "公元前 221"/"前 221"/"-221"/"2024年" → 天文学纪年(整数)
// 失败返回 INT32_MIN
int32_t sxwnl_year_str_to_astro(const char *s);

// 天文学纪年 → "公元X年"/"公元前X年"; full_style=0 时返回简写 "X"/"前X"
// 写入 out, 返回 0 成功; out_size 建议 >= 32
int sxwnl_astro_year_to_str(int32_t aYear, bool full_style, char *out, int out_size);

// "HH"/"HH:MM"/"HH:MM:SS[.fff]"/"12.5" → 浮点小时; 失败返回 NaN
double sxwnl_time_str_to_hour(const char *s);

// ─── 月球与行星天象 ───

// 月相: 0=朔 1=上弦 2=望 3=下弦
typedef struct {
    int32_t phase;        // 0..3
    char    name[8];      // "朔"/"上弦"/"望"/"下弦"
    char    time[32];     // ISO-like 时间字符串
    double  jd;           // 北京时(TD->UTC 已转换)的儒略日
} SxwnlMoonPhase;

// 列出 year 整年内的月相, 写入 out
int sxwnl_get_moon_phases(int year, SxwnlMoonPhase *out, int max_count);

// 月地距极值
typedef struct {
    int32_t kind;         // 0=近地点 1=远地点
    char    name[8];      // "近地"/"远地"
    char    time[32];
    double  jd;
    double  distance_km;  // 地心距(km)
} SxwnlMoonApsis;

int sxwnl_get_moon_apsides(int year, SxwnlMoonApsis *out, int max_count);

// 月交点
typedef struct {
    int32_t kind;         // 0=升 1=降
    char    name[8];      // "升交"/"降交"
    char    time[32];
    double  jd;
} SxwnlMoonNode;

int sxwnl_get_moon_nodes(int year, SxwnlMoonNode *out, int max_count);

// 地日距极值
typedef struct {
    int32_t kind;         // 0=近日 1=远日
    char    name[8];
    char    time[32];
    double  jd;
    double  distance_au;
} SxwnlEarthApsis;

int sxwnl_get_earth_apsides(int year, SxwnlEarthApsis *out, int max_count);

// ─── 行星天象 ───────────────────────────────────────────────
//
//  planet 编号:
//    1=水星  2=金星  3=火星  4=木星  5=土星  6=天王  7=海王
//
//  事件类型:
//    0=合(上合/合日)         适用全部
//    1=冲(外行星)/下合(内行星) 适用全部
//    2=东大距                 仅 1,2
//    3=西大距                 仅 1,2
//    4=顺留                   适用 1..7
//    5=逆留                   适用 1..7
//    6=合月                   适用 1..7
// ──────────────────────────────────────────────────────────

typedef struct {
    int32_t planet;       // 1..7
    int32_t event;        // 上述类型
    char    planet_name[8];
    char    event_name[12];
    char    time[32];
    double  jd;           // 北京时儒略日
    double  value;        // 事件相关值:
                          //   合/冲/合月: 黄经差(度), 接近 0
                          //   大距:       角距(度)
                          //   留:         0
} SxwnlPlanetEvent;

// 列出 year 内全部行星天象(自动按时间排序)
int sxwnl_get_planet_events(int year, SxwnlPlanetEvent *out, int max_count);

// 行星视位置详情(返回 sxwnl 内部格式化字符串; caller 需 sxwnl_string_free)
//   planet: 1..7
//   bj_jd:  北京时儒略日(可用 sxwnl_time_str_to_hour 或 Time→JD 自行构造)
char* sxwnl_calc_planet_position(int planet, double bj_jd,
                                 double longitude, double latitude);

// ─── 地图投影 ───────────────────────────────────────────────
//
//  9 种投影类型(同 sxwnl/eph0.js touY):
//    0 平面正投            1 斜轴等距方位          2 斜轴等积方位
//    3 斜轴等角方位        4 摩尔威特(Mollweide)
//    5 正轴等距圆柱        6 正轴等角圆柱(Mercator)
//    7 多圆锥              8 中国灯笼
//
//  使用流程:
//    handle = sxwnl_projection_create(...);
//    sxwnl_projection_apply(handle, lonlatArray, out);
//    sxwnl_projection_free(handle);
// ──────────────────────────────────────────────────────────

typedef void* SxwnlProjection;

// 创建投影对象
//   type    : 0..8
//   J0_rad  : 基准经度(弧度)
//   W0_rad  : 基准纬度(弧度)
//   win_cx/cy/rx/ry: 可见矩形窗口(投影后空间, 默认 0,0,1,1)
SxwnlProjection sxwnl_projection_create(int type, double J0_rad, double W0_rad,
                                        double win_cx, double win_cy,
                                        double win_rx, double win_ry);
void sxwnl_projection_free(SxwnlProjection proj);

// 单点投影 (J, W 弧度) -> (x, y, z); z<0 表示不可见
int sxwnl_projection_point(SxwnlProjection proj,
                           double J_rad, double W_rad,
                           double *out_x, double *out_y, double *out_z);

// 折线投影
//   input  : [J,W,J,W,...]  1e7 表示 moveto
//   in_count: input 元素数
//   out    : 调用方分配的缓冲区
//   out_max: out 容量(元素数)
//   返回写入元素数(2*点数 + moveto 标记); 不够则返回所需大小的负值
int sxwnl_projection_polyline(SxwnlProjection proj,
                              const double *input, int in_count,
                              double *out, int out_max);

// ─── 世界地图数据 ─────────────────────────────────────────────

// 返回内置 ditu0 (小图) 的元素数; 不超过 out_max 则写入 out
int sxwnl_world_map_get_ditu0(double *out, int out_max);

// 注入大图字符串数据 (ditu1: idx=1, ditu2: idx=2)
//   返回 0 成功
int sxwnl_world_map_set_data(int idx, const char *raw);

// 拷贝已注入的数据
int sxwnl_world_map_get_data(int idx, double *out, int out_max);

// ─── 88 星座 ────────────────────────────────────────────────
typedef struct {
    char    name_abbr[24];   // "仙女座And"
    double  area_sq;         // 面积(平方度)
    char    ra_str[16];      // "0 48.46"
    char    dec_str[16];     // "37 25.91"
    char    quad_family[16]; // "NQ1 英仙"
    char    name_en[32];     // "Andromeda"
} SxwnlConstellation;

// 取出 88 个星座(out 至少 88), 返回写入条数
int sxwnl_get_constellations(SxwnlConstellation *out, int max_count);

// ─── 恒星 ──────────────────────────────────────────────────
typedef struct {
    double ra_rad;       // 赤经 (J2000, 弧度)
    double dec_rad;      // 赤纬 (J2000, 弧度)
    double pm_ra;        // 赤经自行(弧度/世纪)
    double pm_dec;       // 赤纬自行(弧度/世纪)
    double parallax;     // 视差(弧度)
    double mag;          // 星等
    char   name[32];     // 星名
    char   info[32];     // 星座/编号等
} SxwnlStar;

// 注册自定义恒星库 (上游 HXK 字符串格式; # 分隔, * 主星)
int sxwnl_register_star_library(const char *key, const char *raw);

// 取出指定库的恒星 (include_all=1 返回非 * 项)
int sxwnl_get_star_library(const char *key, bool include_all,
                           SxwnlStar *out, int max_count);

// 按关键字检索全部库 (星名/信息/星座缩写)
int sxwnl_search_stars(const char *key, SxwnlStar *out, int max_count);

// hxCalc 多颗恒星位置
//   mode: 0=视位置 1=站心(方位/高度) 2=平位置
//   bj_jd: 北京时儒略日 (内部转 TT 儒略世纪)
//   longitude/latitude (度): 站心模式有效
typedef struct {
    char   name[64];    // 星名 + 附加信息
    double mag;
    double a_rad;       // 视赤经 / 方位角 / 平赤经
    double b_rad;       // 视赤纬 / 高度角 / 平赤纬
} SxwnlStarResult;

int sxwnl_star_hx_calc(const SxwnlStar *in, int in_count,
                       double bj_jd, double nutation_q_days,
                       int mode, double longitude, double latitude,
                       SxwnlStarResult *out, int max_count);

// ─── 老黄历 (Almanac) ──────────────────────────────────────
//
//  设计契约: 老黄历模块"不重算干支/农历/节气" - 所有输入由 libsxwnl
//  内部算出再喂给 almanac 查表器, 保证全项目唯一真源, 避免节气交界
//  日上日历与黄历显示不一致的 bug.
//
//  数据移植自 lunar-javascript (MIT).
// ───────────────────────────────────────────────────────────

// 择日典籍语录 (董公择日要诀 / 玉匣记 / 通胜经 ...)
//   同一日不同典籍各占一条; quote_count 给出实际条数.
typedef struct {
    char source[24];           // 来源 "董公择日要诀"
    char title[40];            // 标题 "正月·开子日"
    char luck[8];              // 整体基调 "吉"/"凶"/"平"/"混"/""
    char body[1024];           // 正文 (单条最长约 ~700 字节)
} SxwnlAlmanacQuote;

#define SXWNL_ALMANAC_QUOTE_MAX 4 // 当前仅董公诀, 预留 4 槽位

// 神煞 (天德/月厌/白虎入中宫/三合 ...)
//   name 留 32 字节 — 最长 "人中三奇(壬癸辛)" 23 字节, 给 32 留余量
typedef struct {
    char name[32];             // "天德"/"月厌大祸"/"人中三奇(壬癸辛)"
    bool is_lucky;             // 吉神=true / 凶神=false
    int32_t weight;            // 权重 1-3 (1一般, 2中, 3大煞)
} SxwnlShenSha;

#define SXWNL_ALMANAC_SHENSHA_MAX 24  // 一日命中神煞通常 5-15 条

// 吉时
typedef struct {
    char name[12];             // "福德"/"凤辇"/"贵人(阳)" 等
    int32_t zhi;               // 0..11
} SxwnlLuckyHour;

#define SXWNL_ALMANAC_LUCKY_HOUR_MAX 8

// 用事择吉建议
typedef struct {
    char event[16];            // "动土"/"上梁"/"安床" 等
    bool suitable;             // 是否适合
    char reason[40];           // 一句话理由 "天赦日"/"杨公忌日"
} SxwnlEventAdvice;

#define SXWNL_ALMANAC_EVENT_MAX 8

// 文本数组 (宜忌/特别提示 用)
#define SXWNL_ALMANAC_TEXT_LIST_ITEM_MAX 5
#define SXWNL_ALMANAC_TEXT_LIST_LEN      16
typedef char SxwnlTextItem[SXWNL_ALMANAC_TEXT_LIST_LEN];
typedef char SxwnlNoteItem[80];
#define SXWNL_ALMANAC_NOTE_MAX 4

typedef struct {
    // ── 二十八宿 ──
    char xiu[8];               // 宿名 "角"
    char xiu_zheng[4];         // 七政 "木"
    char xiu_animal[8];        // 动物 "蛟"
    char xiu_luck[4];          // "吉" / "凶"
    char xiu_gong[4];          // 四象 "东" (东青龙/南朱雀/西白虎/北玄武)

    // ── 黄道黑道 ──
    char twelve_god[8];        // 12 神 "青龙"
    char huang_hei[8];         // "黄道" / "黑道"
    bool is_huang_dao;

    // ── 冲煞 ──
    char chong_sheng_xiao[8];  // "马"
    char chong_gan_zhi[8];     // "戊午"
    char sha[8];               // "南"

    // ── 五吉神方位 ──
    char xi_shen_fang[8];      // 喜神
    char yang_gui_fang[8];     // 阳贵神
    char yin_gui_fang[8];      // 阴贵神
    char fu_shen_fang[8];      // 福神
    char cai_shen_fang[8];     // 财神

    // ── 彭祖百忌 ──
    char peng_zu_gan[64];      // "甲不开仓 财物耗散"
    char peng_zu_zhi[64];      // "子不问卜 自惹祸殃"

    // ── 择日典籍语录 (董公择日要诀, 未来可加玉匣记/通胜经) ─
    SxwnlAlmanacQuote quotes[SXWNL_ALMANAC_QUOTE_MAX];
    int32_t quote_count;

    // ── 神煞 (按权重排序, 吉凶混编) ──────────────────────
    SxwnlShenSha shen_sha[SXWNL_ALMANAC_SHENSHA_MAX];
    int32_t shen_sha_count;

    // ── 宜忌 (神煞投票结果, 各最多 5 条) ─────────────────
    SxwnlTextItem yi[SXWNL_ALMANAC_TEXT_LIST_ITEM_MAX];
    int32_t yi_count;
    SxwnlTextItem ji[SXWNL_ALMANAC_TEXT_LIST_ITEM_MAX];
    int32_t ji_count;

    // ── 吉曜时法 (5吉时 + 2贵人时) ───────────────────────
    SxwnlLuckyHour lucky_hours[SXWNL_ALMANAC_LUCKY_HOUR_MAX];
    int32_t lucky_hour_count;

    // ── 用事择吉建议 (动土/竖柱/上梁/建屋/安灶/安床...) ──
    SxwnlEventAdvice events[SXWNL_ALMANAC_EVENT_MAX];
    int32_t event_count;

    // ── 特别提示 (节气/三煞方位...) ─────────────────────
    SxwnlNoteItem notes[SXWNL_ALMANAC_NOTE_MAX];
    int32_t note_count;
} SxwnlAlmanac;

// 取公历某日的老黄历完整数据. 返回 0 成功, -1 失败.
int sxwnl_get_almanac(int year, int month, int day, SxwnlAlmanac *out);

// 静态知识 (董公总论/口诀/方位 等, 与日历无关).
//   字段尺寸需为最长 UTF-8 字节 +1, 否则中文末尾会被截在多字节中间产生乱码.
typedef struct {
    char category[32];     // "总论"/"基础理论"/"建筑"/"口诀"/"方位"
    char title[64];        // 例 "乌兔太阳推算方法" (8 汉字 = 24 字节)
    char body[1024];
} SxwnlAlmanacTopic;

// 取所有 topic 条目, max_count 容量, 返回真实条数. 失败返 -1.
int sxwnl_get_almanac_topics(SxwnlAlmanacTopic *out, int max_count);

// ─── Eclipse / Astronomy ───

// Returns dynamically allocated string; caller must call sxwnl_string_free()
char* sxwnl_search_eclipses(int year, int month, int count, bool is_solar);
char* sxwnl_calc_eclipse_detail(int year, int month, int day,
                                 int hour, int minute, double second,
                                 bool is_utc, double longitude, double latitude);
char* sxwnl_calc_sun_moon_rise(int year, int month, int day,
                                double longitude, double latitude);

// ─── 日月升降/中天 (对应 JS RTS1, 用于月历下方信息栏) ───────
typedef struct {
    char sun_rise[16];         // 日出   HH:MM:SS
    char sun_set[16];          // 日落
    char sun_meridian[16];     // 日上中天
    char moon_rise[16];        // 月出   "--:--:--" 表示当日无升起
    char moon_set[16];         // 月落
    char moon_meridian[16];    // 月中天
    char civil_dawn[16];       // 晨起天亮 (民用晨, 太阳-6°)
    char civil_dusk[16];       // 晚上天黑 (民用昏)
    char day_length[16];       // 日照时间 (日出→日落, 与 JS sj 对应)
    char light_length[16];     // 白天时间 (民用晨→民用昏, 与 JS ch 对应)
} SxwnlDayRTS;

// 计算指定公历日在指定经纬度的日月升降时刻
//   longitude / latitude: 度, 东经/北纬为正
//   tz_hours:             所在时区偏移 (小时, 东向为正), 北京 = 8
int sxwnl_calc_day_rts(int year, int month, int day,
                       double longitude, double latitude,
                       double tz_hours, SxwnlDayRTS *out);

// ─── 地理目录 (GeoPostion + JWv/SQv 解码) ──────────────────
//   数据完全来自 libsxwnl 内部 (src/geo.cpp), 上层无需重复维护城市表.
typedef struct {
    char   province[64];      // "北京市"/"福建省"/...
    char   district[64];      // "天安门"/"福州市"/...
    double longitude;         // 度, 东+ 西-
    double latitude;          // 度, 北+ 南-
    double timezone;          // 小时, 东+ 西- (国内统一 8)
} SxwnlGeoCity;

typedef struct {
    char province[64];        // 省/直辖市/自治区名
    int  city_count;          // 该省下的城市/区县总数
} SxwnlGeoProvince;

typedef struct {
    char continent[32];       // "亚洲"/"欧洲"/...
    char country[96];         // "中国"/"加拿大东部时区"/...
    double timezone;          // 标准时偏移 (不含 DST)
    char cities[8][32];       // 至多 8 个代表城市名 (展示用)
    int  city_count;          // 实际填充的城市数
} SxwnlTimezoneGroup;

// 列出所有省/直辖市/自治区. 返回实际填充条数, 失败返 -1.
int sxwnl_geo_list_provinces(SxwnlGeoProvince *out, int out_max);
// 列某省内所有城市. 返回实际填充条数, 失败返 -1.
int sxwnl_geo_list_cities(const char *province, SxwnlGeoCity *out, int out_max);
// 按关键词模糊搜索区县/省份名. 返回实际填充条数, 失败返 -1.
int sxwnl_geo_search(const char *keyword, SxwnlGeoCity *out, int out_max);
// 列出国际时区分组. 返回实际填充条数, 失败返 -1.
int sxwnl_geo_list_timezone_groups(SxwnlTimezoneGroup *out, int out_max);
// 取默认地点(北京天安门). 失败返 -1.
int sxwnl_geo_default(SxwnlGeoCity *out);

void  sxwnl_string_free(char *str);

#ifdef __cplusplus
}
#endif

#endif // SXWNL_CAPI_H
