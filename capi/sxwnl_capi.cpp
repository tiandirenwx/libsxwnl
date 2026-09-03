#include "sxwnl_capi.h"
#include "sxtwl.h"
#include "day.h"
#include "bazi.h"
#include "bazi_analysis.h"
#include "geo.h"
#include "JD.h"
#include "SSQ.h"
#include "eph.h"
#include "eph_show.h"
#include "sx_lang_zh.h"
#include "year_utils.h"
#include "const.h"
#include "map_projection.h"
#include "world_map.h"
#include "star_eph.h"
#include "star_catalog.h"
#include "month_calendar.h"
#include "lunar_ob.h"
#include "lunar_month_name.h"
#include "eph_szj.h"
#include "almanac.h"

#include <array>
#include <cmath>
#include <cstring>
#include <memory>
#include <string>
#include <tuple>

// ─── internal helpers ───

static void safe_copy(char *dst, size_t dst_size, const char *src) {
    if (!src) { dst[0] = '\0'; return; }
    std::strncpy(dst, src, dst_size - 1);
    dst[dst_size - 1] = '\0';
}

static void safe_copy(char *dst, size_t dst_size, const std::string &src) {
    safe_copy(dst, dst_size, src.c_str());
}

// ─── exception guards ───
// C++ 异常不允许穿越 extern "C" 边界(否则 std::terminate / UB)。
// 所有会调入 C++ 库的导出函数都用以下 guard 兜底。
template <typename F>
static auto guard(F &&fn, decltype(fn()) fallback) noexcept -> decltype(fn()) {
    try { return fn(); } catch (...) { return fallback; }
}

template <typename F>
static void guard_void(F &&fn) noexcept {
    try { fn(); } catch (...) {}
}

static std::string gz_str(GZ gz) {
    return std::string(Gan[gz.tg]) + Zhi[gz.dz];
}

// 生成阴历月名 —— 统一委托 sxwnl::lunarMonthName(见 src/lunar_month_name.h),
// 与年历/月历/八字各处共用同一份逻辑。month 为月号(1-12, 与 Day::getLunarMonth() 一致),
// style 为 LunarMonthNameStyle(由 Day::getLunarMonthStyle() 给出)。
static std::string lunar_month_str(int year, int month, bool isLeap, int style) {
    return sxwnl::lunarMonthName(year, month, style, isLeap);
}

static std::string lunar_day_str(int day) {
    if (day >= 1 && day <= 30) return Rmc[day - 1];
    return std::to_string(day);
}

static void fill_day_info(Day *d, SxwnlDayInfo *out) {
    std::memset(out, 0, sizeof(SxwnlDayInfo));

    out->solar_year  = d->getSolarYear();
    out->solar_month = d->getSolarMonth();
    out->solar_day   = d->getSolarDay();
    out->lunar_year  = d->getLunarYear(true);
    out->lunar_month = d->getLunarMonth();
    out->lunar_day   = d->getLunarDay();
    out->is_leap_month = d->isLunarLeap();
    out->week_day    = d->getWeek();

    GZ yGZ = d->getYearGZ(false);
    GZ mGZ = d->getMonthGZ();
    GZ dGZ = d->getDayGZ();

    out->year_gan  = yGZ.tg;  out->year_zhi  = yGZ.dz;
    out->month_gan = mGZ.tg;  out->month_zhi = mGZ.dz;
    out->day_gan   = dGZ.tg;  out->day_zhi   = dGZ.dz;

    out->constellation = d->getConstellation();
    out->jian12 = d->getLunarDay12Jian();

    safe_copy(out->year_gz,  sizeof(out->year_gz),  gz_str(yGZ));
    safe_copy(out->month_gz, sizeof(out->month_gz), gz_str(mGZ));
    safe_copy(out->day_gz,   sizeof(out->day_gz),   gz_str(dGZ));
    safe_copy(out->lunar_month_name, sizeof(out->lunar_month_name),
              lunar_month_str(out->lunar_year, out->lunar_month,
                              out->is_leap_month, d->getLunarMonthStyle()));
    safe_copy(out->lunar_day_name, sizeof(out->lunar_day_name),
              lunar_day_str(out->lunar_day));
    safe_copy(out->sheng_xiao, sizeof(out->sheng_xiao), ShengXiao[yGZ.dz]);
    safe_copy(out->constellation_name, sizeof(out->constellation_name),
              XiZ[d->getConstellation()]);
    safe_copy(out->week_name, sizeof(out->week_name), WeekCn[d->getWeek()]);
    safe_copy(out->jian12_name, sizeof(out->jian12_name),
              RiJian12[d->getLunarDay12Jian()]);

    if (d->hasJieQi()) {
        out->jie_qi = d->getJieQi();
        safe_copy(out->jie_qi_name, sizeof(out->jie_qi_name), Jqmc[d->getJieQi()]);
        Time jqTime = JD::JD2DD(d->getJieQiJD());
        safe_copy(out->jie_qi_time, sizeof(out->jie_qi_time), JD::timeStr(jqTime));
    } else {
        out->jie_qi = -1;
        out->jie_qi_name[0] = '\0';
        out->jie_qi_time[0] = '\0';
    }

    // 历谱口径节气(整日表+QB, 对齐 sxwnl 网页版 ob.Ljq): 端上日历格子标签用此值,
    // 古代与天文口径 jie_qi 可能差 1 天; 精确时刻仍取 jie_qi_time(天文)。
    if (d->hasLiPuJieQi()) {
        out->lipu_jie_qi = d->getLiPuJieQi();
        safe_copy(out->lipu_jie_qi_name, sizeof(out->lipu_jie_qi_name),
                  Jqmc[d->getLiPuJieQi()]);
    } else {
        out->lipu_jie_qi = -1;
        out->lipu_jie_qi_name[0] = '\0';
    }

    if (d->hasYueXiang()) {
        out->yue_xiang = d->getYueXiang();
        safe_copy(out->yue_xiang_name, sizeof(out->yue_xiang_name),
                  YueXiangName[d->getYueXiang()]);
        // 与 jie_qi_time 格式保持一致(含完整日期+时间)
        Time yxTime = JD::JD2DD(d->getYueXiangJD());
        safe_copy(out->yue_xiang_time, sizeof(out->yue_xiang_time),
                  JD::timeStr(yxTime));
    } else {
        out->yue_xiang = -1;
        out->yue_xiang_name[0] = '\0';
        out->yue_xiang_time[0] = '\0';
    }

    // 节日 / 杂节 / 放假
    festival::DayInfo finfo = d->getFestivalInfo();
    safe_copy(out->holiday, sizeof(out->holiday), finfo.holiday);
    safe_copy(out->major,   sizeof(out->major),   finfo.major);
    safe_copy(out->minor,   sizeof(out->minor),   finfo.minor);
    safe_copy(out->misc,    sizeof(out->misc),    finfo.misc);
    out->is_off_day = finfo.isOffDay;

    // 农历纪年扩展
    out->lunar_jun_year    = d->getLunarYear(true);
    out->lunar_lichun_year = d->getLunarYear(false);
    out->huangdi_year      = d->getHuangdiYear();

    // 回历
    out->moslem_year  = d->getMoslemYear();
    out->moslem_month = d->getMoslemMonth();
    out->moslem_day   = d->getMoslemDay();

    // 纳音
    safe_copy(out->year_nayin,  sizeof(out->year_nayin),  bazi::naYin(yGZ.tg, yGZ.dz));
    safe_copy(out->month_nayin, sizeof(out->month_nayin), bazi::naYin(mGZ.tg, mGZ.dz));
    safe_copy(out->day_nayin,   sizeof(out->day_nayin),   bazi::naYin(dGZ.tg, dGZ.dz));

    // 儒略日(整数 + 小数), 与 JS 版 Bd0+J2000+jdF 等价
    out->julian_day = (double)d->getJulianDate();
}

static int days_in_month(int year, int month) {
    static const int table[] = {31,28,31,30,31,30,31,31,30,31,30,31};
    if (month < 1 || month > 12) return 0;
    int d = table[month - 1];
    if (month == 2 && ((year % 4 == 0 && year % 100 != 0) || year % 400 == 0))
        d = 29;
    return d;
}

// ─── Bazi handle: wraps BaziBase + cached strings ───

struct SxwnlBaziHandle {
    std::unique_ptr<BaziBase> bazi;
    std::string userName, gender, solarBirth, lunarBirth;
    std::string dateOfBirth, shengXiao, age, lifa;
    std::string dingQiType, ast, jieQi, qiYun, jiaoYun;

    void cacheStrings() {
        userName   = bazi->getUserName();
        gender     = bazi->getUserGender();
        solarBirth = bazi->getSolarBirth();
        shengXiao  = bazi->getShenXiao();
        age        = bazi->getAge();
        lifa       = bazi->getLifa();
        dingQiType = bazi->getDingQiType();  // 定气方式(历法), 用于简洁版"依据..."文案
        ast        = bazi->getAst();
        jieQi      = bazi->getJieQiTerms();   // 仅节气交接, 不含经纬度/真太阳时块
        qiYun      = bazi->getQiYun();
        jiaoYun    = bazi->getJiaoYun();
        // getLunarInfo() 返回 (年号, 生肖, 农历生日)
        auto lunarInfo = bazi->getLunarInfo();
        dateOfBirth = std::get<0>(lunarInfo);  // 出生年代(年号)
        lunarBirth  = std::get<2>(lunarInfo);  // 农历生日
    }
};

// ═══════════════════════════════════════════════════════════
//  Calendar API implementation
// ═══════════════════════════════════════════════════════════

int sxwnl_get_day_info(int year, int month, int day, SxwnlDayInfo *out) {
    return guard([&]() -> int {
        if (!out) return -1;
        std::unique_ptr<Day> d(sxtwl::fromSolar(year, (uint8_t)month, day));
        if (!d) return -1;
        fill_day_info(d.get(), out);
        return 0;
    }, -1);
}

// ═══════════════════════════════════════════════════════════
//  Almanac (老黄历)
// ═══════════════════════════════════════════════════════════
//
//  实现策略: libsxwnl 算干支/农历/星期等 -> 喂入 sxwnl::almanac 查表器.
//  almanac 模块本身不做任何天文计算, 严格保证 "干支唯一真源" 在 libsxwnl.

#include "almanac_topics.h"

int sxwnl_get_almanac_topics(SxwnlAlmanacTopic *out, int max_count) {
    return guard([&]() -> int {
        if (!out || max_count <= 0) return -1;
        const auto& topics = sxwnl::almanac::topics::kTopics;
        const int n = std::min(static_cast<int>(topics.size()), max_count);
        for (int i = 0; i < n; ++i) {
            safe_copy(out[i].category, sizeof(out[i].category), std::string(topics[i].category));
            safe_copy(out[i].title,    sizeof(out[i].title),    std::string(topics[i].title));
            safe_copy(out[i].body,     sizeof(out[i].body),     std::string(topics[i].body));
        }
        return n;
    }, -1);
}

int sxwnl_get_almanac(int year, int month, int day, SxwnlAlmanac *out) {
    return guard([&]() -> int {
        if (!out) return -1;
        std::unique_ptr<Day> d(sxtwl::fromSolar(year, (uint8_t)month, day));
        if (!d) return -1;

        const GZ yGZ = d->getYearGZ(false);
        const GZ mGZ = d->getMonthGZ();
        const GZ dGZ = d->getDayGZ();

        sxwnl::almanac::DayContext ctx;
        ctx.year_gan      = yGZ.tg;  ctx.year_zhi  = yGZ.dz;
        ctx.month_gan     = mGZ.tg;  ctx.month_zhi = mGZ.dz;
        ctx.day_gan       = dGZ.tg;  ctx.day_zhi   = dGZ.dz;
        ctx.week_day      = d->getWeek();
        ctx.lunar_month   = d->getLunarMonth();
        ctx.lunar_day     = d->getLunarDay();
        ctx.is_leap_month = d->isLunarLeap();
        ctx.julian_day    = static_cast<double>(d->getJulianDate());
        // 历谱口径(整日表+QB): 与日历"节气日"同源同日, 保证黄历特别提示与日历一致
        ctx.today_jie_qi  = d->hasLiPuJieQi() ? d->getLiPuJieQi() : -1;
        // 探查"明日"是否为节气: 用儒略日 +1 反推日期再查
        {
            int yn = year, mn = month, dn = day;
            // 简易日期加一天 (足够 1900-2100 区间)
            static const int kDaysPerMonth[13] = {0,31,28,31,30,31,30,31,31,30,31,30,31};
            auto isLeap = [](int y){ return (y%4==0 && y%100!=0) || y%400==0; };
            int dim = kDaysPerMonth[mn] + (mn == 2 && isLeap(yn) ? 1 : 0);
            if (++dn > dim) { dn = 1; if (++mn > 12) { mn = 1; ++yn; } }
            std::unique_ptr<Day> dn_obj(sxtwl::fromSolar(yn, (uint8_t)mn, dn));
            if (dn_obj && dn_obj->hasLiPuJieQi()) ctx.tomorrow_jie_qi = dn_obj->getLiPuJieQi();
        }

        const auto a = sxwnl::almanac::query(ctx);

        std::memset(out, 0, sizeof(*out));
        safe_copy(out->xiu,              sizeof(out->xiu),              a.xiu);
        safe_copy(out->xiu_zheng,        sizeof(out->xiu_zheng),        a.xiu_zheng);
        safe_copy(out->xiu_animal,       sizeof(out->xiu_animal),       a.xiu_animal);
        safe_copy(out->xiu_luck,         sizeof(out->xiu_luck),         a.xiu_luck);
        safe_copy(out->xiu_gong,         sizeof(out->xiu_gong),         a.xiu_gong);
        safe_copy(out->twelve_god,       sizeof(out->twelve_god),       a.twelve_god);
        safe_copy(out->huang_hei,        sizeof(out->huang_hei),        a.huang_hei);
        out->is_huang_dao = a.is_huang_dao;
        safe_copy(out->chong_sheng_xiao, sizeof(out->chong_sheng_xiao), a.chong_sheng_xiao);
        safe_copy(out->chong_gan_zhi,    sizeof(out->chong_gan_zhi),    a.chong_gan_zhi);
        safe_copy(out->sha,              sizeof(out->sha),              a.sha);
        safe_copy(out->xi_shen_fang,     sizeof(out->xi_shen_fang),     a.xi_shen_fang);
        safe_copy(out->yang_gui_fang,    sizeof(out->yang_gui_fang),    a.yang_gui_fang);
        safe_copy(out->yin_gui_fang,     sizeof(out->yin_gui_fang),     a.yin_gui_fang);
        safe_copy(out->fu_shen_fang,     sizeof(out->fu_shen_fang),     a.fu_shen_fang);
        safe_copy(out->cai_shen_fang,    sizeof(out->cai_shen_fang),    a.cai_shen_fang);
        safe_copy(out->peng_zu_gan,      sizeof(out->peng_zu_gan),      a.peng_zu_gan);
        safe_copy(out->peng_zu_zhi,      sizeof(out->peng_zu_zhi),      a.peng_zu_zhi);

        const int n_quote = std::min(
            static_cast<int>(a.quotes.size()),
            static_cast<int>(SXWNL_ALMANAC_QUOTE_MAX));
        for (int i = 0; i < n_quote; ++i) {
            const auto &q = a.quotes[i];
            safe_copy(out->quotes[i].source, sizeof(out->quotes[i].source), q.source);
            safe_copy(out->quotes[i].title,  sizeof(out->quotes[i].title),  q.title);
            safe_copy(out->quotes[i].luck,   sizeof(out->quotes[i].luck),   q.luck);
            safe_copy(out->quotes[i].body,   sizeof(out->quotes[i].body),   q.body);
        }
        out->quote_count = n_quote;

        // 神煞
        const int n_ss = std::min(
            static_cast<int>(a.shen_sha.size()),
            static_cast<int>(SXWNL_ALMANAC_SHENSHA_MAX));
        for (int i = 0; i < n_ss; ++i) {
            safe_copy(out->shen_sha[i].name, sizeof(out->shen_sha[i].name), a.shen_sha[i].name);
            out->shen_sha[i].is_lucky = a.shen_sha[i].is_lucky;
            out->shen_sha[i].weight   = a.shen_sha[i].weight;
        }
        out->shen_sha_count = n_ss;

        // 宜忌
        const int n_yi = std::min(
            static_cast<int>(a.yi.size()),
            static_cast<int>(SXWNL_ALMANAC_TEXT_LIST_ITEM_MAX));
        for (int i = 0; i < n_yi; ++i)
            safe_copy(out->yi[i], SXWNL_ALMANAC_TEXT_LIST_LEN, a.yi[i]);
        out->yi_count = n_yi;

        const int n_ji = std::min(
            static_cast<int>(a.ji.size()),
            static_cast<int>(SXWNL_ALMANAC_TEXT_LIST_ITEM_MAX));
        for (int i = 0; i < n_ji; ++i)
            safe_copy(out->ji[i], SXWNL_ALMANAC_TEXT_LIST_LEN, a.ji[i]);
        out->ji_count = n_ji;

        // 吉时
        const int n_lh = std::min(
            static_cast<int>(a.lucky_hours.size()),
            static_cast<int>(SXWNL_ALMANAC_LUCKY_HOUR_MAX));
        for (int i = 0; i < n_lh; ++i) {
            safe_copy(out->lucky_hours[i].name, sizeof(out->lucky_hours[i].name),
                      a.lucky_hours[i].name);
            out->lucky_hours[i].zhi = a.lucky_hours[i].zhi;
        }
        out->lucky_hour_count = n_lh;

        // 用事择吉
        const int n_ev = std::min(
            static_cast<int>(a.event_advices.size()),
            static_cast<int>(SXWNL_ALMANAC_EVENT_MAX));
        for (int i = 0; i < n_ev; ++i) {
            safe_copy(out->events[i].event,  sizeof(out->events[i].event),
                      a.event_advices[i].event);
            out->events[i].suitable = a.event_advices[i].suitable;
            safe_copy(out->events[i].reason, sizeof(out->events[i].reason),
                      a.event_advices[i].reason);
        }
        out->event_count = n_ev;

        // 特别提示
        const int n_no = std::min(
            static_cast<int>(a.notes.size()),
            static_cast<int>(SXWNL_ALMANAC_NOTE_MAX));
        for (int i = 0; i < n_no; ++i)
            safe_copy(out->notes[i], sizeof(out->notes[i]), a.notes[i]);
        out->note_count = n_no;

        return 0;
    }, -1);
}

int sxwnl_get_month_days(int year, int month, SxwnlDayInfo *out, int max_count) {
    return guard([&]() -> int {
        if (!out || month < 1 || month > 12) return -1;
        int total = days_in_month(year, month);
        if (total > max_count) total = max_count;

        for (int i = 0; i < total; i++) {
            std::unique_ptr<Day> d(sxtwl::fromSolar(year, (uint8_t)month, i + 1));
            if (!d) return -1;
            fill_day_info(d.get(), &out[i]);
        }
        return total;
    }, -1);
}

int sxwnl_lunar_to_solar(int year, int month, int day, bool is_leap,
                          int *out_year, int *out_month, int *out_day) {
    return guard([&]() -> int {
        std::unique_ptr<Day> d(sxtwl::fromLunar(year, (uint8_t)month, day, is_leap));
        if (!d) return -1;
        if (out_year)  *out_year  = d->getSolarYear();
        if (out_month) *out_month = d->getSolarMonth();
        if (out_day)   *out_day   = d->getSolarDay();
        return 0;
    }, -1);
}

int sxwnl_solar_to_lunar(int year, int month, int day, SxwnlDayInfo *out) {
    return sxwnl_get_day_info(year, month, day, out);
}

int sxwnl_solar_to_moslem(int year, int month, int day,
                          int *out_h_year, int *out_h_month, int *out_h_day) {
    return guard([&]() -> int {
        std::unique_ptr<Day> d(sxtwl::fromSolar(year, (uint8_t)month, day));
        if (!d) return -1;
        if (out_h_year)  *out_h_year  = d->getMoslemYear();
        if (out_h_month) *out_h_month = d->getMoslemMonth();
        if (out_h_day)   *out_h_day   = d->getMoslemDay();
        return 0;
    }, -1);
}

int sxwnl_moslem_to_solar(int h_year, int h_month, int h_day,
                          int *out_year, int *out_month, int *out_day) {
    return guard([&]() -> int {
        if (h_month < 1 || h_month > 12 || h_day < 1 || h_day > 30) return -1;
        std::unique_ptr<Day> d(Day::fromMoslem(h_year, (uint8_t)h_month, h_day));
        if (!d) return -1;
        if (out_year)  *out_year  = d->getSolarYear();
        if (out_month) *out_month = d->getSolarMonth();
        if (out_day)   *out_day   = d->getSolarDay();
        return 0;
    }, -1);
}

int sxwnl_get_day_info_by_moslem(int h_year, int h_month, int h_day,
                                 SxwnlDayInfo *out) {
    return guard([&]() -> int {
        if (!out) return -1;
        if (h_month < 1 || h_month > 12 || h_day < 1 || h_day > 30) return -1;
        std::unique_ptr<Day> d(Day::fromMoslem(h_year, (uint8_t)h_month, h_day));
        if (!d) return -1;
        fill_day_info(d.get(), out);
        return 0;
    }, -1);
}

int sxwnl_get_calendar_month(int year, int month,
                             SxwnlCalendarMonthHeader *out_header,
                             SxwnlDayInfo *out_days, int max_days) {
    return guard([&]() -> int {
        if (month < 1 || month > 12) return -1;
        sxwnl::MonthCalendar cal(year, month);
        const auto &dv = cal.days();
        int dn = (int)dv.size();

        if (out_header) {
            out_header->year             = cal.year();
            out_header->month            = cal.month();
            out_header->first_week_day   = cal.firstWeek();
            out_header->day_count        = dn;
            out_header->total_weeks      = cal.totalWeeks();
            out_header->first_julian_day = cal.firstJulianDay();
            out_header->year_gan         = cal.yearGZ().tg;
            out_header->year_zhi         = cal.yearGZ().dz;
            std::string yg = std::string(Gan[cal.yearGZ().tg]) + Zhi[cal.yearGZ().dz];
            safe_copy(out_header->year_gz,    sizeof(out_header->year_gz),    yg);
            safe_copy(out_header->sheng_xiao, sizeof(out_header->sheng_xiao), cal.shengXiao());
            safe_copy(out_header->nianhao,    sizeof(out_header->nianhao),    cal.nianHao());
        }

        if (out_days && max_days > 0) {
            int n = std::min(dn, max_days);
            for (int i = 0; i < n; ++i) {
                // 直接复用 MonthCalendar 内部已经构造好的 Day
                fill_day_info(dv[i].day.get(), &out_days[i]);
            }
        }
        return dn;
    }, -1);
}

int sxwnl_get_year_leap_month(int year) {
    return guard([&] { return sxtwl::getRunMonth(year); }, 0);
}

int sxwnl_get_lunar_month_days(int year, int month, bool is_leap, bool is_spec) {
    return guard([&] { return (int)sxtwl::getLunarMonthNum(year, (uint8_t)month, is_leap, is_spec); }, 0);
}

int sxwnl_get_lunar_day_name(int day, char *out, int cap) {
    return guard([&]() -> int {
        if (!out || cap <= 0) return 0;
        std::string name = lunar_day_str(day);   // 复用 Rmc[] 表, 越界回退数字
        safe_copy(out, (size_t)cap, name);
        return (int)std::strlen(out);
    }, 0);
}

int sxwnl_get_solar_month_valid_days(int year, int month, int *out, int max_count) {
    return guard([&]() -> int {
        if (!out || max_count <= 0 || month < 1 || month > 12) return 0;
        int count = 0;
        // 用本库历法做"往返校验": 不存在的日期(改历缺口/超出月长)换算回来 != 原值。
        for (int d = 1; d <= 31 && count < max_count; d++) {
            Time t = JD::JD2DD(JD::DD2JD(year, (uint8_t)month, (long double)d + 0.5L));
            if (t.Y == year && t.M == month && t.D == d) out[count++] = d;
        }
        return count;
    }, 0);
}

// 枚举某节气年窗口的逐月(原始顺序), 名称含闰/特殊/后九月
static void enum_lunar_window(int year, std::vector<SxwnlLunarMonth> &v) {
    SSQ &ssq = SSQ::getInstance();
    ssq.calcY(int2((year - 2000.0) * 365.2422 + 180));
    std::vector<long double> vecZQ = ssq.getZhongQi();
    std::vector<int> vecHS = ssq.getHS();
    std::vector<bool> vecSpec = ssq.getSpecificLunarMonth();
    std::vector<int> vecYueName = ssq.getYueName();
    std::vector<int> vecStyle = ssq.getMonthDisplayStyle();
    int leap = ssq.getLeap();

    if (vecZQ.size() <= 24) return;  // 需要 vecZQ[24] 作为窗口上界

    for (int i = 0; i < 14; i++) {
        if ((int)vecHS.size() <= i + 1 || vecHS[i + 1] > vecZQ[24]) break;
        // 并行数组边界保护
        if (i >= (int)vecSpec.size() || i >= (int)vecYueName.size() ||
            i >= (int)vecStyle.size()) break;

        bool isLeap = (leap && i == leap);
        int yn = vecYueName[i];
        // 统一走 sxwnl::lunarMonthName: 后九月/十三月(闰)、连续同名月(SYmc)、
        // 689-701 建寅(一月) 全部由 (year, yn, style, isLeap) 决定
        std::string name = sxwnl::lunarMonthName(year, yn, vecStyle[i], isLeap);

        SxwnlLunarMonth e{};
        e.month = yn;
        e.is_leap = isLeap ? 1 : 0;
        e.is_spec = vecSpec[i] ? 1 : 0;
        safe_copy(e.name, sizeof(e.name), name);
        v.push_back(e);
    }
}

int sxwnl_get_lunar_months(int year, SxwnlLunarMonth *out, int max_count) {
    return guard([&]() -> int {
        if (!out || max_count <= 0) return 0;
        // 农历年Y: 正..十月取自窗口Y, 冬腊月取自窗口Y+1(因 month>10 内部 year+1)
        std::vector<SxwnlLunarMonth> wy, wy1;
        enum_lunar_window(year, wy);
        enum_lunar_window(year + 1, wy1);

        int count = 0;
        for (int m = 1; m <= 12 && count < max_count; m++) {
            const std::vector<SxwnlLunarMonth> &src = (m <= 10) ? wy : wy1;
            // 先正常月, 再闰月(紧随其后)
            for (int pass = 0; pass < 2 && count < max_count; pass++) {
                for (size_t i = 0; i < src.size() && count < max_count; i++) {
                    if (src[i].month != m) continue;
                    if (pass == 0 && src[i].is_leap) continue;
                    if (pass == 1 && !src[i].is_leap) continue;
                    out[count++] = src[i];
                }
            }
        }
        return count;
    }, 0);
}

// ═══════════════════════════════════════════════════════════
//  JieQi API implementation
// ═══════════════════════════════════════════════════════════

int sxwnl_get_jieqi_list(int year, SxwnlJieQiItem *out, int max_count) {
    return guard([&]() -> int {
        if (!out) return -1;
        int count = 0;

        for (int m = 1; m <= 12 && count < max_count; m++) {
            int total = days_in_month(year, m);
            for (int d = 1; d <= total && count < max_count; d++) {
                std::unique_ptr<Day> day(sxtwl::fromSolar(year, (uint8_t)m, d));
                if (day && day->hasJieQi()) {
                    SxwnlJieQiItem &item = out[count];
                    std::memset(&item, 0, sizeof(item));
                    item.idx = day->getJieQi();
                    item.solar_month = m;
                    item.solar_day = d;
                    safe_copy(item.name, sizeof(item.name), Jqmc[day->getJieQi()]);
                    Time jqTime = JD::JD2DD(day->getJieQiJD());
                    safe_copy(item.time, sizeof(item.time), JD::timeStr(jqTime));
                    count++;
                }
            }
        }
        return count;
    }, -1);
}

// ═══════════════════════════════════════════════════════════
//  Year Calendar API: 按农历月聚合 (对应 sxwnl/lunar.js nianLi2HTML)
//   - 每月输出: 月名/大小/闰、朔日干支与公历日期、本月范围内出现的节气
//   - 实现思路: 复用 SSQ 单例的 HS(合朔)/ZQ(中气节气)/ZQpe(补气) 表,
//     再借助 sxtwl::fromSolar + getDayGZ 取干支, 避免重复实现历法逻辑
// ═══════════════════════════════════════════════════════════

int sxwnl_get_year_calendar(int year, SxwnlYearCalMonth *out, int max_count) {
    return guard([&]() -> int {
        if (!out || max_count <= 0) return -1;

        SSQ &ssq = SSQ::getInstance();
        const int seedJd = int2((year - 2000) * 365.2422 + 180);
        ssq.calcY(seedJd);

        const auto HS  = ssq.getHS();          // 15 个合朔 (J2000 起)
        const auto ZQ  = ssq.getZhongQi();     // 25 个气
        const auto Pe  = ssq.getZQPe();        // 2 个补气 (置于 ZQ[0] 之前)
        const auto Ym  = ssq.getYm();          // 14 个月建序号 (0=子/冬月, 1=丑/腊月, 2=寅/正月, ..., 11=亥/十月)
        const auto Dx  = ssq.getDx();          // 14 个月大小
        const auto Spc = ssq.getSpecificLunarMonth();
        const auto Sty = ssq.getMonthDisplayStyle(); // 14 个月显示风格
        const int  leap = ssq.getLeap();

        if (HS.size() < 15 || ZQ.size() < 25 || Pe.size() < 2)
            return -1;

        const long double zqEnd = ZQ[24];      // 下一年冬至, 用作上限
        int count = 0;

        for (int i = 0; i < (int)HS.size() - 1 && count < max_count; ++i) {
            if (HS[i + 1] > zqEnd) break;      // 已含下一年冬至, 停止

            SxwnlYearCalMonth &m = out[count];
            std::memset(&m, 0, sizeof(m));

            // SSQ.getYm() 给出月建序号 mc, 而 Ymc[]/SYmc[] 的索引是 (月号 - 1):
            //   Ymc[0]="正"(寅) Ymc[1]="二" ... Ymc[9]="十" Ymc[10]="冬" Ymc[11]="腊"
            // 故月建→月名索引为 ymi = (mc + 10) mod 12
            int mc  = (i < (int)Ym.size()) ? Ym[i] : 0;
            int ymi = ((mc % 12) + 10) % 12;
            bool isLeap = (leap && i == leap);
            bool isSpec = (i < (int)Spc.size()) ? Spc[i] : false;
            int style = (i < (int)Sty.size()) ? Sty[i] : LUNAR_MONTH_NORMAL;
            // 统一走 sxwnl::lunarMonthName: 后九月/十三月(闰)、连续同名月(SYmc)、
            // 689-701 建寅(一月) 全部由 (year, ymi+1, style, isLeap) 决定
            std::string monthName = sxwnl::lunarMonthName(year, ymi + 1, style, isLeap);

            m.month_idx = ymi;
            m.is_leap   = isLeap ? 1 : 0;
            m.is_spec   = isSpec ? 1 : 0;
            m.day_count = (i < (int)Dx.size()) ? Dx[i] : 0;
            safe_copy(m.month_name, sizeof(m.month_name), monthName);

            // 朔日 → 公历日期 + 干支. HS[i] 是 J2000 相对整数日, +J2000 得绝对 JD
            Time shuoT = JD::JD2DD((long double)HS[i] + J2000);
            m.solar_year  = shuoT.Y;
            m.solar_month = shuoT.M;
            m.solar_day   = (int)shuoT.D;

            std::unique_ptr<Day> shuoDay(sxtwl::fromSolar(shuoT.Y, (uint8_t)shuoT.M, (int)shuoT.D));
            if (shuoDay) {
                GZ g = shuoDay->getDayGZ();
                safe_copy(m.shuo_gz, sizeof(m.shuo_gz), gz_str(g));
            }

            // 本月范围 [HS[i], HS[i+1]) 内的节气. ZQ/Pe 同样是 J2000 相对值
            int jqCount = 0;
            constexpr int kMaxJq = sizeof(m.jieqi) / sizeof(m.jieqi[0]);
            for (int j = -2; j < 24 && jqCount < kMaxJq; ++j) {
                long double qiJdRel;   // J2000 相对儒略日(整数)
                int qiIdx;
                if (j >= 0) {
                    qiJdRel = ZQ[j];
                    qiIdx   = j;        // 0=冬至..23=大雪
                } else if (j == -1) {
                    qiJdRel = (long double)Pe[0];
                    qiIdx   = 23;       // 大雪 (前一年)
                } else { // j == -2
                    qiJdRel = (long double)Pe[1];
                    qiIdx   = 22;       // 小雪
                }

                if (qiJdRel < (long double)HS[i] || qiJdRel >= (long double)HS[i + 1])
                    continue;

                SxwnlYearCalJieQi &jq = m.jieqi[jqCount++];
                std::memset(&jq, 0, sizeof(jq));

                Time qt = JD::JD2DD(qiJdRel + J2000);
                jq.idx         = qiIdx;
                jq.solar_month = qt.M;
                jq.solar_day   = (int)qt.D;
                jq.day_offset  = (int)qiJdRel - HS[i];   // 距月首天数
                safe_copy(jq.name, sizeof(jq.name), Jqmc[qiIdx]);
                if (jq.day_offset >= 0 && jq.day_offset < 30) {
                    safe_copy(jq.day_name, sizeof(jq.day_name), Rmc[jq.day_offset]);
                }

                std::unique_ptr<Day> qiDay(sxtwl::fromSolar(qt.Y, (uint8_t)qt.M, (int)qt.D));
                if (qiDay) {
                    GZ g = qiDay->getDayGZ();
                    safe_copy(jq.gz, sizeof(jq.gz), gz_str(g));
                }
                // 精确交气时刻: 直接对历谱气日做 qi_accurate2 精化(对齐 sxwnl
                // nianLiHTML / LunarYear::getNianLiStr), 而非依赖 Day::hasJieQi().
                //   古代(1645 年前)平气历谱日与天文定气日可能差 1 天, 此时
                //   hasJieQi() 在历谱日上为 false, 会漏掉时刻; 直接精化则恒有值.
                long double qiAcc = qi_accurate2(qiJdRel);  // J2000 相对精确交气(北京时)
                Time accT = JD::JD2DD(qiAcc + (long double)J2000);
                jq.acc_month = accT.M;
                jq.acc_day   = (int)accT.D;
                safe_copy(jq.time, sizeof(jq.time),
                          JD::timeStr(qiAcc + (long double)J2000));
            }
            m.jq_count = jqCount;
            ++count;
        }
        return count;
    }, -1);
}

// ═══════════════════════════════════════════════════════════
//  Bazi API implementation
// ═══════════════════════════════════════════════════════════

SxwnlBazi sxwnl_bazi_create(const SxwnlBaziParams *params) {
    if (!params) return nullptr;

    SBaziInputPara input{};
    input.name     = params->name;
    input.gender   = params->gender;
    input.isAst    = params->is_ast;
    input.calendar = params->is_lunar ? CalendarLunar : CalendarSolar;
    input.isRun    = params->is_leap;
    input.isSpec   = params->is_spec;
    input.lifa     = static_cast<LiFaType>(params->lifa);
    if (input.lifa == LifaUnknown) input.lifa = XianDaiNongLifa_DingQiFa;

    input.birthDayTime = Time(params->year, params->month, params->day,
                              params->hour, params->minute, params->second);
    input.jw = {
        params->longitude > 0 ? params->longitude : 116.3833333,
        params->latitude > 0  ? params->latitude  : 39.9,
        "", ""
    };

    // extern "C" 边界: 不允许 C++ 异常逃逸(否则 std::terminate);
    // 用 unique_ptr 保证抛出时不泄漏 handle。
    try {
        auto handle = std::make_unique<SxwnlBaziHandle>();
        handle->bazi = std::make_unique<BaziBase>(input);
        handle->bazi->calcBaziPaiPan();
        handle->cacheStrings();
        return handle.release();
    } catch (...) {
        return nullptr;
    }
}

void sxwnl_bazi_free(SxwnlBazi bazi) {
    delete static_cast<SxwnlBaziHandle*>(bazi);
}

#define BAZI_GETTER(FUNC, FIELD) \
    const char* FUNC(SxwnlBazi bazi) { \
        auto *h = static_cast<SxwnlBaziHandle*>(bazi); \
        return h ? h->FIELD.c_str() : ""; \
    }

BAZI_GETTER(sxwnl_bazi_get_user_name,    userName)
BAZI_GETTER(sxwnl_bazi_get_gender,        gender)
BAZI_GETTER(sxwnl_bazi_get_solar_birth,   solarBirth)
BAZI_GETTER(sxwnl_bazi_get_lunar_birth,   lunarBirth)
BAZI_GETTER(sxwnl_bazi_get_date_of_birth, dateOfBirth)
BAZI_GETTER(sxwnl_bazi_get_sheng_xiao,    shengXiao)
BAZI_GETTER(sxwnl_bazi_get_age,           age)
BAZI_GETTER(sxwnl_bazi_get_lifa,          lifa)
BAZI_GETTER(sxwnl_bazi_get_ding_qi_type,  dingQiType)
BAZI_GETTER(sxwnl_bazi_get_ast,           ast)
BAZI_GETTER(sxwnl_bazi_get_jie_qi,        jieQi)
BAZI_GETTER(sxwnl_bazi_get_qi_yun,        qiYun)
BAZI_GETTER(sxwnl_bazi_get_jiao_yun,      jiaoYun)

#undef BAZI_GETTER

int sxwnl_bazi_get_si_zhu(SxwnlBazi bazi, SxwnlPillar out[4]) {
    return guard([&]() -> int {
        auto *h = static_cast<SxwnlBaziHandle*>(bazi);
        if (!h) return -1;

        auto szData = h->bazi->getSiZhu();
        size_t cols = szData.empty() ? 0 : szData[0].size();

        std::memset(out, 0, sizeof(SxwnlPillar) * 4);
        // szData columns: [0]=header/label, [1]=year, [2]=month, [3]=day, [4]=hour
        // 逐行校验列数(各行长度未必一致), 避免越界读
        auto cell = [&](size_t row, size_t col) -> const std::string& {
            static const std::string kEmpty;
            if (row < szData.size() && col < szData[row].size()) return szData[row][col];
            return kEmpty;
        };
        for (size_t c = 0; c < 4 && c + 1 < cols; c++) {
            safe_copy(out[c].shi_shen,  8, cell(0, c + 1));
            safe_copy(out[c].tian_gan,  8, cell(1, c + 1));
            safe_copy(out[c].di_zhi,    8, cell(2, c + 1));
            safe_copy(out[c].cang_gan1, 8, cell(3, c + 1));
            safe_copy(out[c].cang_gan2, 8, cell(4, c + 1));
            safe_copy(out[c].cang_gan3, 8, cell(5, c + 1));
        }
        return 4;
    }, -1);
}

int sxwnl_bazi_get_da_yun_count(SxwnlBazi bazi) {
    return guard([&]() -> int {
        auto *h = static_cast<SxwnlBaziHandle*>(bazi);
        if (!h) return 0;
        return (int)h->bazi->getDaYunList().size();
    }, 0);
}

int sxwnl_bazi_get_da_yun(SxwnlBazi bazi, SxwnlDaYun *out, int max_count) {
    return guard([&]() -> int {
        auto *h = static_cast<SxwnlBaziHandle*>(bazi);
        if (!h || !out) return -1;

        auto daYun  = h->bazi->getDaYunList();
        auto starts = h->bazi->getStartYearList();
        auto ends   = h->bazi->getEndYearList();

        int count = (int)daYun.size();
        if (count > max_count) count = max_count;

        for (int i = 0; i < count; i++) {
            std::memset(&out[i], 0, sizeof(SxwnlDaYun));
            safe_copy(out[i].gan_zhi, sizeof(out[i].gan_zhi), daYun[i]);
            out[i].start_year = (i < (int)starts.size()) ? starts[i] : 0;
            out[i].end_year   = (i < (int)ends.size())   ? ends[i]   : 0;
        }
        return count;
    }, -1);
}

int sxwnl_bazi_get_fleeting_years(SxwnlBazi bazi, char (*out)[8], int max_count) {
    return guard([&]() -> int {
        auto *h = static_cast<SxwnlBaziHandle*>(bazi);
        if (!h || !out) return -1;

        auto years = h->bazi->getFleetingYearList(false);
        int count = (int)years.size();
        if (count > max_count) count = max_count;

        for (int i = 0; i < count; i++) {
            safe_copy(out[i], 8, years[i]);
        }
        return count;
    }, -1);
}

void sxwnl_get_hour_gz(int day_gan, int hour, char out[8]) {
    guard_void([&] {
        GZ gz = sxtwl::getShiGz((uint8_t)day_gan, (uint8_t)hour);
        safe_copy(out, 8, gz_str(gz));
    });
}

// ═══════════════════════════════════════════════════════════
//  Bazi 完整排盘列信息
// ═══════════════════════════════════════════════════════════

// 填充一柱完整信息. main_star 为主星文字(日柱传性别). start_year 仅大运/流年用
static void fill_column(SxwnlBaziColumn *c, const std::array<int, 8> &sz,
                        int day_gan, int gan, int zhi,
                        const std::string &main_star, int start_year) {
    std::memset(c, 0, sizeof(*c));
    safe_copy(c->gan, sizeof(c->gan), Gan[gan]);
    safe_copy(c->zhi, sizeof(c->zhi), Zhi[zhi]);
    safe_copy(c->gan_shi_shen, sizeof(c->gan_shi_shen), main_star);
    safe_copy(c->nayin, sizeof(c->nayin), bazi::naYin(gan, zhi));
    safe_copy(c->xing_yun, sizeof(c->xing_yun), bazi::changSheng(day_gan, zhi));
    safe_copy(c->zi_zuo, sizeof(c->zi_zuo), bazi::changSheng(gan, zhi));
    safe_copy(c->kong_wang, sizeof(c->kong_wang), bazi::kongWangStr(gan, zhi));

    auto cg = bazi::cangGan(zhi);
    c->cang_gan_count = (int)cg.size();
    for (size_t i = 0; i < cg.size() && i < 3; i++) {
        safe_copy(c->cang_gan[i], 8, Gan[cg[i]]);
        safe_copy(c->cang_gan_shi_shen[i], 8, bazi::shiShenShort(day_gan, cg[i]));
    }

    auto ss = bazi::shenSha(sz, gan, zhi);
    int n = (int)ss.size();
    if (n > 12) n = 12;
    c->shen_sha_count = n;
    for (int i = 0; i < n; i++) safe_copy(c->shen_sha[i], 20, ss[i]);

    c->start_year = start_year;
}

// 计算大运干支索引(8步)及对应起始年
static void calc_da_yun(BaziBase *b, std::array<int, 8> sz,
                        std::vector<std::pair<int, int>> &gz, std::vector<int> &startYears) {
    bool female = b->getGenderIsFemale();
    bool forward = !((sz[0] % 2) ^ (female ? 1 : 0)); // 阳男阴女顺行
    int dGan = forward ? 1 : 9;   // +1 顺 / -1(=+9) 逆
    int dZhi = forward ? 1 : 11;
    int j = sz[2], k = sz[3];
    int startYear = b->getStartYear();
    for (int i = 0; i < 8; i++) {
        j = (j + dGan) % 10;
        k = (k + dZhi) % 12;
        gz.push_back({j, k});
        startYears.push_back(startYear + i * 10);
    }
}

// 流年干支(以立春为界, 公历年Y)
static std::pair<int, int> liu_nian_gz(int year) {
    int g = ((year - 4) % 10 + 10) % 10;
    int z = ((year - 4) % 12 + 12) % 12;
    return {g, z};
}

int sxwnl_bazi_get_columns(SxwnlBazi bazi, SxwnlBaziColumn out[4]) {
    return guard([&]() -> int {
        auto *h = static_cast<SxwnlBaziHandle *>(bazi);
        if (!h || !out) return -1;
        auto sz = h->bazi->getSiZhuIndex();
        int dayGan = sz[4];
        bool female = h->bazi->getGenderIsFemale();

        const int ganIdx[4] = {0, 2, 4, 6};
        const int zhiIdx[4] = {1, 3, 5, 7};
        for (int p = 0; p < 4; p++) {
            int gan = sz[ganIdx[p]], zhi = sz[zhiIdx[p]];
            std::string mainStar;
            if (p == 2) mainStar = female ? "女" : "男";  // 日柱标性别(日主)
            else mainStar = bazi::shiShenShort(dayGan, gan);
            fill_column(&out[p], sz, dayGan, gan, zhi, mainStar, 0);
        }
        return 4;
    }, -1);
}

int sxwnl_bazi_get_da_yun_columns(SxwnlBazi bazi, SxwnlBaziColumn *out, int max_count) {
    return guard([&]() -> int {
        auto *h = static_cast<SxwnlBaziHandle *>(bazi);
        if (!h || !out) return -1;
        auto sz = h->bazi->getSiZhuIndex();
        int dayGan = sz[4];

        std::vector<std::pair<int, int>> gz;
        std::vector<int> starts;
        calc_da_yun(h->bazi.get(), sz, gz, starts);

        int count = (int)gz.size();
        if (count > max_count) count = max_count;
        for (int i = 0; i < count; i++) {
            std::string star = bazi::shiShenShort(dayGan, gz[i].first);
            fill_column(&out[i], sz, dayGan, gz[i].first, gz[i].second, star, starts[i]);
        }
        return count;
    }, -1);
}

int sxwnl_bazi_get_current_da_yun(SxwnlBazi bazi, SxwnlBaziColumn *out) {
    return guard([&]() -> int {
        auto *h = static_cast<SxwnlBaziHandle *>(bazi);
        if (!h || !out) return -1;
        auto sz = h->bazi->getSiZhuIndex();
        int dayGan = sz[4];
        int curYear = JD::getNowTime().getYear();

        std::vector<std::pair<int, int>> gz;
        std::vector<int> starts;
        calc_da_yun(h->bazi.get(), sz, gz, starts);
        if (gz.empty()) return -1;

        int sel = 0;
        if (curYear < starts[0]) sel = 0; // 未起运: 取首步
        else {
            sel = (int)gz.size() - 1;
            for (int i = 0; i < (int)gz.size(); i++)
                if (curYear >= starts[i] && curYear < starts[i] + 10) { sel = i; break; }
        }
        std::string star = bazi::shiShenShort(dayGan, gz[sel].first);
        fill_column(out, sz, dayGan, gz[sel].first, gz[sel].second, star, starts[sel]);
        return 0;
    }, -1);
}

int sxwnl_bazi_get_current_liu_nian(SxwnlBazi bazi, SxwnlBaziColumn *out) {
    return guard([&]() -> int {
        auto *h = static_cast<SxwnlBaziHandle *>(bazi);
        if (!h || !out) return -1;
        auto sz = h->bazi->getSiZhuIndex();
        int dayGan = sz[4];
        int curYear = JD::getNowTime().getYear();
        auto gz = liu_nian_gz(curYear);
        std::string star = bazi::shiShenShort(dayGan, gz.first);
        fill_column(out, sz, dayGan, gz.first, gz.second, star, curYear);
        return 0;
    }, -1);
}

const char *sxwnl_bazi_get_si_ling(SxwnlBazi bazi) {
    return guard([&]() -> const char * {
        auto *h = static_cast<SxwnlBaziHandle *>(bazi);
        if (!h) return "";
        auto sz = h->bazi->getSiZhuIndex();
        int g = bazi::siLing(sz[3], h->bazi->getDaysAfterJie());
        return (g >= 0) ? Gan[g] : "";
    }, "");
}

void sxwnl_bazi_get_wuxing_count(SxwnlBazi bazi, int out[5], bool include_cang_gan) {
    guard_void([&] {
        if (!out) return;
        auto *h = static_cast<SxwnlBaziHandle *>(bazi);
        if (!h) { for (int i = 0; i < 5; i++) out[i] = 0; return; }
        auto cnt = bazi::wuXingCount(h->bazi->getSiZhuIndex(), include_cang_gan);
        for (int i = 0; i < 5; i++) out[i] = cnt[i];
    });
}

void sxwnl_bazi_get_wuxing_status(SxwnlBazi bazi, char out[5][8]) {
    guard_void([&] {
        auto *h = static_cast<SxwnlBaziHandle *>(bazi);
        if (!h || !out) return;
        auto st = bazi::wuXingStatus(h->bazi->getSiZhuIndex()[3]);
        for (int i = 0; i < 5; i++) safe_copy(out[i], 8, st[i]);
    });
}

int sxwnl_bazi_get_liu_nian(SxwnlBazi bazi, int start_year,
                            SxwnlLiuNianItem *out, int max_count) {
    return guard([&]() -> int {
    auto *h = static_cast<SxwnlBaziHandle *>(bazi);
    if (!h || !out) return -1;
    auto sz = h->bazi->getSiZhuIndex();
    int dayGan = sz[4];
    bool female = h->bazi->getGenderIsFemale();
    bool forward = !((sz[0] % 2) ^ (female ? 1 : 0));
    int firstYear = h->bazi->getStartYear(); // 1虚岁对应公历年(birth year≈起运基准)

    // 小运自时柱起, 顺/逆与大运同向; 1虚岁=时柱进一位
    int hourGan = sz[6], hourZhi = sz[7];
    int birthYear = h->bazi->getSolarYearOfBirth();

    int count = 0;
    for (int i = 0; i < max_count; i++) {
        int year = start_year + i;
        SxwnlLiuNianItem &it = out[count];
        std::memset(&it, 0, sizeof(it));
        it.year = year;
        int age = year - birthYear + 1; // 虚岁(出生年为1岁)
        it.age = age;

        auto gz = liu_nian_gz(year);
        safe_copy(it.gan_zhi, 8, std::string(Gan[gz.first]) + Zhi[gz.second]);
        safe_copy(it.gan_shi_shen, 8, bazi::shiShenShort(dayGan, gz.first));
        int bq = bazi::benQiGan(gz.second);
        if (bq >= 0)
            safe_copy(it.zhi_shi_shen, 8, bazi::shiShenShort(dayGan, bq));

        // 小运: 步数=该年虚岁(时柱进一位为1虚岁)
        int step = age > 0 ? age : 1;
        int dirZhi = forward ? step : (-step);
        int xg = ((hourGan + dirZhi) % 10 + 10) % 10;
        int xz = ((hourZhi + dirZhi) % 12 + 12) % 12;
        safe_copy(it.xiao_yun, 8, std::string(Gan[xg]) + Zhi[xz]);
        safe_copy(it.xiao_yun_shi_shen, 8, bazi::shiShenShort(dayGan, xg));

        (void)firstYear;
        count++;
    }
    return count;
    }, -1);
}

// ═══════════════════════════════════════════════════════════
//  四柱反推
// ═══════════════════════════════════════════════════════════

int sxwnl_bazi_reverse(const int sz[8], int start_year, int end_year,
                       SxwnlReverseItem *out, int max_count) {
    return guard([&]() -> int {
    if (!sz || !out || max_count <= 0) return 0;
    int yg = sz[0], yz = sz[1], mg = sz[2], mz = sz[3];
    int dg = sz[4], dz = sz[5], hg = sz[6], hz = sz[7];
    bool hourKnown = (hg >= 0 && hz >= 0);

    int nDay = bazi::jiaZiIndex(dg, dz);
    int nYear = bazi::jiaZiIndex(yg, yz);
    if (nDay < 0 || nYear < 0) return 0;

    // 日柱对应 d0 残差: D=d0-6+9000000, D%60==nDay, 9000000%60==0 → d0 ≡ nDay+6 (mod 60)
    int rDay = ((nDay + 6) % 60 + 60) % 60;

    int count = 0;
    // 以年柱(立春界)锁定候选年: (Y-1984)%60 == nYear (mod 60)
    int yStep = ((nYear - (start_year - 1984)) % 60 + 60) % 60;
    int firstY = start_year + yStep;

    for (int Y = firstY; Y <= end_year && count < max_count; Y += 60) {
        // 该立春年覆盖区间约 [Y-01-01, (Y+1)-02-20]
        long double jd1 = JD::toJD(Time(Y, 1, 1, 12, 0, 0));
        long double jd2 = JD::toJD(Time(Y + 1, 2, 20, 12, 0, 0));
        int d0a = int2(jd1) - J2000;
        int d0b = int2(jd2) - J2000;
        // 对齐到日柱残差
        int d0 = d0a + (((rDay - (d0a % 60)) % 60) + 60) % 60;

        for (; d0 <= d0b && count < max_count; d0 += 60) {
            Time t = JD::JD2DD((long double)(d0 + J2000) + 0.0001);
            std::unique_ptr<Day> day(sxtwl::fromSolar(t.Y, (uint8_t)t.M, t.D));
            if (!day) continue;

            GZ gzY = day->getYearGZ(false);
            if (gzY.tg != yg || gzY.dz != yz) continue;
            GZ gzM = day->getMonthGZ();
            if (gzM.tg != mg || gzM.dz != mz) continue;
            GZ gzD = day->getDayGZ();
            if (gzD.tg != dg || gzD.dz != dz) continue;

            int hour = -1;
            if (hourKnown) {
                bool ok = false;
                for (int h = 0; h < 24; h++) {
                    GZ s = sxtwl::getShiGz(gzD.tg, (uint8_t)h);
                    if (s.tg == hg && s.dz == hz) { hour = h; ok = true; break; }
                }
                if (!ok) continue; // 该日不可能出现此时干支
            }

            SxwnlReverseItem &it = out[count];
            std::memset(&it, 0, sizeof(it));
            it.year = t.Y; it.month = t.M; it.day = t.D; it.hour = hour;
            safe_copy(it.ganzhi[0], 8, std::string(Gan[gzY.tg]) + Zhi[gzY.dz]);
            safe_copy(it.ganzhi[1], 8, std::string(Gan[gzM.tg]) + Zhi[gzM.dz]);
            safe_copy(it.ganzhi[2], 8, std::string(Gan[gzD.tg]) + Zhi[gzD.dz]);
            if (hourKnown)
                safe_copy(it.ganzhi[3], 8, std::string(Gan[hg]) + Zhi[hz]);
            count++;
        }
    }
    return count;
    }, 0);
}

// ═══════════════════════════════════════════════════════════
//  Eclipse / Astronomy API implementation
// ═══════════════════════════════════════════════════════════

static char* dup_string(const std::string &s) {
    char *p = (char*)std::malloc(s.size() + 1);
    if (p) std::memcpy(p, s.c_str(), s.size() + 1);
    return p;
}

char* sxwnl_search_eclipses(int year, int month, int count, bool is_solar) {
    return guard([&]() -> char * {
        std::string result = rs_search(year, month, count, !is_solar);
        return dup_string(result);
    }, nullptr);
}

char* sxwnl_calc_eclipse_detail(int year, int month, int day,
                                 int hour, int minute, double second,
                                 bool is_utc, double longitude, double latitude) {
    return guard([&]() -> char * {
        Time t(year, month, day, hour, minute, second);
        JINGWEI jw{longitude, latitude, "", ""};
        std::string result = rysCalc(t, is_utc, false, jw);
        return dup_string(result);
    }, nullptr);
}

char* sxwnl_calc_sun_moon_rise(int year, int month, int day,
                                double longitude, double latitude) {
    return guard([&]() -> char * {
        JINGWEI jw{longitude, latitude, "", ""};
        std::string result = shengjiang(year, month, day, jw);
        return dup_string(result);
    }, nullptr);
}

// 日月升降/中天 (单天版本, 与 JS SZJ.calcRTS(jd,1,...) 等价)
int sxwnl_calc_day_rts(int year, int month, int day,
                       double longitude, double latitude,
                       double tz_hours, SxwnlDayRTS *out) {
    return guard([&]() -> int {
        if (!out) return -1;
        std::memset(out, 0, sizeof(SxwnlDayRTS));

        // 默认值: 月相在跨日时可能空缺, 与 JS RTS1 一致
        safe_copy(out->moon_rise,     sizeof(out->moon_rise),     "--:--:--");
        safe_copy(out->moon_set,      sizeof(out->moon_set),      "--:--:--");
        safe_copy(out->moon_meridian, sizeof(out->moon_meridian), "--:--:--");

        // 当地正午对应的"伪 UT"JD (与 JS curJD 同源):
        //   t.h=12 视为 UTC, 实际是本地正午
        Time t{};
        t.Y = year; t.M = (uint8_t)month; t.D = day;
        t.h = 12; t.m = 0; t.s = 0.0;
        long double jd_local = JD::toJD(t) - J2000;       // 本地正午对应"伪 UT"JD
        long double sq       = -tz_hours / 24.0L;          // 与 JS 同号约定: 北京 sq = -8/24

        SunMoonRiseSet &szj = SunMoonRiseSet::getInstance();
        szj.setL((long double)longitude / radd);
        szj.setFa((long double)latitude / radd);

        // 太阳: 单天
        SJ rs = szj.St(jd_local + sq);
        safe_copy(out->sun_rise,     sizeof(out->sun_rise),     JD::timeStr(rs.s - sq));
        safe_copy(out->sun_set,      sizeof(out->sun_set),      JD::timeStr(rs.j - sq));
        safe_copy(out->sun_meridian, sizeof(out->sun_meridian), JD::timeStr(rs.z - sq));
        safe_copy(out->civil_dawn,   sizeof(out->civil_dawn),   JD::timeStr(rs.c - sq));
        safe_copy(out->civil_dusk,   sizeof(out->civil_dusk),   JD::timeStr(rs.h - sq));
        // JS: ch = h - c - 0.5  (timeStr +0.5, 所以这里 -0.5)
        safe_copy(out->light_length, sizeof(out->light_length), JD::timeStr(rs.h - rs.c - 0.5L));
        // JS: sj = j - s - 0.5
        safe_copy(out->day_length,   sizeof(out->day_length),   JD::timeStr(rs.j - rs.s - 0.5L));

        // 月亮: 扫描 i=-1..1 三天, 抓取真正落在当天的 升/中/降
        for (int i = -1; i <= 1; ++i) {
            SJ rm = szj.Mt(jd_local + i + sq);
            long long ofs;

            ofs = (long long)std::floor((double)(rm.s - sq + 0.5L)) - (long long)jd_local;
            if (ofs == 0) safe_copy(out->moon_rise,     sizeof(out->moon_rise),     JD::timeStr(rm.s - sq));

            ofs = (long long)std::floor((double)(rm.z - sq + 0.5L)) - (long long)jd_local;
            if (ofs == 0) safe_copy(out->moon_meridian, sizeof(out->moon_meridian), JD::timeStr(rm.z - sq));

            ofs = (long long)std::floor((double)(rm.j - sq + 0.5L)) - (long long)jd_local;
            if (ofs == 0) safe_copy(out->moon_set,      sizeof(out->moon_set),      JD::timeStr(rm.j - sq));
        }
        return 0;
    }, -1);
}

void sxwnl_string_free(char *str) {
    std::free(str);
}

// ═══════════════════════════════════════════════════════════════
//  Geo Directory API (薄包装 src/geo.cpp 的 GeoPostion 单例)
// ═══════════════════════════════════════════════════════════════
namespace {
inline void copy_city(const JINGWEI &jw, SxwnlGeoCity &out) {
    safe_copy(out.province, sizeof(out.province), jw.s);
    safe_copy(out.district, sizeof(out.district), jw.x);
    out.longitude = jw.J;
    out.latitude  = jw.W;
    out.timezone  = jw.tz;
}
} // namespace

int sxwnl_geo_list_provinces(SxwnlGeoProvince *out, int out_max) {
    return guard([&]() -> int {
        if (!out || out_max <= 0) return -1;
        auto &g = GeoPostion::getInstance();
        auto names = g.listProvinces();
        int n = 0;
        for (const auto &name : names) {
            if (n >= out_max) break;
            std::memset(&out[n], 0, sizeof(SxwnlGeoProvince));
            safe_copy(out[n].province, sizeof(out[n].province), name);
            out[n].city_count = (int)g.listCitiesIn(name).size();
            ++n;
        }
        return n;
    }, -1);
}

int sxwnl_geo_list_cities(const char *province, SxwnlGeoCity *out, int out_max) {
    return guard([&]() -> int {
        if (!province || !out || out_max <= 0) return -1;
        auto &g = GeoPostion::getInstance();
        auto cities = g.listCitiesIn(province);
        int n = 0;
        for (const auto &c : cities) {
            if (n >= out_max) break;
            std::memset(&out[n], 0, sizeof(SxwnlGeoCity));
            copy_city(c, out[n]);
            ++n;
        }
        return n;
    }, -1);
}

int sxwnl_geo_search(const char *keyword, SxwnlGeoCity *out, int out_max) {
    return guard([&]() -> int {
        if (!keyword || !out || out_max <= 0) return -1;
        auto &g = GeoPostion::getInstance();
        auto matches = g.search(keyword, out_max);
        int n = 0;
        for (const auto &c : matches) {
            if (n >= out_max) break;
            std::memset(&out[n], 0, sizeof(SxwnlGeoCity));
            copy_city(c, out[n]);
            ++n;
        }
        return n;
    }, -1);
}

int sxwnl_geo_list_timezone_groups(SxwnlTimezoneGroup *out, int out_max) {
    return guard([&]() -> int {
        if (!out || out_max <= 0) return -1;
        auto &g = GeoPostion::getInstance();
        const auto &groups = g.timezoneGroups();
        int n = 0;
        for (const auto &grp : groups) {
            if (n >= out_max) break;
            std::memset(&out[n], 0, sizeof(SxwnlTimezoneGroup));
            safe_copy(out[n].continent, sizeof(out[n].continent), grp.continent);
            safe_copy(out[n].country,   sizeof(out[n].country),   grp.country);
            out[n].timezone = grp.timezone;
            int cn = 0;
            for (const auto &city : grp.cities) {
                if (cn >= 8) break;
                safe_copy(out[n].cities[cn], sizeof(out[n].cities[cn]), city);
                ++cn;
            }
            out[n].city_count = cn;
            ++n;
        }
        return n;
    }, -1);
}

int sxwnl_geo_default(SxwnlGeoCity *out) {
    return guard([&]() -> int {
        if (!out) return -1;
        std::memset(out, 0, sizeof(SxwnlGeoCity));
        auto &g = GeoPostion::getInstance();
        JINGWEI jw = g.getDefaultGeoPos();
        copy_city(jw, *out);
        return 0;
    }, -1);
}

// ═══════════════════════════════════════════════════════════════
//  Festival API
// ═══════════════════════════════════════════════════════════════
int sxwnl_get_festivals(int year, int month, int day, SxwnlFestival *out) {
    return guard([&]() -> int {
        if (!out) return -1;
        std::unique_ptr<Day> d(sxtwl::fromSolar(year, (uint8_t)month, day));
        if (!d) return -1;
        std::memset(out, 0, sizeof(SxwnlFestival));
        festival::DayInfo info = d->getFestivalInfo();
        safe_copy(out->holiday, sizeof(out->holiday), info.holiday);
        safe_copy(out->major,   sizeof(out->major),   info.major);
        safe_copy(out->minor,   sizeof(out->minor),   info.minor);
        safe_copy(out->misc,    sizeof(out->misc),    info.misc);
        out->is_off_day = info.isOffDay;
        return 0;
    }, -1);
}

// ═══════════════════════════════════════════════════════════════
//  Year / Time string utilities
// ═══════════════════════════════════════════════════════════════
int32_t sxwnl_year_str_to_astro(const char *s) {
    if (!s) return INT32_MIN;
    return guard([&]() -> int32_t {
        return (int32_t)year_utils::year2Ayear(std::string(s));
    }, INT32_MIN);
}

int sxwnl_astro_year_to_str(int32_t aYear, bool full_style, char *out, int out_size) {
    return guard([&]() -> int {
        if (!out || out_size <= 0) return -1;
        std::string s = year_utils::Ayear2year((int)aYear, full_style);
        if ((int)s.size() + 1 > out_size) return -1;
        safe_copy(out, (size_t)out_size, s);
        return 0;
    }, -1);
}

double sxwnl_time_str_to_hour(const char *s) {
    if (!s) return std::nan("");
    return guard([&]() -> double {
        return year_utils::timeStr2hour(std::string(s));
    }, std::nan(""));
}

// ═══════════════════════════════════════════════════════════════
//  天象事件工具
//
//  内部约定:
//   - XL::* 接受/返回的 t 单位是"自 J2000 起的儒略世纪数"(jc = (jd-J2000)/36525)
//   - dt_T(d) 中 d 是"自 J2000 起的天数"
//   - 我们最终输出: jd = TT - dt_T (得到 UTC), 再 +8/24 转北京时
// ═══════════════════════════════════════════════════════════════

static double jc_to_bj_jd(long double jc) {
    long double d_tt = jc * 36525.0L;          // TT, days from J2000
    long double d_utc = d_tt - dt_T(d_tt);     // 转 UTC
    long double jd_bj = d_utc + J2000 + 8.0L / 24.0L;
    return (double)jd_bj;
}

static std::string fmt_time_from_jd(double jd_bj) {
    Time t = JD::JD2DD(jd_bj);
    char buf[40];
    std::snprintf(buf, sizeof(buf), "%04d-%02d-%02d %02d:%02d:%02d",
                  (int)t.Y, (int)t.M, (int)t.D, (int)t.h, (int)t.m, (int)t.s);
    return std::string(buf);
}

// 求 [year, year+1) 内年首/年尾的"自 J2000 起的天数"区间
static void year_window_d(int year, long double &d0, long double &d1) {
    Time a{}; a.Y = year;     a.M = 1; a.D = 1; a.h = 0; a.m = 0; a.s = 0;
    Time b{}; b.Y = year + 1; b.M = 1; b.D = 1; b.h = 0; b.m = 0; b.s = 0;
    d0 = JD::toJD(a) - J2000;
    d1 = JD::toJD(b) - J2000;
}

// ── 月相 ──
int sxwnl_get_moon_phases(int year, SxwnlMoonPhase *out, int max_count) {
    return guard([&]() -> int {
        if (!out || max_count <= 0) return 0;
        static const char *NAMES[4] = {"朔", "上弦", "望", "下弦"};

        long double d0, d1;
        year_window_d(year, d0, d1);

        // 估算起始相位序号 N: N * 7.38264 天 ≈ d0
        // 7.38264 = 朔望月 29.5306 / 4
        long double avg = 29.5306L / 4.0L;
        long long N = (long long)std::floor((d0 - 6) / avg) - 2;

        int count = 0;
        // 多估 4 期防止漏取
        for (long long k = N; count < max_count; ++k) {
            long double W = (long double)k * (3.141592653589793L / 2.0L);
            long double jc = XL::MS_aLon_t(W);    // TT, julian centuries
            long double d_tt = jc * 36525.0L;
            if (d_tt < d0 - 1) continue;
            if (d_tt >= d1)    break;

            SxwnlMoonPhase e{};
            e.phase = ((int)((k % 4) + 4)) % 4;
            safe_copy(e.name, sizeof(e.name), NAMES[e.phase]);
            e.jd = jc_to_bj_jd(jc);
            safe_copy(e.time, sizeof(e.time), fmt_time_from_jd(e.jd));
            out[count++] = e;
        }
        return count;
    }, 0);
}

// ── 月地近/远点 ──
int sxwnl_get_moon_apsides(int year, SxwnlMoonApsis *out, int max_count) {
    return guard([&]() -> int {
        if (!out || max_count <= 0) return 0;
        long double d0, d1;
        year_window_d(year, d0, d1);

        long double tc0 = d0 / 36525.0L;
        long double tc1 = d1 / 36525.0L;

        // 同时跟踪近点与远点, 周期约 27.55 天
        int count = 0;
        const long double STEP = 27.55454988L / 36525.0L;
        long double tc = tc0 - STEP * 2;
        while (tc < tc1 + STEP * 2 && count < max_count) {
            for (int kind = 0; kind < 2 && count < max_count; ++kind) {
                bool isMin = (kind == 0);
                Vector2 r = XL::moonMinR(tc, isMin);
                long double t_tt = r[0] * 36525.0L;
                if (t_tt < d0 || t_tt >= d1) continue;

                SxwnlMoonApsis e{};
                e.kind = kind;
                safe_copy(e.name, sizeof(e.name), kind == 0 ? "近地" : "远地");
                e.jd = jc_to_bj_jd(r[0]);
                safe_copy(e.time, sizeof(e.time), fmt_time_from_jd(e.jd));
                // r[1] 单位为 km (sxwnl 内部约定: XL1_calc returns 距离单位是 km)
                e.distance_km = (double)r[1];
                // 去重: 仅在新条目时序大于上次时记录
                if (count == 0 || e.jd > out[count - 1].jd + 0.5) {
                    out[count++] = e;
                } else if (e.jd > out[count - 1].jd) {
                    out[count - 1] = e;
                }
            }
            tc += STEP;
        }
        // 排序
        std::sort(out, out + count, [](const SxwnlMoonApsis &a, const SxwnlMoonApsis &b) {
            return a.jd < b.jd;
        });
        return count;
    }, 0);
}

// ── 月交点 ──
int sxwnl_get_moon_nodes(int year, SxwnlMoonNode *out, int max_count) {
    return guard([&]() -> int {
        if (!out || max_count <= 0) return 0;
        long double d0, d1;
        year_window_d(year, d0, d1);

        const long double STEP = 27.21222082L / 36525.0L;
        long double tc = d0 / 36525.0L - STEP * 2;
        long double tcEnd = d1 / 36525.0L + STEP * 2;

        int count = 0;
        while (tc < tcEnd && count < max_count) {
            for (int kind = 0; kind < 2 && count < max_count; ++kind) {
                bool asc = (kind == 0);
                Vector2 r = XL::moonNode(tc, asc);
                long double t_tt = r[0] * 36525.0L;
                if (t_tt < d0 || t_tt >= d1) continue;

                SxwnlMoonNode e{};
                e.kind = kind;
                safe_copy(e.name, sizeof(e.name), kind == 0 ? "升交" : "降交");
                e.jd = jc_to_bj_jd(r[0]);
                safe_copy(e.time, sizeof(e.time), fmt_time_from_jd(e.jd));
                if (count == 0 || e.jd > out[count - 1].jd + 0.5) {
                    out[count++] = e;
                } else if (e.jd > out[count - 1].jd) {
                    out[count - 1] = e;
                }
            }
            tc += STEP;
        }
        std::sort(out, out + count, [](const SxwnlMoonNode &a, const SxwnlMoonNode &b) {
            return a.jd < b.jd;
        });
        return count;
    }, 0);
}

// ── 行星天象 ──
static const char* PLANET_NAMES[8] = {
    "", "水星", "金星", "火星", "木星", "土星", "天王", "海王"
};

static const long double PLANET_SYNODIC[8] = {
    0, 115.8774, 583.9213, 779.9362, 398.8840, 378.0916, 369.6562, 367.4867
};

// 收集一类天象到 out
static void collect_event(int planet, int evt, long double jc_event,
                          long double d_lo, long double d_hi,
                          double value,
                          SxwnlPlanetEvent *out, int &count, int max_count) {
    if (count >= max_count) return;
    long double d_tt = jc_event * 36525.0L;
    if (d_tt < d_lo || d_tt >= d_hi) return;

    static const char* EVT_NAMES[7] = {
        "合", "冲", "东大距", "西大距", "顺留", "逆留", "合月"
    };
    // 内行星的"冲"实际上是"下合"
    const char* evtName = EVT_NAMES[evt];
    if (evt == 1 && (planet == 1 || planet == 2)) evtName = "下合";

    SxwnlPlanetEvent &e = out[count++];
    e.planet = planet;
    e.event  = evt;
    safe_copy(e.planet_name, sizeof(e.planet_name), PLANET_NAMES[planet]);
    safe_copy(e.event_name,  sizeof(e.event_name),  evtName);
    e.jd = jc_to_bj_jd(jc_event);
    safe_copy(e.time, sizeof(e.time), fmt_time_from_jd(e.jd));
    e.value = value;
}

int sxwnl_get_planet_events(int year, SxwnlPlanetEvent *out, int max_count) {
    return guard([&]() -> int {
        if (!out || max_count <= 0) return 0;
        long double d0, d1;
        year_window_d(year, d0, d1);

        int count = 0;
        for (int xt = 1; xt <= 7 && count < max_count; ++xt) {
            long double hh = PLANET_SYNODIC[xt] / 36525.0L;

            // 合 / 冲(下合): 沿会合周期遍历
            for (int f = 0; f < 2; ++f) {
                long double tc = d0 / 36525.0L - hh;
                long double tcEnd = d1 / 36525.0L + hh;
                long double lastT = -1e18L;
                while (tc < tcEnd && count < max_count) {
                    Vector2 r = XL::xingHR(xt, tc, f == 1);
                    if (std::abs((double)(r[0] - lastT)) > 1.0 / 36525.0) {
                        long double yangVal = r[1]; // 黄纬差(弧度)
                        collect_event(xt, f, r[0], d0, d1,
                                      (double)(yangVal * 180.0 / M_PI),
                                      out, count, max_count);
                        lastT = r[0];
                    }
                    tc += hh;
                }
            }

            // 大距: 仅 1,2
            if (xt == 1 || xt == 2) {
                for (int dx = 0; dx < 2; ++dx) {
                    long double tc = d0 / 36525.0L - hh;
                    long double tcEnd = d1 / 36525.0L + hh;
                    long double lastT = -1e18L;
                    while (tc < tcEnd && count < max_count) {
                        Vector2 r = XL::daJu(xt, tc, dx == 0);  // dx=0 -> 东
                        if (std::abs((double)(r[0] - lastT)) > 1.0 / 36525.0) {
                            collect_event(xt, dx == 0 ? 2 : 3, r[0], d0, d1,
                                          (double)(r[1] * 180.0 / M_PI),
                                          out, count, max_count);
                            lastT = r[0];
                        }
                        tc += hh;
                    }
                }
            }

            // 留: sn=1 顺留, sn=0 逆留
            for (int sn = 0; sn < 2; ++sn) {
                long double tc = d0 / 36525.0L - hh;
                long double tcEnd = d1 / 36525.0L + hh;
                long double lastT = -1e18L;
                while (tc < tcEnd && count < max_count) {
                    long double t = XL::xingLiu(xt, tc, sn == 0);
                    if (std::abs((double)(t - lastT)) > 1.0 / 36525.0) {
                        collect_event(xt, sn == 0 ? 4 : 5, t, d0, d1,
                                      0.0, out, count, max_count);
                        lastT = t;
                    }
                    tc += hh;
                }
            }

            // 合月: 约每月一次, 用朔望月步进
            {
                long double tc = d0 / 36525.0L - 29.5306L / 36525.0L;
                long double tcEnd = d1 / 36525.0L + 29.5306L / 36525.0L;
                long double lastT = -1e18L;
                while (tc < tcEnd && count < max_count) {
                    Vector2 r = XL::xingHY(xt, tc);
                    if (std::abs((double)(r[0] - lastT)) > 1.0 / 36525.0) {
                        collect_event(xt, 6, r[0], d0, d1,
                                      (double)(r[1] * 180.0 / M_PI),
                                      out, count, max_count);
                        lastT = r[0];
                    }
                    tc += 29.5306L / 36525.0L;
                }
            }
        }

        // 按 jd 排序
        std::sort(out, out + count,
                  [](const SxwnlPlanetEvent &a, const SxwnlPlanetEvent &b) {
                      return a.jd < b.jd;
                  });

        // 去重(同一事件可能因步进重复落入)
        int unique = 0;
        for (int i = 0; i < count; ++i) {
            if (unique > 0 &&
                out[unique - 1].planet == out[i].planet &&
                out[unique - 1].event  == out[i].event  &&
                std::abs(out[unique - 1].jd - out[i].jd) < 0.5) {
                continue;
            }
            out[unique++] = out[i];
        }
        return unique;
    }, 0);
}

// ═══════════════════════════════════════════════════════════════
//  Map projection & world map
// ═══════════════════════════════════════════════════════════════
SxwnlProjection sxwnl_projection_create(int type, double J0_rad, double W0_rad,
                                        double win_cx, double win_cy,
                                        double win_rx, double win_ry) {
    return guard([&]() -> SxwnlProjection {
        auto *p = new map_projection::Projector();
        p->setlx(type, (long double)J0_rad, (long double)W0_rad,
                 {(long double)win_cx, (long double)win_cy,
                  (long double)win_rx, (long double)win_ry});
        return (SxwnlProjection)p;
    }, (SxwnlProjection)nullptr);
}

void sxwnl_projection_free(SxwnlProjection proj) {
    if (!proj) return;
    delete reinterpret_cast<map_projection::Projector*>(proj);
}

int sxwnl_projection_point(SxwnlProjection proj,
                           double J_rad, double W_rad,
                           double *out_x, double *out_y, double *out_z) {
    return guard([&]() -> int {
        if (!proj || !out_x || !out_y || !out_z) return -1;
        auto *p = reinterpret_cast<map_projection::Projector*>(proj);
        auto pt = p->toxy((long double)J_rad, (long double)W_rad);
        *out_x = (double)pt.x; *out_y = (double)pt.y; *out_z = (double)pt.z;
        return 0;
    }, -1);
}

int sxwnl_projection_polyline(SxwnlProjection proj,
                              const double *input, int in_count,
                              double *out, int out_max) {
    return guard([&]() -> int {
        if (!proj || !input || in_count <= 0) return 0;
        auto *p = reinterpret_cast<map_projection::Projector*>(proj);
        std::vector<long double> in;
        in.reserve(in_count);
        for (int i = 0; i < in_count; ++i) in.push_back((long double)input[i]);
        auto result = p->lineArr(in);
        int sz = (int)result.size();
        if (out_max < sz) return -sz;
        if (out) for (int i = 0; i < sz; ++i) out[i] = (double)result[i];
        return sz;
    }, 0);
}

int sxwnl_world_map_get_ditu0(double *out, int out_max) {
    return guard([&]() -> int {
        const auto &v = world_map::ditu0();
        int sz = (int)v.size();
        if (out_max < sz) return -sz;
        if (out) for (int i = 0; i < sz; ++i) out[i] = (double)v[i];
        return sz;
    }, 0);
}

int sxwnl_world_map_set_data(int idx, const char *raw) {
    return guard([&]() -> int {
        if (!raw) return -1;
        return world_map::setMapData(idx, std::string(raw)) ? 0 : -1;
    }, -1);
}

int sxwnl_world_map_get_data(int idx, double *out, int out_max) {
    return guard([&]() -> int {
        const auto &v = world_map::getMapData(idx);
        int sz = (int)v.size();
        if (out_max < sz) return -sz;
        if (out) for (int i = 0; i < sz; ++i) out[i] = (double)v[i];
        return sz;
    }, 0);
}

// ═══════════════════════════════════════════════════════════════
//  星座 + 恒星库
// ═══════════════════════════════════════════════════════════════
int sxwnl_get_constellations(SxwnlConstellation *out, int max_count) {
    return guard([&]() -> int {
        if (!out || max_count <= 0) return 0;
        const auto &v = star_catalog::list88();
        int n = std::min<int>((int)v.size(), max_count);
        for (int i = 0; i < n; ++i) {
            SxwnlConstellation &o = out[i];
            safe_copy(o.name_abbr,   sizeof(o.name_abbr),   v[i].nameAbbr);
            o.area_sq = (double)v[i].areaSq;
            safe_copy(o.ra_str,      sizeof(o.ra_str),      v[i].raStr);
            safe_copy(o.dec_str,     sizeof(o.dec_str),     v[i].decStr);
            safe_copy(o.quad_family, sizeof(o.quad_family), v[i].quadFamily);
            safe_copy(o.name_en,     sizeof(o.name_en),     v[i].nameEn);
        }
        return n;
    }, 0);
}

static void fill_star(SxwnlStar &dst, const star_catalog::Star &src) {
    dst.ra_rad   = (double)src.ra;
    dst.dec_rad  = (double)src.dec;
    dst.pm_ra    = (double)src.pmRa;
    dst.pm_dec   = (double)src.pmDec;
    dst.parallax = (double)src.parallax;
    dst.mag      = (double)src.mag;
    safe_copy(dst.name, sizeof(dst.name), src.name);
    safe_copy(dst.info, sizeof(dst.info), src.info);
}

static star_catalog::Star to_star(const SxwnlStar &s) {
    star_catalog::Star o;
    o.ra       = (long double)s.ra_rad;
    o.dec      = (long double)s.dec_rad;
    o.pmRa     = (long double)s.pm_ra;
    o.pmDec    = (long double)s.pm_dec;
    o.parallax = (long double)s.parallax;
    o.mag      = (long double)s.mag;
    o.name     = s.name;
    o.info     = s.info;
    return o;
}

int sxwnl_register_star_library(const char *key, const char *raw) {
    return guard([&]() -> int {
        if (!key || !raw) return -1;
        star_catalog::registerLibrary(key, raw);
        return 0;
    }, -1);
}

int sxwnl_get_star_library(const char *key, bool include_all,
                           SxwnlStar *out, int max_count) {
    return guard([&]() -> int {
        if (!key || !out || max_count <= 0) return 0;
        auto v = star_catalog::getLibrary(key, include_all);
        int n = std::min<int>((int)v.size(), max_count);
        for (int i = 0; i < n; ++i) fill_star(out[i], v[i]);
        return n;
    }, 0);
}

int sxwnl_search_stars(const char *key, SxwnlStar *out, int max_count) {
    return guard([&]() -> int {
        if (!key || !out || max_count <= 0) return 0;
        auto v = star_catalog::search(key);
        int n = std::min<int>((int)v.size(), max_count);
        for (int i = 0; i < n; ++i) fill_star(out[i], v[i]);
        return n;
    }, 0);
}

int sxwnl_star_hx_calc(const SxwnlStar *in, int in_count,
                       double bj_jd, double nutation_q_days,
                       int mode, double longitude, double latitude,
                       SxwnlStarResult *out, int max_count) {
    return guard([&]() -> int {
        if (!in || !out || in_count <= 0 || max_count <= 0) return 0;
        // 北京时 jd → TT jd → 儒略世纪
        long double d_utc = (long double)bj_jd - J2000 - 8.0L / 24.0L;
        long double d_tt  = d_utc + dt_T(d_utc);
        long double t     = d_tt / 36525.0L;

        std::vector<star_catalog::Star> stars;
        stars.reserve(in_count);
        for (int i = 0; i < in_count; ++i) stars.push_back(to_star(in[i]));

        auto results = star_catalog::hxCalc(
            stars, t, (long double)nutation_q_days, mode,
            (long double)longitude * (long double)PI / 180.0L,
            (long double)latitude  * (long double)PI / 180.0L);

        int n = std::min<int>((int)results.size(), max_count);
        for (int i = 0; i < n; ++i) {
            safe_copy(out[i].name, sizeof(out[i].name), results[i].name);
            out[i].mag   = (double)results[i].mag;
            out[i].a_rad = (double)results[i].a;
            out[i].b_rad = (double)results[i].b;
        }
        return n;
    }, 0);
}

char* sxwnl_calc_planet_position(int planet, double bj_jd,
                                 double longitude, double latitude) {
    return guard([&]() -> char* {
        if (planet < 1 || planet > 7) return nullptr;
        // 转 TT: bj_jd 是北京时, JD::JD2DD - 8h -> UTC; UT -> TT 加 dt_T
        long double d_utc = (long double)bj_jd - J2000 - 8.0L / 24.0L;
        long double d_tt = d_utc + dt_T(d_utc);
        std::string s = XL::xingX(planet, d_tt, longitude * M_PI / 180.0,
                                  latitude * M_PI / 180.0);
        return dup_string(s);
    }, nullptr);
}

// ── 地球近/远日 ──
int sxwnl_get_earth_apsides(int year, SxwnlEarthApsis *out, int max_count) {
    return guard([&]() -> int {
        if (!out || max_count <= 0) return 0;
        long double d0, d1;
        year_window_d(year, d0, d1);

        // 一年内最多 1 近 1 远, 直接两次调用
        int count = 0;
        long double tc_mid = (d0 + d1) * 0.5L / 36525.0L;
        for (int kind = 0; kind < 2 && count < max_count; ++kind) {
            bool isMin = (kind == 0);
            Vector2 r = XL::earthMinR(tc_mid, isMin);
            long double t_tt = r[0] * 36525.0L;
            if (t_tt < d0 || t_tt >= d1) {
                // 跨年情况, 试相邻平点
                Vector2 r2 = XL::earthMinR(tc_mid - 1.0L, isMin);
                long double t2 = r2[0] * 36525.0L;
                if (t2 >= d0 && t2 < d1) r = r2;
                else continue;
            }
            SxwnlEarthApsis e{};
            e.kind = kind;
            safe_copy(e.name, sizeof(e.name), kind == 0 ? "近日" : "远日");
            e.jd = jc_to_bj_jd(r[0]);
            safe_copy(e.time, sizeof(e.time), fmt_time_from_jd(e.jd));
            e.distance_au = (double)r[1];
            out[count++] = e;
        }
        std::sort(out, out + count, [](const SxwnlEarthApsis &a, const SxwnlEarthApsis &b) {
            return a.jd < b.jd;
        });
        return count;
    }, 0);
}
