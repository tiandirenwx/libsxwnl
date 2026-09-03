// Quick build verification for sxwnl_capi
// Build: g++ -std=c++20 -I../src -I. test_capi.cpp sxwnl_capi.cpp ../src/*.cpp -o test_capi
#include "sxwnl_capi.h"
#include <cstdio>
#include <cstring>
#include <cstdlib>

int main() {
    printf("=== sxwnl_capi build test ===\n\n");

    // 1. Day info
    SxwnlDayInfo day;
    if (sxwnl_get_day_info(2026, 5, 20, &day) == 0) {
        printf("[Calendar] 2026-05-20\n");
        printf("  Lunar: %s%s  GanZhi: %s年 %s月 %s日\n",
               day.lunar_month_name, day.lunar_day_name,
               day.year_gz, day.month_gz, day.day_gz);
        printf("  ShengXiao: %s  Constellation: %s  %s\n",
               day.sheng_xiao, day.constellation_name, day.week_name);
        if (day.jie_qi >= 0)
            printf("  JieQi: %s @ %s\n", day.jie_qi_name, day.jie_qi_time);
    }

    // 2. Month data
    SxwnlDayInfo month[31];
    int n = sxwnl_get_month_days(2026, 5, month, 31);
    printf("\n[Month] 2026-05: %d days\n", n);

    // 3. JieQi list
    SxwnlJieQiItem jq[30];
    int jqn = sxwnl_get_jieqi_list(2026, jq, 30);
    printf("\n[JieQi] 2026: %d items\n", jqn);
    for (int i = 0; i < jqn; i++)
        printf("  %s: %d/%d %s\n", jq[i].name, jq[i].solar_month, jq[i].solar_day, jq[i].time);

    // 4. Bazi
    SxwnlBaziParams bp{};
    strcpy(bp.name, "Test");
    bp.gender = false;
    bp.year = 1990; bp.month = 6; bp.day = 15;
    bp.hour = 10; bp.minute = 30; bp.second = 0;
    bp.longitude = 116.38; bp.latitude = 39.9;

    SxwnlBazi bz = sxwnl_bazi_create(&bp);
    if (bz) {
        printf("\n[Bazi] %s\n", sxwnl_bazi_get_solar_birth(bz));
        printf("  Lunar: %s\n", sxwnl_bazi_get_lunar_birth(bz));
        printf("  ShengXiao: %s  Age: %s\n",
               sxwnl_bazi_get_sheng_xiao(bz), sxwnl_bazi_get_age(bz));

        SxwnlPillar pls[4];
        sxwnl_bazi_get_si_zhu(bz, pls);
        printf("  SiZhu:");
        for (int i = 0; i < 4; i++)
            printf(" %s%s", pls[i].tian_gan, pls[i].di_zhi);
        printf("\n");

        printf("  QiYun: %s  JiaoYun: %s\n",
               sxwnl_bazi_get_qi_yun(bz), sxwnl_bazi_get_jiao_yun(bz));

        sxwnl_bazi_free(bz);
    }

    // 5. Eclipse search
    char *ecl = sxwnl_search_eclipses(2026, 1, 3, true);
    if (ecl) {
        printf("\n[Eclipse] Solar eclipses:\n%.200s...\n", ecl);
        sxwnl_string_free(ecl);
    }

    // 6. Festivals
    printf("\n[Festival] 春节 2026-02-17:\n");
    SxwnlFestival ft{};
    if (sxwnl_get_festivals(2026, 2, 17, &ft) == 0) {
        printf("  holiday=%s\n  major=%s\n  minor=%s\n  misc=%s\n  off=%d\n",
               ft.holiday, ft.major, ft.minor, ft.misc, (int)ft.is_off_day);
    }
    printf("[Festival] 2026-10-01 (国庆):\n");
    if (sxwnl_get_festivals(2026, 10, 1, &ft) == 0) {
        printf("  holiday=%s off=%d\n", ft.holiday, (int)ft.is_off_day);
    }
    printf("[Festival] 2026-01-15 (数九):\n");
    if (sxwnl_get_festivals(2026, 1, 15, &ft) == 0) {
        printf("  misc=%s\n", ft.misc);
    }

    // 7. Year string utils
    printf("\n[YearUtils]\n");
    printf("  '公元前221' -> astro %d\n", sxwnl_year_str_to_astro("公元前221"));
    printf("  '-221'      -> astro %d\n", sxwnl_year_str_to_astro("-221"));
    printf("  '2024年'    -> astro %d\n", sxwnl_year_str_to_astro("2024年"));
    char yb[64];
    sxwnl_astro_year_to_str(-220, true,  yb, sizeof(yb)); printf("  astro -220 -> %s\n", yb);
    sxwnl_astro_year_to_str( 2024, true, yb, sizeof(yb)); printf("  astro 2024 -> %s\n", yb);
    printf("  '12:30:00' -> %.6f h\n", sxwnl_time_str_to_hour("12:30:00"));
    printf("  '8.5'      -> %.6f h\n", sxwnl_time_str_to_hour("8.5"));

    // 8. Moon phases
    printf("\n[MoonPhase] 2026 (first 5):\n");
    SxwnlMoonPhase mp[60];
    int mpn = sxwnl_get_moon_phases(2026, mp, 60);
    for (int i = 0; i < mpn && i < 5; ++i)
        printf("  %s @ %s\n", mp[i].name, mp[i].time);
    printf("  total = %d\n", mpn);

    // 9. Earth apsides
    printf("\n[EarthApsis] 2026:\n");
    SxwnlEarthApsis ea[4];
    int ean = sxwnl_get_earth_apsides(2026, ea, 4);
    for (int i = 0; i < ean; ++i)
        printf("  %s @ %s  r=%.6f AU\n", ea[i].name, ea[i].time, ea[i].distance_au);

    // 10. Moon nodes
    printf("\n[MoonNode] 2026 (first 4):\n");
    SxwnlMoonNode mn[60];
    int mnn = sxwnl_get_moon_nodes(2026, mn, 60);
    for (int i = 0; i < mnn && i < 4; ++i)
        printf("  %s @ %s\n", mn[i].name, mn[i].time);

    // 11. Planet events
    printf("\n[PlanetEvent] 2026 (first 12):\n");
    SxwnlPlanetEvent pe[200];
    int pen = sxwnl_get_planet_events(2026, pe, 200);
    for (int i = 0; i < pen && i < 12; ++i)
        printf("  %s %s @ %s  v=%.3f\n",
               pe[i].planet_name, pe[i].event_name, pe[i].time, pe[i].value);
    printf("  total = %d\n", pen);

    // 12. Map projection
    printf("\n[MapProjection]\n");
    SxwnlProjection mp4 = sxwnl_projection_create(4, 0, 0, 0, 0, 1, 1); // Mollweide
    double x, y, z;
    sxwnl_projection_point(mp4, 0, 0, &x, &y, &z);
    printf("  Mollweide(0, 0)        -> (%.4f, %.4f) z=%.1f\n", x, y, z);
    sxwnl_projection_point(mp4, 3.14159/2, 0.5, &x, &y, &z);
    printf("  Mollweide(90E, 28.6N)  -> (%.4f, %.4f) z=%.1f\n", x, y, z);
    sxwnl_projection_free(mp4);

    SxwnlProjection mp6 = sxwnl_projection_create(6, 0, 0, 0, 0, 1, 1); // Mercator
    sxwnl_projection_point(mp6, 2.04, 0.7, &x, &y, &z);  // 116E,40N
    printf("  Mercator(116E, 40N)    -> (%.4f, %.4f)\n", x, y);
    sxwnl_projection_free(mp6);

    // 13. World map ditu0
    printf("\n[WorldMap] ditu0\n");
    int ditu0_sz = sxwnl_world_map_get_ditu0(nullptr, 0);
    printf("  total elements = %d\n", -ditu0_sz);
    int n0 = -ditu0_sz;
    if (n0 > 0) {
        double *buf = (double*)malloc(sizeof(double) * n0);
        sxwnl_world_map_get_ditu0(buf, n0);
        // 计算有多少条折线
        int lines = 0;
        for (int i = 0; i < n0; ++i) if (buf[i] >= 1e6) lines++;
        printf("  polylines = %d\n", lines);
        printf("  first 4 entries: %.3f %.3f %.3f %.3f\n",
               buf[0], buf[1], buf[2], buf[3]);
        free(buf);
    }

    // 14. 88 constellations
    printf("\n[Constellation] (sample 5):\n");
    SxwnlConstellation xz[88];
    int nxz = sxwnl_get_constellations(xz, 88);
    printf("  total = %d\n", nxz);
    for (int i = 0; i < 5 && i < nxz; ++i)
        printf("  %s area=%.1f RA=%s DEC=%s [%s] %s\n",
               xz[i].name_abbr, xz[i].area_sq,
               xz[i].ra_str, xz[i].dec_str,
               xz[i].quad_family, xz[i].name_en);

    // 15. Star library + search
    printf("\n[Star] library 库0 first 3:\n");
    SxwnlStar st[64];
    int nst = sxwnl_get_star_library("库0", false, st, 64);
    for (int i = 0; i < nst && i < 3; ++i)
        printf("  %s | %s mag=%.2f RA=%.4f rad DEC=%.4f rad\n",
               st[i].name, st[i].info, st[i].mag, st[i].ra_rad, st[i].dec_rad);
    printf("  搜索 'Lyr': \n");
    nst = sxwnl_search_stars("Lyr", st, 64);
    for (int i = 0; i < nst; ++i)
        printf("    %s | %s mag=%.2f\n", st[i].name, st[i].info, st[i].mag);

    // 16. hxCalc 视位置 (J2000=2026.5)
    printf("\n[hxCalc] 织女星 视位置 @ 2026-06-29 21:00 北京:\n");
    nst = sxwnl_search_stars("α A0", st, 64);  // 织女星 Lyr α A0
    if (nst > 0) {
        // 北京时 2026-06-29 21:00 → JD = 2461220 + (21-12)/24 ≈
        // 简单写法: JD ~ 2461221.0417
        SxwnlStarResult sr[8];
        int nsr = sxwnl_star_hx_calc(st, 1, 2461221.0417, 0.0, 0,
                                     116.38, 39.9, sr, 8);
        for (int i = 0; i < nsr; ++i)
            printf("  %s mag=%.2f  视α=%.6f rad  视δ=%.6f rad\n",
                   sr[i].name, sr[i].mag, sr[i].a_rad, sr[i].b_rad);
    }

    // 17. 回历互转
    printf("\n[Moslem] 公历 ↔ 回历:\n");
    int hy, hm, hd;
    sxwnl_solar_to_moslem(2026, 6, 29, &hy, &hm, &hd);
    printf("  2026-06-29 → 回历 %d-%d-%d (期望 1448-1-13)\n", hy, hm, hd);
    int sy, sm, sd;
    sxwnl_moslem_to_solar(1448, 1, 13, &sy, &sm, &sd);
    printf("  回历 1448-01-13 → %d-%d-%d (期望 2026-6-29)\n", sy, sm, sd);
    sxwnl_moslem_to_solar(1, 1, 1, &sy, &sm, &sd);
    printf("  回历 1-01-01 (元年)  → %d-%d-%d (期望 622-7-16)\n", sy, sm, sd);

    SxwnlDayInfo dh;
    sxwnl_get_day_info_by_moslem(1448, 1, 13, &dh);
    printf("  回历查询: 公历=%d-%d-%d 农历=%s%s 星期=%d 干支日=%s\n",
           dh.solar_year, dh.solar_month, dh.solar_day,
           dh.lunar_month_name, dh.lunar_day_name, dh.week_day, dh.day_gz);
    printf("  同日回历(双向一致): %d-%d-%d\n",
           dh.moslem_year, dh.moslem_month, dh.moslem_day);

    // 18. 综合月历(三历 + 干支 + 节气 + 月相 + 节日)
    printf("\n[CalendarMonth] 2026-06 综合月历:\n");
    SxwnlCalendarMonthHeader head;
    SxwnlDayInfo cm[31];
    int dn = sxwnl_get_calendar_month(2026, 6, &head, cm, 31);
    printf("  %d 年 %d 月 首周=%d 总周=%d 干支=%s(%s年) 年号=%s\n",
           head.year, head.month, head.first_week_day, head.total_weeks,
           head.year_gz, head.sheng_xiao, head.nianhao);
    printf("  天数=%d 月首JD=%d\n", dn, head.first_julian_day);
    for (int i = 0; i < dn; ++i) {
        if (i % 7 != 0 && cm[i].jie_qi < 0 && cm[i].yue_xiang < 0 &&
            !cm[i].holiday[0] && !cm[i].misc[0]) continue;
        printf("  %02d周%s | 农:%s%s 回:%d-%d-%d 干支:%s%s%s 节气:%s%s 月相:%s%s 节日:%s%s%s%s\n",
               cm[i].solar_day, cm[i].week_name + 6,
               cm[i].lunar_month_name, cm[i].lunar_day_name,
               cm[i].moslem_year, cm[i].moslem_month, cm[i].moslem_day,
               cm[i].year_gz, cm[i].month_gz, cm[i].day_gz,
               cm[i].jie_qi_name, cm[i].jie_qi_time,
               cm[i].yue_xiang_name, "",
               cm[i].holiday,
               cm[i].holiday[0] && cm[i].misc[0]?"/":"",
               cm[i].misc,
               cm[i].is_off_day?"(休)":"");
    }

    printf("\n=== All tests passed ===\n");
    return 0;
}
