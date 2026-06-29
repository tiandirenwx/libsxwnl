#include "napi_bazi.h"
#include "napi_utils.h"
#include "sxwnl_capi.h"

// 将一柱(SxwnlBaziColumn)转换为 napi 对象
static napi_value columnToNapi(napi_env env, const SxwnlBaziColumn &c) {
    NArr cg(env, c.cang_gan_count);
    for (int i = 0; i < c.cang_gan_count; i++)
        cg.push(NObj(env).s("gan", c.cang_gan[i]).s("shiShen", c.cang_gan_shi_shen[i]));
    NArr ss(env, c.shen_sha_count);
    for (int i = 0; i < c.shen_sha_count; i++) ss.push(c.shen_sha[i]);
    return NObj(env)
        .s("gan", c.gan).s("zhi", c.zhi).s("ganShiShen", c.gan_shi_shen)
        .s("nayin", c.nayin).s("xingYun", c.xing_yun).s("ziZuo", c.zi_zuo)
        .s("kongWang", c.kong_wang)
        .v("cangGan", cg).v("shenSha", ss)
        .i("startYear", c.start_year);
}

napi_value NapiCalcBazi(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 1);

    SxwnlBaziParams p{};
    auto name = a.objStr(0, "name");
    strncpy(p.name, name.c_str(), sizeof(p.name) - 1);
    p.gender    = a.objInt(0, "gender") != 0; // 0男 1女
    p.is_ast    = a.objBool(0, "astEnabled");
    p.is_lunar  = a.objBool(0, "isLunar");
    p.is_leap   = a.objBool(0, "isLeapMonth");
    p.is_spec   = a.objBool(0, "isSpecMonth");
    p.year      = a.objInt(0, "year");
    p.month     = a.objInt(0, "month");
    p.day       = a.objInt(0, "day");
    p.hour      = a.objInt(0, "hour");
    p.minute    = a.objInt(0, "minute");
    p.second    = a.objDbl(0, "second");
    p.longitude = a.objDbl(0, "longitude");
    p.latitude  = a.objDbl(0, "latitude");
    p.lifa      = a.objInt(0, "lifa");

    SxwnlBazi bz = sxwnl_bazi_create(&p);
    if (!bz) return js_null(env);

    NObj result(env);
    result
        .s("userName",    sxwnl_bazi_get_user_name(bz))
        .s("gender",      sxwnl_bazi_get_gender(bz))
        .s("solarBirth",  sxwnl_bazi_get_solar_birth(bz))
        .s("lunarBirth",  sxwnl_bazi_get_lunar_birth(bz))
        .s("dateOfBirth", sxwnl_bazi_get_date_of_birth(bz))
        .s("shengXiao",   sxwnl_bazi_get_sheng_xiao(bz))
        .s("age",         sxwnl_bazi_get_age(bz))
        .s("lifa",        sxwnl_bazi_get_lifa(bz))
        .s("dingQiType",  sxwnl_bazi_get_ding_qi_type(bz))
        .s("ast",         sxwnl_bazi_get_ast(bz))
        .s("jieQi",       sxwnl_bazi_get_jie_qi(bz))
        .s("qiYun",       sxwnl_bazi_get_qi_yun(bz))
        .s("jiaoYun",     sxwnl_bazi_get_jiao_yun(bz));

    // ── 完整排盘列信息 ──
    SxwnlBaziColumn cols[4];
    if (sxwnl_bazi_get_columns(bz, cols) == 4) {
        NArr colArr(env, 4);
        for (int i = 0; i < 4; i++) colArr.push(columnToNapi(env, cols[i]));
        result.v("columns", colArr);
    }

    SxwnlBaziColumn curDy{};
    int curDyStart = 0;
    if (sxwnl_bazi_get_current_da_yun(bz, &curDy) == 0) {
        result.v("currentDaYun", columnToNapi(env, curDy));
        curDyStart = curDy.start_year;
    }
    SxwnlBaziColumn curLn{};
    if (sxwnl_bazi_get_current_liu_nian(bz, &curLn) == 0)
        result.v("currentLiuNian", columnToNapi(env, curLn));

    // 大运完整列 + 每个大运对应的 10 年流年/小运
    SxwnlBaziColumn dyCols[12];
    int dyN = sxwnl_bazi_get_da_yun_columns(bz, dyCols, 12);
    if (dyN > 0) {
        NArr a(env, dyN);
        NArr allLn(env, dyN);
        for (int i = 0; i < dyN; i++) {
            a.push(columnToNapi(env, dyCols[i]));
            SxwnlLiuNianItem lns[10];
            int ln = sxwnl_bazi_get_liu_nian(bz, dyCols[i].start_year, lns, 10);
            NArr bucket(env, ln > 0 ? ln : 0);
            for (int j = 0; j < ln; j++) {
                bucket.push(NObj(env)
                    .i("year", lns[j].year).i("age", lns[j].age)
                    .s("ganZhi", lns[j].gan_zhi)
                    .s("ganShiShen", lns[j].gan_shi_shen)
                    .s("zhiShiShen", lns[j].zhi_shi_shen)
                    .s("xiaoYun", lns[j].xiao_yun)
                    .s("xiaoYunShiShen", lns[j].xiao_yun_shi_shen));
            }
            allLn.push(bucket);
        }
        result.v("daYunColumns", a);
        result.v("liuNianAll", allLn);
    }

    // 五行统计 + 旺衰 + 司令
    int wx[5];
    sxwnl_bazi_get_wuxing_count(bz, wx, true);
    NArr wxArr(env, 5);
    for (int i = 0; i < 5; i++) wxArr.push(wx[i]);
    result.v("wuXingCount", wxArr);

    char st[5][8];
    sxwnl_bazi_get_wuxing_status(bz, st);
    NArr stArr(env, 5);
    for (int i = 0; i < 5; i++) stArr.push(st[i]);
    result.v("wuXingStatus", stArr);

    result.s("siLing", sxwnl_bazi_get_si_ling(bz));

    // 当前大运区间的流年 + 小运(10年)
    if (curDyStart > 0) {
        SxwnlLiuNianItem lns[10];
        int ln = sxwnl_bazi_get_liu_nian(bz, curDyStart, lns, 10);
        if (ln > 0) {
            NArr a(env, ln);
            for (int i = 0; i < ln; i++) {
                a.push(NObj(env)
                    .i("year", lns[i].year).i("age", lns[i].age)
                    .s("ganZhi", lns[i].gan_zhi)
                    .s("ganShiShen", lns[i].gan_shi_shen)
                    .s("zhiShiShen", lns[i].zhi_shi_shen)
                    .s("xiaoYun", lns[i].xiao_yun)
                    .s("xiaoYunShiShen", lns[i].xiao_yun_shi_shen));
            }
            result.v("liuNian", a);
        }
    }

    sxwnl_bazi_free(bz);
    return result;
}

// 四柱反推: 入参 (yg,yz,mg,mz,dg,dz,hg,hz, startYear, endYear)
// 时柱未知时 hg=hz=-1; 返回 {year,month,day,hour,yearGZ,monthGZ,dayGZ,hourGZ}[]
napi_value NapiBaziReverse(napi_env env, napi_callback_info info) {
    NArgs a(env, info, 10);
    int sz[8];
    for (int i = 0; i < 8; i++) sz[i] = a.intAt(i);
    int startY = a.intAt(8), endY = a.intAt(9);

    const int kMax = 60;
    SxwnlReverseItem out[kMax];
    int n = sxwnl_bazi_reverse(sz, startY, endY, out, kMax);

    NArr arr(env, n > 0 ? n : 0);
    for (int i = 0; i < n; i++) {
        arr.push(NObj(env)
            .i("year", out[i].year).i("month", out[i].month)
            .i("day", out[i].day).i("hour", out[i].hour)
            .s("yearGZ", out[i].ganzhi[0]).s("monthGZ", out[i].ganzhi[1])
            .s("dayGZ", out[i].ganzhi[2]).s("hourGZ", out[i].ganzhi[3]));
    }
    return arr;
}
