#include <exception>
#include <csignal>
#include <execinfo.h>
#include <unistd.h>

#include "sx_qs_test.h"
#include "sx_date_test.h"
#include "sx_xingli_test.h"
#include "sx_tianxiang_test.h"
#include "geo.h"
#include "eph_show.h"
#include "bazi.h"
#include "sxtwl.h"
#include "day.h"
#include "festival_zh.h"
#include "year_utils.h"
#include "map_projection.h"
#include "world_map.h"
#include "star_eph.h"
#include "star_catalog.h"
#include "month_calendar.h"
#include "JD.h"
#include "sx_lang_zh.h"
#include <iostream>
#include <iomanip>

//===============================================================
#ifdef _WIN32
#include <windows.h>
#endif

//===============================================================
void terminate_handler()
{
    std::cerr << "Unhandled exception caught in terminate handler." << std::endl;

    // 获取调用栈信息
    void *array[20];
    int size = backtrace(array, 20);

    // 输出调用栈信息
    std::cerr << "Backtrace:" << std::endl;
    backtrace_symbols_fd(array, size, STDERR_FILENO);

    // 杀死进程并生成core文件
    signal(SIGABRT, SIG_DFL);
    std::abort();
}

int main(int argc, char *argv[])
{
    // 设置异常处理函数
    std::set_terminate(terminate_handler);
#ifdef _WIN32
    // windows中在控制台中输出包含 Unicode 字符的文本
    SetConsoleOutputCP(65001);
    // system("@chcp 65001");
#endif

    // 综合月历: cal [year] [month]
    if (argc >= 2 && "cal" == std::string(argv[1]))
    {
        int y = (argc > 2) ? std::atoi(argv[2]) : 2026;
        int m = (argc > 3) ? std::atoi(argv[3]) : 6;
        sxwnl::MonthCalendar cal(y, m);
        std::cout << "════ " << cal.year() << "年" << cal.month() << "月  "
                  << Gan[cal.yearGZ().tg] << Zhi[cal.yearGZ().dz]
                  << "(" << cal.shengXiao() << "年)  "
                  << "月首星期:" << (int)cal.firstWeek()
                  << "  天数:" << cal.dayCount()
                  << "  总周:" << (int)cal.totalWeeks() << std::endl;
        std::cout << "  年号: " << cal.nianHao() << std::endl;
        for (const auto &cd : cal.days())
        {
            Day &d = *cd.day;
            GZ ygz = d.getYearGZ(true);
            GZ mgz = d.getMonthGZ();
            GZ dgz = d.getDayGZ();
            std::cout << "  " << std::setw(2) << d.getSolarDay()
                      << " 周" << (int)d.getWeek()
                      << "  农:" << d.getLunarMonthName() << d.getLunarDayName()
                      << "  回:" << d.getMoslemYear() << "-"
                      << (int)d.getMoslemMonth() << "-" << d.getMoslemDay()
                      << "  干支:" << Gan[ygz.tg] << Zhi[ygz.dz]
                                   << Gan[mgz.tg] << Zhi[mgz.dz]
                                   << Gan[dgz.tg] << Zhi[dgz.dz];
            if (d.hasJieQi())    std::cout << "  ✦" << d.getJieQiName()    << " " << d.getJieQiTimeStr();
            if (d.hasYueXiang()) std::cout << "  ●" << d.getYueXiangName() << " " << d.getYueXiangTimeStr();
            auto info = d.getFestivalInfo();
            if (!info.holiday.empty()) std::cout << "  ★" << info.holiday;
            if (!info.misc.empty())    std::cout << "  "  << info.misc;
            if (info.isOffDay)         std::cout << " [休]";
            std::cout << std::endl;
        }
        return 0;
    }

    if (2 == argc)
    {
        // for full test
        if ("full" == std::string(argv[1]))
        {
            lunar2solarFullTest();
            return 0;
        }
        else if ("qs" == std::string(argv[1]))
        {
            dingQi_v();
            dingSuo_v();
            return 0;
        }
        else if ("xl" == std::string(argv[1]))
        {
            // 星历计算
            Time t{};
            t.h = 0, t.m = 0, t.s = 0;
            t.Y = 2008;
            t.M = 1;
            t.D = 1;
            for (size_t i = 1; i < 11; i++)
            {
                pCalc(i, t, 3, 1, true);
            }
            return 0;
        }
        else if ("tx" == std::string(argv[1]))
        {
            Time t{};
            t.h = 0, t.m = 0, t.s = 0;
            t.Y = 2008;
            t.M = 1;
            t.D = 1;
            for (size_t i = 1; i < 16; i++)
            {
                tianXiang(i, 2, t, 10);
            }
            return 0;
        }
        else if ("rs" == std::string(argv[1]))
        {
            // 2008-08-01 10:18:21 TD
            //  日月食
            GeoPostion &gep = GeoPostion::getInstance();
            JINGWEI jw = gep.getCityGeoPos();
            Time t = {2008, 8, 1, 18, 17, 15.0};
            std::cout << rysCalc(t, true, false,jw) << std::endl;
            std::cout << rs_search(2008, 8, 200, true) << std::endl; // 日食粗搜索
            std::cout << rs2_calc(5, 0, 29.5306) << std::endl;
            std::cout << rs2_jxb() << std::endl;
            return 0;
        }
        else if ("sj" == std::string(argv[1]))
        {
            Time t = JD::getNowTime();
            GeoPostion &gep = GeoPostion::getInstance();
            JINGWEI jw = gep.getCityGeoPos();
            std::cout << shengjiang(t.Y, t.M, t.D, jw) << std::endl;
            std::cout << "---------太阳升降----------" << std::endl;
            std::cout << shengjiang2(t.Y, jw) << std::endl;
            std::cout << "----------时差表-----------" << std::endl;
            std::cout << shengjiang3(t.Y) << std::endl;
            return 0;
        }
        else if ("ty" == std::string(argv[1]))
        {
            // 真太阳时
            auto nowTime = JD::getNowTime();
            // auto nowTime = Time{1984,2,10,7,35,10};
            auto pty = JD::timeStr(nowTime);
            GeoPostion &gep = GeoPostion::getInstance();
            JINGWEI jw = gep.getCityGeoPos();
            auto astTime = JD::calcAST(nowTime, jw.J);
            auto zty = std::get<2>(astTime);
            std::cout << jw.s << jw.x << "平太阳时:" << pty << ",真太阳时:" << zty << std::endl;
            return 0;
        }
        else if ("ft" == std::string(argv[1]))
        {
            // 节日演示: 输出今日 + 几个典型节日
            auto demo = [](int y, int m, int d, const char *label) {
                std::unique_ptr<Day> day(sxtwl::fromSolar(y, m, d));
                festival::DayInfo f = day->getFestivalInfo();
                std::cout << label << " " << y << "-" << m << "-" << d
                          << " holiday=" << f.holiday
                          << " major="   << f.major
                          << " minor="   << f.minor
                          << " misc="    << f.misc
                          << " off="     << f.isOffDay << std::endl;
            };
            demo(2026, 1, 1,  "[元旦]");
            demo(2026, 2, 17, "[春节]");
            demo(2026, 4, 5,  "[清明]");
            demo(2026, 10, 1, "[国庆]");
            demo(2026, 12, 22,"[冬至]");
            return 0;
        }
        else if ("yu" == std::string(argv[1]))
        {
            // 纪年/时间字符串工具演示
            std::cout << "公元前 221 -> "  << year_utils::year2Ayear("公元前 221") << std::endl;
            std::cout << "-221       -> "  << year_utils::year2Ayear("-221") << std::endl;
            std::cout << "2024年     -> "  << year_utils::year2Ayear("2024年") << std::endl;
            std::cout << "astro -220 -> "  << year_utils::Ayear2year(-220) << std::endl;
            std::cout << "astro 2024 -> "  << year_utils::Ayear2year(2024) << std::endl;
            std::cout << "12:30:00   -> "  << year_utils::timeStr2hour("12:30:00") << " h" << std::endl;
            std::cout << "8:15:30.5  -> "  << year_utils::timeStr2hour("8:15:30.5") << " h" << std::endl;
            return 0;
        }
        else if ("mp" == std::string(argv[1]))
        {
            // 9 种地图投影演示
            map_projection::Projector p;
            const char *names[9] = {
                "0平面正投","1斜轴等距方位","2斜轴等积方位","3斜轴等角方位",
                "4摩尔威特","5正轴等距圆柱","6正轴等角圆柱","7多圆锥","8中国灯笼"
            };
            for (int t = 0; t <= 8; ++t)
            {
                p.setlx(t, 0.0L, 0.0L);
                // (北京 116E,40N)
                auto pt = p.toxy(2.0245L, 0.6981L);
                std::cout << names[t]
                          << " 北京 -> (" << std::fixed << std::setprecision(6)
                          << (double)pt.x << ", " << (double)pt.y << ") z="
                          << (double)pt.z << std::endl;
            }
            const auto &v = world_map::ditu0();
            std::cout << "ditu0 points = " << v.size() / 2 << std::endl;
            return 0;
        }
        else if ("hx" == std::string(argv[1]))
        {
            using namespace star_catalog;
            const auto &xz = list88();
            std::cout << "88 星座共 " << xz.size() << " 个" << std::endl;
            std::cout << "[0] " << xz[0].nameAbbr << "  英文: " << xz[0].nameEn
                      << "  RA=" << xz[0].raStr << "  DEC=" << xz[0].decStr
                      << std::endl;

            auto stars = getLibrary("库1");
            std::cout << "库1 主星数 = " << stars.size() << std::endl;
            for (const auto &s : stars)
                std::cout << "  " << s.name << " | " << s.info
                          << "  ra=" << (double)s.ra
                          << " dec=" << (double)s.dec
                          << " mag=" << (double)s.mag << std::endl;

            // 织女星 视位置
            auto vega = search(std::string("α A0"));
            if (!vega.empty())
            {
                long double jd_bj = 2461221.0417L;          // 2026-06-29 21:00 北京
                long double d_utc = jd_bj - J2000 - 8.0L / 24.0L;
                long double t     = (d_utc + dt_T(d_utc)) / 36525.0L;
                auto r = hxCalc(vega, t, 0.0L, 0,
                                116.38L * (long double)PI / 180.0L,
                                39.9L  * (long double)PI / 180.0L);
                for (const auto &x : r)
                    std::cout << "视位置 " << x.name
                              << "  α=" << (double)x.a
                              << "  δ=" << (double)x.b << std::endl;
            }
            return 0;
        }
        else if ("bz" == std::string(argv[1]))
        {
            // 真太阳时
            auto nowTime = JD::getNowTime();
            GeoPostion &gep = GeoPostion::getInstance();
            JINGWEI jw = gep.getCityGeoPos();
            SBaziInputPara sBZ;
            sBZ.birthDayTime = nowTime;
            sBZ.calendar = CalendarSolar;
            sBZ.gender = false;
            sBZ.isAst = true;
            sBZ.isRun = false;
            sBZ.isSpec = false;
            sBZ.jw = jw;
            sBZ.lifa = YuWuWeiZiPingLifa_DingDongZhi;
            sBZ.name = "无名";
            BaziBase obj(sBZ);
            obj.calcBaziPaiPan();
            auto s = obj.printBazi();
            std::cout << s << std::endl;
            return 0;
        }
        else
        {
            lunarYearTest(std::stoi(argv[1]), 1);
        }

        return 0;
    }

    if (3 == argc && "jqlist" == std::string(argv[1]))
    {
        // 测试getJieQiList函数：通过输入年份计算全年节气列表
        auto year = std::stoi(argv[2]);
        getJieQiList(year);
        return 0;
    }

    if (3 == argc)
    {
        lunarYearTest(std::stoi(argv[1]), std::stoi(argv[2]));
        return 0;
    }

    if (4 == argc)
    {
        if("jq" == std::string(argv[1]))
        {
            auto year = std::stoi(argv[2]),month = std::stoi(argv[3]);
            jqCalc(year,month);
        }
        else
        {
            auto y = std::stoi(argv[1]), m = std::stoi(argv[2]), d = std::stoi(argv[3]);
            lunar2solarSingleTest(y, m, d);
        }

        return 0;
    }

    if (argc == 8 && "bz" == std::string(argv[1]))
    {
        Time birthTime = {std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]), std::stod(argv[6]), std::stod(argv[7]), 0};
        // GeoPostion &gep = GeoPostion::getInstance();
        // JINGWEI jw = gep.getCityGeoPos();
        JINGWEI jw = {120, 39.9, "默认", "北京"};
        SBaziInputPara sBZ;
        sBZ.birthDayTime = birthTime;
        sBZ.calendar = CalendarLunar;
        sBZ.gender = false;
        if (std::stoi(argv[2]) == 1)
        {
            sBZ.gender = true;
        }
        sBZ.isAst = false;
        sBZ.isRun = false;
        sBZ.isSpec = false;
        sBZ.jw = jw;
        sBZ.lifa = XianDaiNongLifa_DingQiFa;
        sBZ.name = "无名";
        BaziBase obj(sBZ);
        obj.calcBaziPaiPan();
        auto s = obj.printBazi();
        std::cout << s << std::endl;
        return 0;
    }

    if (argc == 10 && "bz" == std::string(argv[1]))
    {
        Time birthTime = {std::stoi(argv[5]), std::stoi(argv[6]), std::stoi(argv[7]), std::stod(argv[8]), std::stod(argv[9]), 0};
        // GeoPostion &gep = GeoPostion::getInstance();
        // JINGWEI jw = gep.getCityGeoPos();
        JINGWEI jw={120, 39.9, "默认", "北京"};
        SBaziInputPara sBZ;
        sBZ.birthDayTime = birthTime;
        sBZ.calendar = CalendarLunar;
        sBZ.gender = (std::stoi(argv[2]) == 1);
        sBZ.isAst = false;
        sBZ.isRun = (std::stoi(argv[3]) == 1);
        sBZ.isSpec = (std::stoi(argv[4]) == 1);
        sBZ.jw = jw;
        sBZ.lifa = YuWuWeiZiPingLifa_DingDongZhi;
        sBZ.name = "无名";
        BaziBase obj(sBZ);
        obj.calcBaziPaiPan();
        auto s = obj.printBazi();
        std::cout << s << std::endl;
        return 0;
    }

    // lunar2solarSingleTest(-721, 1, 1);
    /*
    //  星历计算
    Time t{};
    t.h = 0, t.m = 0, t.s = 0;
    t.Y = 2008;
    t.M = 1;
    t.D = 1;
    pCalc(1, t, 10, 1, true);


    Time t = {2008, 8, 1, 18, 17, 15.0};
    std::cout << rysCalc(t, true, false) << std::endl;
    std::cout << rs_search(2008, 8, 200, 1) << std::endl; // 日食粗搜索
    std::cout << rs2_calc(5, 0, 29.5306) << std::endl;
    std::cout << rs2_jxb() << std::endl;
   */

      
    
    // 命理八字
    /*
    auto nowTime = JD::getNowTime();
    ML_calc(nowTime);
    return 0;
    */
    
    Time birthTime = {1978, 11, 8, 14, 32, 0};
    //GeoPostion &gep = GeoPostion::getInstance();
    //JINGWEI jw = gep.getCityGeoPos();
    JINGWEI jw={120, 39.9, "深圳", "蛇常堂"};
    SBaziInputPara sBZ;
    sBZ.birthDayTime = birthTime;
    sBZ.calendar = CalendarSolar;
    sBZ.gender = false;
    sBZ.isAst = false;
    sBZ.isRun = false;
    sBZ.isSpec = false;
    sBZ.jw = jw;
    sBZ.lifa = YuWuWeiZiPingLifa_DingXiaZhi;
    sBZ.name = "禚主";
    BaziBase obj(sBZ);
    obj.calcBaziPaiPan();
    auto s = obj.printBazi();
    std::cout << s << std::endl;

    return 0;

    // // 定气节令获取

    // for (int i = 1900; i < 2300; i++)
    // {
    //     for (int j = 1; j <= 12; j++)
    //     {
    //         jqCalc(i, j);
    //
    //

    return 0;
}
