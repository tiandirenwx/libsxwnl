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
        sBZ.calendar = CalendarSolar;
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
        sBZ.lifa = YuWuWeiZiPingLifa_DingXiaZhi;
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
    /*
    Time birthTime = {2017, 9, 2, 12, 37, 0};
    //GeoPostion &gep = GeoPostion::getInstance();
    //JINGWEI jw = gep.getCityGeoPos();
    JINGWEI jw={120, 39.9, "默认", "默认"};
    SBaziInputPara sBZ;
    sBZ.birthDayTime = birthTime;
    sBZ.calendar = CalendarLunar;
    sBZ.gender = true;
    sBZ.isAst = false;
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
    */

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
