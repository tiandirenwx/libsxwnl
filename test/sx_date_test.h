#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <chrono>
#include "sxtwl.h"
#include "sx_lang_zh.h"
#include "lunar_ob.h"
#include "geo.h"

int round_double(double number)
{
    return (number > 0.0) ? (number + 0.5) : (number - 0.5);
}

GZ getGZ(std::string tgStr, std::string dzStr)
{
    int tg = -1;
    int dz = -1;
    for (size_t i = 0; i < 10; i++)
    {
        if (std::string(Gan[i]) == tgStr)
        {
            tg = i;
            break;
        }
    }

    for (size_t i = 0; i < 12; i++)
    {
        if (std::string(Zhi[i]) == dzStr)
        {
            dz = i;
            break;
        }
    }
    return GZ(tg, dz);
}

std::string printLunarMonth(Day &day)
{
    std::string strTemp1, strTemp2;
    int year = day.getSolarYear();
    bool bdLeapYear = false;
    if (day.isLunarLeap())
    {
        if (year >= -721 && year < -220)
        {
            strTemp1 = BDLeapYueName[0];
            bdLeapYear = true;
        }
        else if (year >= -220 && year <= -104)
        {
            strTemp1 = BDLeapYueName[1];
            bdLeapYear = true;
        }
        else
        {
            strTemp1 = "闰";
        }
    }
    if (bdLeapYear)
    {
        return strTemp1;
    }

    if (day.isSpecNextMonth())
    {
        strTemp2 = SYmc[day.getLunarMonth() - 1];
    }
    else
    {
        strTemp2 = Ymc[day.getLunarMonth() - 1];
    }

    return strTemp1 + strTemp2 + "月";
}

void printDay(Day &day)
{

    std::cout << "\n===================================================" << std::endl;

    // 公历
    int solarYear = day.getSolarYear();
    std::cout << "公历：" << solarYear << "年"
              << (int)day.getSolarMonth() << "月"
              << day.getSolarDay() << "日" << std::endl;

    // 农历
    std::cout << "农历：" << day.getLunarYear() << "年"
              << printLunarMonth(day)
              << Rmc[day.getLunarDay() - 1] << "日"
              << ",黄帝纪年:" << day.getHuangdiYear() << "年"
              << ",儒略日:" << std::fixed << std::setprecision(6) << day.getJulianDate()
              << ",执日:" << RiJian12[day.getLunarDay12Jian()]
              << std::endl;

    // 回历
    std::cout << "回历" << day.getMoslemYear() << "年"
              << (int)day.getMoslemMonth() << "月"
              << day.getMoslemDay() << "日" << std::endl;

    // 年天二地支
    std::cout << "干支纪年:"
              << Gan[day.getYearGZ().tg] << Zhi[day.getYearGZ().dz] << "年"
              << Gan[day.getMonthGZ().tg] << Zhi[day.getMonthGZ().dz] << "月"
              << Gan[day.getDayGZ().tg] << Zhi[day.getDayGZ().dz] << "日"
              << std::endl;

    // 星期几
    std::cout << "星期:" << WeekCn[day.getWeek()]
              << ",本月第一天星期:" << WeekCn[day.getFirstWeekDayOfMonth()]
              << ",周序数:" << NumCn[(int)day.getWeekIndex()]
              << ",本月总周数:" << NumCn[(int)day.getTotalWeekNumsOfMonth()]
              << std::endl;
    std::cout << "星座:" << XiZ[day.getConstellation()] << std::endl;

    // 节气
    if (day.hasJieQi())
    {
        auto jd = day.getJieQiJD();
        auto t = sxtwl::JD2DD(jd);

        jd = sxtwl::toJD(t);
        std::cout << Jqmc[day.getJieQi()] << ":公历" << t.getYear() << "-" << t.getMonth() << "-"
                  << t.getDay() << " " << int(t.getHour()) << ":" << int(t.getMin()) << ":" << round_double(t.getSec())
                  << ",儒略日:" << std::fixed << std::setprecision(6) << jd << std::endl;
    }

    // 月相
    if (day.hasYueXiang())
    {
        auto jd = day.getYueXiangJD();
        auto t = sxtwl::JD2DD(jd);

        jd = sxtwl::toJD(t);
        std::cout << YueXiangName[day.getYueXiang()] << ":公历" << t.getYear() << "-" << t.getMonth() << "-"
                  << t.getDay() << " " << int(t.getHour()) << ":" << int(t.getMin()) << ":" << round_double(t.getSec())
                  << ",儒略日:" << std::fixed << std::setprecision(6) << jd << std::endl;
    }

    int lunMonth = sxtwl::getRunMonth(solarYear);
    if (lunMonth > 0)
    {
        std::cout << solarYear << "年"
                  << "闰" << lunMonth << "月" << std::endl;
    }
}

bool checkDayConvert(Day &day, std::string &s,bool isFullTest = false)
{
    if (!isFullTest)
    {
        printDay(day);
    }

    auto sYear = day.getSolarYear();
    auto sMonth = day.getSolarMonth();
    auto sDay = day.getSolarDay();
    auto lYear = day.getLunarYear();
    auto lMonth = day.getLunarMonth();
    auto lDay = day.getLunarDay();
    auto isRun = day.isLunarLeap();
    auto isSpec = day.isSpecNextMonth();

    auto lunarDay = Day::fromLunar(lYear, lMonth, lDay, isRun, isSpec);
    if (nullptr == lunarDay)
    {
        s = "程序异常1:输入公历: " + std::to_string(sYear) + "年" + std::to_string(sMonth) + "月" + std::to_string(sDay) + "日";
        return false;
    }
    if (lunarDay->getSolarYear() != sYear || lunarDay->getSolarMonth() != sMonth || lunarDay->getSolarDay() != sDay)
    {
        s = "异常例子:公历" + std::to_string(sYear) + "年" + std::to_string(sMonth) + "月" + std::to_string(sDay) + "日, " +
            "转换后农历:" + std::to_string(lunarDay->getLunarYear()) + "年" + std::to_string(lunarDay->getLunarMonth()) + "月" +
            std::to_string(lunarDay->getLunarDay()) + "日, " +
            "农历再转公历:" + std::to_string(lunarDay->getSolarYear()) + "年" +
            std::to_string(lunarDay->getSolarMonth()) + "月" + std::to_string(lunarDay->getSolarDay()) + "日";
        delete lunarDay;
        return false;
    }

    delete lunarDay;

    return true;
}

void lunar2solarFullTest()
{
    auto start = std::chrono::high_resolution_clock::now();
    // for full test
    Day *dayStart = sxtwl::fromSolar(-721, 1, 1);
    int i = 0, j = 0;
    std::vector<std::string> vecErrStr;
    for (;;)
    {
        std::string strTmp;

        if (checkDayConvert(*dayStart, strTmp, true))
        {
            Day *prevDay = dayStart;
            dayStart = dayStart->after(1);
            delete prevDay;
            j++;
        }
        else
        {
            vecErrStr.push_back(strTmp);
            i++;
            break;
        }

        if (dayStart->getSolarYear() > 10000)
        {
            break;
        }
    }
    delete dayStart;

    for (const auto &it : vecErrStr)
    {
        std::cout << it << std::endl;
    }

    std::cout << "程序共有:" << i << "例子未通过" << std::endl;
    std::cout << "程序共有:" << j << "例子通过" << std::endl;

    auto end = std::chrono::high_resolution_clock::now();
    // 计算函数的耗时
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "耗时：" << duration.count() << " 微秒" << std::endl;
}

// for one case
void lunar2solarSingleTest(int year, uint8_t month, int day)
{
    auto start = std::chrono::high_resolution_clock::now();

    Day *dayTmp = sxtwl::fromSolar(year, month, day);
    std::string errStr;
    auto ret = checkDayConvert(*dayTmp, errStr, false);
    if (!ret)
    {
        std::cout << errStr << std::endl;
    }
    delete dayTmp;
    auto end = std::chrono::high_resolution_clock::now();

    // 计算函数的耗时
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << "耗时：" << duration.count() << " 微秒" << std::endl;
}

void lunarYearTest(int year, int n)
{
    auto y = year;
    for (int i = 0; i < n; i++)
    {
        LunarYear ob(y);
        auto nh = ob.getNianHao();
        std::cout << y << "的年号是:" << nh << std::endl;
        auto nianli = ob.getNianLiStr();
        std::cout << y << "的年历是:"
                  << "\n"
                  << nianli << std::endl;
        y++;
    }
}

struct MLBZ
{
    std::string bz_jn;
    std::string bz_jy;
    std::string bz_jr;
    std::string bz_js;
    std::string bz_JS;
    std::string bz_zty;
};

void mingLiBaZi(double jd, double J, MLBZ &ob)
{ // 命理八字计算。jd为格林尼治UT(J2000起算),J为本地经度,返回在物件ob中
    int i, v;
    long double c;
    long double jd2 = jd + dt_T(jd);                         // 力学时
    long double w = XL::S_aLon(jd2 / 36525.0, -1);           // 此刻太阳视黄经
    int k = int2((w / pi2 * 360 + 45 + 15 * 360) / 30); // 1984年立春起算的节气数(不含中气)
    jd += pty_zty2(jd2 / 36525) + J / M_PI / 2;         // 本地真太阳时(使用低精度算法计算时差)
    ob.bz_zty = JD::timeStr(jd);

    jd += 13.0 / 24;                             // 转为前一日23点起算(原jd为本日中午12点起算)
    int D = floor(jd), SC = int2((jd - D) * 12); // 日数与时辰

    v = int2(k / 12.0 + 6000000);
    ob.bz_jn = std::string(Gan[v % 10]) + Zhi[v % 12];
    v = k + 2 + 60000000;
    ob.bz_jy = std::string(Gan[v % 10]) + Zhi[v % 12];
    v = D - 6 + 9000000;
    ob.bz_jr = std::string(Gan[v % 10]) + Zhi[v % 12];
    v = (D - 1) * 12 + 90000000 + SC;
    ob.bz_js = std::string(Gan[v % 10]) + Zhi[v % 12];

    v -= SC, ob.bz_JS = ""; // 全天纪时表
    for (i = 0; i < 13; i++)
    {                                                                       // 一天中包含有13个纪时
        std::string c = std::string(Gan[(v + i) % 10]) + Zhi[(v + i) % 12]; // 各时辰的八字
        if (SC == i)
        {
            ob.bz_js = c, c = "\033[31m" + c + "\033[0m"; // 红色显示这时辰
        }
        ob.bz_JS += (i ? " " : "") + c;
    }
}

void ML_calc(Time &dat)
{
    MLBZ ob = {};
    GeoPostion &gep = GeoPostion::getInstance();
    JINGWEI jw = gep.getCityGeoPos();
    double jd = JD::toJD({dat.Y, dat.M, dat.D, dat.h, dat.m, dat.s});
    mingLiBaZi(jd + (-8.0) / 24 - J2000, jw.J / radd, ob); // 八字计算

    std::cout << "\033[31;1m[日标]：\033[0m"
              << "公历 " << dat.Y << "-" << dat.M << "-" << dat.D << " 儒略日数 " << int2(jd + 0.5) << " 距2000年首" << int2(jd + 0.5 - J2000) << "日\n"
              << "\033[31;1m[八字]：\033[0m" << ob.bz_jn << "年 " << ob.bz_jy << "月 " << ob.bz_jr << "日 " << ob.bz_js << "时 真太阳 \033[31m" << ob.bz_zty << "\033[0m"
              << "\n\033[1;32m[纪时]：\033[0m" << ob.bz_JS << "\n"
              << "\033[1;32m[时标]：\033[0;1;35m"
              << "23　 01　 03　 05　 07　 09　 11　 13　 15　 17　 19　 21　 23\033[0m" << std::endl;
}

bool almostEqual(long double a, long double b)
{
    return std::abs(a - b) < 1e-12;
}

void jqCalc(int year, int month)
{
    Time t0 = {year, month, 1, 12, 0, 0.1};
    auto Bd0 = int2(JD::toJD(t0)) - J2000; // 公历月首,中午
    t0.M++;
    if (t0.M > 12)
    {
        t0.Y++;
        t0.M = 1;
    }
    auto Bdn = int2(JD::toJD(t0)) - J2000 - Bd0; // 本月天数(公历)
    int d0 = Bd0;
    SSQ &ssq = SSQ::getInstance();
    std::vector<long double> vecZQ = ssq.getZhongQi();
    std::string Ljq;
    bool hasJie = false;
    for (size_t i = 0, j = 0; i < Bdn; i++)
    {
        d0 = Bd0 + i;

        // 如果d0已在计算农历范围内则不再计算
        if (vecZQ.empty() || d0 < vecZQ[0] || d0 >= vecZQ[24])
        {
            ssq.calcY(d0);
            vecZQ = ssq.getZhongQi();
        }

        auto qk = int2((d0 - vecZQ[0] - 7) / 15.2184);
        if (qk < 23 && d0 >= vecZQ[qk + 1])
        {
            qk++; // 节气的取值范围是0-23
        }

        //奇数对应是节，偶数对应就是气
        if (d0 == vecZQ[qk] && (qk & 1) == 1)
        {
            Ljq = Jqmc[qk];
            hasJie = true;
            break;
        }
    }

    if (hasJie)
    {
        auto jqd = qi_accurate2(d0, false, 120) + J2000;
        auto jqt = JD::JD2DD(jqd);
        std::cout << year << "年" << month << "月节令是:" << Ljq << ",时间:" << JD::timeStr(jqt) << std::endl;
    }
}