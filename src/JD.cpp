#include "JD.h"
#include "const.h"
#include <cstring>
#include <cstdio>
#include <cmath>
#include <ctime>
#include "eph.h"

// https://github.com/guolisen/YiEngine/blob/2ce67dc91fd5fea8e394a5af60dc1e56c5044452/src/DateTime/JulianDay.cpp
/*
《儒略历》是现今国际通用的公历的前身。西方国家16世纪大多采用它。公元前46年，罗马统帅儒略·恺撒在埃及亚历山大的希腊数学家兼天文学家索西琴尼的帮助下制订的，
 并在公元前46年1月1日起执行实行，取代旧罗马历法的一种历法。所以人们就把这一历法称为《儒略历》。
《儒略历》以回归年为基本单位，是一部纯粹的阳历。它将全年分设为12个月，单数月是大月，长31日，双月是小月，长为30日，只有2月平年是29日，闰年30日。每年设365日，
 每四年一闰，闰年366日，每年平均长度是365.25日。《儒略历》编制好后，儒略·恺撒的继承人奥古斯都又从2月减去一日加上8月（8月的拉丁名即他的名字奥古斯都》，
 又把9月、11月改为小月，10月、12月改为大月。
《儒略历》比回归年365.2422日长0.0078日，400年要多出3.12日。从公元325年定春分为3月21日提早到了3月11日。1500年后由于误差较大，
 被罗马教皇格里高利十三世于1582年进行改善与修订，变为格里历（Gregorian calendar），即沿用至今的世界通用的公历。
格里高利历 1582年10月4日之前，应用的是儒略历。
罗马教皇格里高利十三世设立了改革历法的专门委员会，比较了各种方案后， 决定采用意大利医生利里奥的方案，在400年中去掉儒略历多出的三个闰年。
1582年3月1日，格里高利颁发了改历命令，内容是：
一、1582年10月4日后的一天是10月15日，而不是10月5日，但星期序号仍然连续计算，10月4日是星期四，第二天10月15日是星期五。这样，就把从公元325年以来积累的老账一笔勾销了。
二、为避免以后再发生春分飘离的现象，改闰年方法为： 凡公元年数能被4整除的是闰年，但当公元年数后边是带两个“0”的“世纪年”时，必须能被400整除的年才是闰年。
格里高利历的历年平均长度为365日5时49分12秒，比回归年长26秒
 */
// 公历转儒略日
long double JD::DD2JD(int y, uint8_t m, long double d)
{
    int n = 0, G = 0;
    // 判断是否为格里高利历日1582*372+10*31+15
    if (y * 372 + m * 31 + static_cast<int64_t>(d) >= 588829)
    {
        G = 1;
    }
    if (m <= 2)
    {
        m += 12, y--;
    }
    // 加百年闰
    if (G)
    {
        n = int2(y / 100), n = 2 - n + int(n / 4);
    }

    return int2(365.25 * (y + 4716)) + int2(30.6001 * (m + 1)) + d + n - 1524.5;
}

// 儒略日数转公历
Time JD::JD2DD(long double jd)
{
    Time r;
    int D = int2(jd + 0.5);
    long double F = jd + 0.5 - D, c; // 取得日数的整数部份A及小数部分F
    // 1582-10-15,[明]神宗 朱翊钧 万历10年 农历壬午年【马年】九月十九，甲戌日, 格里高利日
    if (D >= 2299161)
    {
        c = int((D - 1867216.25) / 36524.25), D += 1 + c - int2(c / 4);
    }

    // 年数
    D += 1524;
    r.Y = int2((D - 122.1) / 365.25);

    // 月数
    D -= int2(365.25 * r.Y);
    r.M = int2(D / 30.601);

    // 日数
    D -= int2(30.601 * r.M);
    r.D = D;

    if (r.M > 13)
    {
        r.M -= 13, r.Y -= 4715;
    }
    else
    {
        r.M -= 1, r.Y -= 4716;
    }
    // 日的小数转为时分秒
    F *= 24;
    r.h = int2(F);
    F -= r.h;
    F *= 60;
    r.m = int2(F);
    F -= r.m;
    F *= 60;
    r.s = F;
    return r;
}

long double JD::toJD(const Time &time)
{
    return JD::DD2JD(time.Y, time.M, time.D + ((time.s / 60 + time.m) / 60 + time.h) / 24);
}

// 提取jd中的时间(去除日期);
std::string JD::timeStr(long double jd)
{
    int h, m, s;
    jd += 0.5;
    jd = (jd - int2(jd));
    s = int2(jd * 86400 + 0.5);
    h = int2(s / 3600);
    s -= h * 3600;
    m = int2(s / 60);
    s -= m * 60;
    std::string ret;
    char buff[11];
    memset(buff, 0, 11);
    snprintf(buff, 11, "0%d", h);
    ret.append(buff + strlen(buff) - 2);
    ret += ":";

    memset(buff, 0, 11);
    snprintf(buff, 11, "0%d", m);
    ret.append(buff + strlen(buff) - 2);
    ret += ":";

    memset(buff, 0, 11);
    snprintf(buff, 11, "0%d", s);
    ret.append(buff + strlen(buff) - 2);

    return ret;
}

std::string JD::timeStr(const Time &t)
{
    std::string
        Y = "     " + std::to_string(t.Y),
        M = "0" + std::to_string(t.M),
        D = "0" + std::to_string(t.D);

    int h = int(t.h), m = int(t.m), s = int2(t.s + .5);
    if (s >= 60)
    {
        s -= 60, m++;
    }
    if (m >= 60)
    {
        m -= 60, h++;
    }

    std::string _h, _m, _s;
    _h = "0" + std::to_string(h);
    _m = "0" + std::to_string(m);
    _s = "0" + std::to_string(s);
    Y = Y.substr(Y.length() - 5, 5);
    M = M.substr(M.length() - 2, 2);
    D = D.substr(D.length() - 2, 2);
    _h = _h.substr(_h.length() - 2, 2);
    _m = _m.substr(_m.length() - 2, 2);
    _s = _s.substr(_s.length() - 2, 2);

    return Y + "-" + M + "-" + D + " " + _h + ":" + _m + ":" + _s;
}

std::string JD::formatStr(long double jd)
{
    Time t0 = JD2DD(jd);
    auto Y = t0.getYear();
    auto M = t0.getMonth();
    auto D = t0.getDay();
    auto h = t0.getHour();
    auto m = t0.getMin();
    auto s = t0.getSec();
    if (t0.getHour() < 13)
    {
        h = t0.getHour();
    }
    else
    {
        h = t0.getHour() - 12;
    }

    auto hm = t0.getHour() * 100 + m;
    std::string sMeridiem;
    if (hm < 600)
    {
        sMeridiem = "凌晨";
    }
    else if (hm < 900)
    {
        sMeridiem = "早上";
    }
    else if (hm < 1130)
    {
        sMeridiem = "上午";
    }
    else if (hm < 1230)
    {
        sMeridiem = "中午";
    }
    else if (hm < 1800)
    {
        sMeridiem = "下午";
    }
    else
    {
        sMeridiem = "晚上";
    }

    std::string res = std::to_string(Y) + "年" + std::to_string(M) + "月" + std::to_string(D) + "日" +
            sMeridiem + std::to_string(int(h)) + "点" +
            std::to_string(int(m)) + "分" +
            std::to_string(int(s)) + "秒";
    return res;
}

std::string JD::JD2str(long double jd)
{
    Time r = JD2DD(jd);
    return timeStr(r);
}

Time JD::getNowTime()
{
    struct tm *bjs;
    time_t time0;
    time0 = time(nullptr);
    bjs = localtime(&time0);
    Time t =
        {
            bjs->tm_year + 1900,
            bjs->tm_mon + 1,
            bjs->tm_mday,
            (double)bjs->tm_hour,
            (double)bjs->tm_min,
            (double)bjs->tm_sec};
    return t;
}

// 计算真太阳时
std::tuple<long double, Time, std::string> JD::calcAST(Time &dt, long double lon)
{
    auto utcJd = toJD(dt) - 8.0 / 24.0;              // utc时间
    auto tdJd = dt_T(utcJd - J2000) + utcJd - J2000; // 力学时
    utcJd += pty_zty2(tdJd / 36525) + lon / (360 - 0.0);
    auto dtAst = JD::JD2DD(utcJd);
    auto zty = JD::timeStr(utcJd);
    return std::make_tuple(utcJd, dtAst, zty);
}

long double JD::calcAST(long double jd, long double lon)
{
    auto utcJd = jd - 8.0 / 24.0;                    // utc时间
    auto tdJd = dt_T(utcJd - J2000) + utcJd - J2000; // 力学时
    utcJd += pty_zty2(tdJd / 36525) + lon / (360 - 0.0);
    return utcJd;
}

double mod(double x, double y)
{
    return x - y * std::floor(x / y);
}

std::tuple<double, double> JD::gcal2jd(int year, int month, int day)
{

    int64_t a = static_cast<int64_t>((month - 14) / 12.0);
    int64_t jd0 = static_cast<int64_t>((1461 * (year + 4800 + a)) / 4.0);
    jd0 += static_cast<int64_t>((367 * (month - 2 - 12 * a)) / 12.0);
    int64_t x = static_cast<int64_t>((year + 4900 + a) / 100.0);
    jd0 -= static_cast<int>((3 * x) / 4.0);
    double jd = jd0 + day - 2432075.5;

    jd -= 0.5;

    return std::make_tuple(MJD_0, jd);
}

std::tuple<int, int, int, double> JD::jd2gcal(double jd1, double jd2)
{
    double jd1_f, jd2_f;
    double jd1_i, jd2_i;

    jd1_f = modf(jd1, &jd1_i);
    jd2_f = modf(jd2, &jd2_i);

    double jd_i = jd1_i + jd2_i;

    double f = jd1_f + jd2_f;

    if (-0.5 < f && f < 0.5)
    {
        f += 0.5;
    }
    else if (f >= 0.5)
    {
        jd_i += 1;
        f -= 0.5;
    }
    else if (f <= -0.5)
    {
        jd_i -= 1;
        f += 1.5;
    }

    double ell = jd_i + 68569;
    int64_t n = static_cast<int64_t>((4 * ell) / 146097.0);
    ell -= static_cast<int64_t>((146097 * n + 3) / 4.0);
    int64_t i = static_cast<int64_t>((4000 * (ell + 1)) / 1461001);
    ell -= static_cast<int64_t>((1461 * i) / 4.0) - 31;
    int64_t j = static_cast<int64_t>((80 * ell) / 2447.0);
    int day = static_cast<int>(ell - static_cast<int64_t>((2447 * j) / 80.0));
    ell = static_cast<int64_t>(j / 11.0);
    int month = static_cast<int>(j + 2 - (12 * ell));
    int year = static_cast<int>(100 * (n - 49) + i + ell);

    return std::make_tuple(year, month, day, f);
}

std::tuple<int, int, int, double> JD::jd2jcal(double jd1, double jd2)
{
    double jd1_f, jd2_f;
    double jd1_i, jd2_i;

    jd1_f = modf(jd1, &jd1_i);
    jd2_f = modf(jd2, &jd2_i);

    double jd_i = jd1_i + jd2_i;

    double f = jd1_f + jd2_f;

    if (-0.5 < f && f < 0.5)
    {
        f += 0.5;
    }
    else if (f >= 0.5)
    {
        jd_i += 1;
        f -= 0.5;
    }
    else if (f <= -0.5)
    {
        jd_i -= 1;
        f += 1.5;
    }

    double j = jd_i + 1402.0;
    int64_t k = static_cast<int64_t>((j - 1) / 1461.0);
    double ell = j - (1461.0 * k);
    int64_t n = static_cast<int64_t>((ell - 1) / 365.0) - static_cast<int64_t>(ell / 1461.0);
    double i = ell - (365.0 * n) + 30.0;
    double j2 = static_cast<int64_t>((80.0 * i) / 2447.0);
    double day = i - static_cast<int64_t>((2447.0 * j2) / 80.0);
    double i2 = static_cast<int64_t>(j2 / 11.0);
    double month = j2 + 2 - (12.0 * i2);
    double year = (4 * k) + n + i2 - 4716.0;

    return std::make_tuple(static_cast<int>(year), static_cast<int>(month), static_cast<int>(day), f);
}

std::tuple<double, double> JD::jcal2jd(int year, int month, int day)
{
    double jd = 367 * year;
    int x = static_cast<int>((month - 9) / 7.0);
    jd -= static_cast<int>((7 * (year + 5001 + x)) / 4.0);
    jd += static_cast<int>((275 * month) / 9.0);
    jd += day;
    jd += 1729777 - 2400000.5;

    jd -= 0.5;

    return std::make_tuple(MJD_0, jd);
}
