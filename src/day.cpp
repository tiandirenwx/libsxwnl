#include "day.h"
#include "eph.h"
#include <memory>

namespace sxtwl
{
    GZ getShiGz(uint8_t dayTg, uint8_t hour, bool isZaoWanZiShi = true);
};


Day::Day(int d,double jf)
{
    this->d0 = d;
    this->jdF = jf;
    this->bd0 = -2;
    this->bdn = 0;

    this->lunar_total_days_ = 0;
    this->solar_month_ = 0;
    this->jieling_ = -2;
    this->lunar_lichun_year_ = 0;
    this->lunar_jun_year_ = 0;
    this->huangdi_year_ = 0;

    this->moslem_month_ = 0;

    this->gz_lichun_year_ = nullptr;
    this->gz_jan_year_ = nullptr;
    this->gz_month_ = nullptr;
    this->gz_day_ = nullptr;

    this->week_ = 0xFF;
    this->week0_ = 0xFF;
    this->weeki_ = 0xFF;
    this->weekn_ = 0xFF;

    this->xingzuo_ = 0xFF;
    this->jieqi_jd_ = 0;

    this->yx_idx_ = -2;
    this->yxjd_ = 0;
}

Day::~Day()
{

    if (this->gz_lichun_year_)
    {
        delete this->gz_lichun_year_;
        this->gz_lichun_year_ = nullptr;
    }

    if (this->gz_jan_year_)
    {
        delete this->gz_jan_year_;
        this->gz_jan_year_ = nullptr;
    }

    if (this->gz_month_)
    {
        delete this->gz_month_;
        this->gz_month_ = nullptr;
    }

    if (this->gz_day_)
    {
        delete this->gz_day_;
        this->gz_day_ = nullptr;
    }
}

void Day::checkSSQ()
{
    SSQ &SSQPtr = SSQ::getInstance();
    std::vector<long double> vecZQ = SSQPtr.getZhongQi();
    if (vecZQ.empty() || this->d0 < vecZQ[0] || this->d0 >= vecZQ[24])
    {
        SSQPtr.calcY(this->d0);
    }
}

/**
 * 确定已经计算过阴历信息
 */
void Day::checkLunarData()
{
    // 已经计算过了
    if (this->lunar_total_days_ != 0)
    {
        return;
    }
    this->checkSSQ();
    SSQ &SSQPtr = SSQ::getInstance();
    std::vector<long double> vecZQ = SSQPtr.getZhongQi();
    std::vector<int> vecHS = SSQPtr.getHS();
    std::vector<int> vecYm = SSQPtr.getYm();
    std::vector<int> vecDx = SSQPtr.getDx();
    std::vector<bool> vecSpecMonth = SSQPtr.getSpecificLunarMonth();
    int mk = int2((this->d0 - vecHS[0]) / 30);
    if (mk < 13 && vecHS[mk + 1] <= this->d0)
    {
        mk++; // 农历所在月的序数
    }

    // 阴历所处的日
    this->lunar_day_ = this->d0 - vecHS[mk];

    this->cur_dz_ = int(this->d0 - vecZQ[0]);  // 距冬至的天数
    this->cur_xz_ = int(this->d0 - vecZQ[12]); // 距夏至的天数
    this->cur_lq_ = int(this->d0 - vecZQ[15]); // 距立秋的天数
    this->cur_mz_ = int(this->d0 - vecZQ[11]); // 距芒种的天数
    this->cur_xs_ = int(this->d0 - vecZQ[13]); // 距小暑的天数

    // 月的信息
    this->checkExtSolarData();
    if (this->d0 == vecHS[mk] || this->d0 == this->bd0)
    {
        int leap = SSQPtr.getLeap();
        this->lunar_month_ = vecYm[mk];                         // 月名称
        this->lunar_total_days_ = vecDx[mk];                    // 月大小
        this->is_lunar_leap_month_ = (leap != 0 && leap == mk); // 闰状况
        this->is_lunar_spec_month_ = vecSpecMonth[mk];          // 是否是特殊的重月
    }
    else
    {
        auto d = this->d0 - 1;
        auto yesterday = new Day(d);
        yesterday->checkSSQ();
        SSQ &ssq = SSQ::getInstance();
        std::vector<int> vecYesHS = ssq.getHS();
        std::vector<int> vecYesYm = ssq.getYm();
        std::vector<int> vecYesDx = ssq.getDx();
        std::vector<bool> vecYesSpecMonth = ssq.getSpecificLunarMonth();
        int ymk = int2((d - vecYesHS[0]) / 30);
        if (ymk < 13 && vecYesHS[ymk + 1] <= d)
        {
            ymk++; // 农历所在月的序数
        }
        int yesleap = ssq.getLeap();
        this->lunar_month_ = vecYesYm[ymk];                            // 月名称
        this->lunar_total_days_ = vecYesDx[ymk];                       // 月大小
        this->is_lunar_leap_month_ = (yesleap != 0 && yesleap == ymk); // 闰状况
        this->is_lunar_spec_month_ = vecYesSpecMonth[ymk];             // 是否是特殊的重月
        delete yesterday;
    }
}

void Day::checkSolarData()
{
    if (this->solar_month_ != 0)
    {
        return;
    }

    Time t = JD::JD2DD(this->d0 + J2000);
    this->solar_year_ = t.Y;
    this->solar_day_ = t.D;
    this->solar_month_ = t.M;
}

void Day::checkExtSolarData()
{
    if (this->bd0 != -2)
    {
        return;
    }

    Time t{};
    t.h = 12, t.m = 0, t.s = 0.1;
    t.Y = this->getSolarYear();
    t.M = this->getSolarMonth();
    t.D = 1;
    this->bd0 = int2(JD::toJD(t)) - J2000;
    t.M++;
    if (t.M > 12)
    {
        t.Y++;
        t.M = 1;
    }

    this->bdn = int2(JD::toJD(t)) - J2000 - this->bd0;
}

void Day::checkMoslemData()
{
    if (this->moslem_month_ != 0)
    {
        return;
    }

    int d = d0 + 503105;
    int z = int2(d / 10631); // 10631为一周期(30年)
    d -= z * 10631;
    int y = int2((d + 0.5) / 354.366); // 加0.5的作用是保证闰年正确(一周中的闰年是第2,5,7,10,13,16,18,21,24,26,29年)
    d -= int2(y * 354.366 + 0.5);
    int m = int2((d + 0.11) / 29.51); // 分子加0.11,分母加0.01的作用是第354或355天的的月分保持为12月(m=11)
    d -= int2(m * 29.5 + 0.5);

    this->moslem_year_ = z * 30 + y + 1;
    this->moslem_month_ = m + 1;
    this->moslem_day_ = d + 1;
}

/**
 * 计算节气数据
 */
void Day::checkJQData()
{
    if (this->jieling_ != -2)
    {
        return;
    }

    this->jieling_ = -1;
    this->getJieQiJD();
}
void Day::checkYXData()
{
    if (this->yx_idx_ != -2)
    {
        return;
    }

    this->yx_idx_ = -1;
    this->getYueXiangJD();
}

Day *Day::after(int day)
{
    return new Day(this->d0 + day);
}

Day *Day::before(int day)
{
    return new Day(this->d0 - day);
}

double Day::getJulianDate()
{
    return this->d0 + this->jdF + J2000;
}

/**
 * 获取阴历日期
 */
int Day::getLunarDay()
{
    this->checkLunarData();
    return this->lunar_day_ + 1;
}

/**
 * 获取阴历月
 */
uint8_t Day::getLunarMonth()
{
    this->checkLunarData();
    static const int yueIndex[13] = {11, 12, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1};
    return yueIndex[this->lunar_month_];
}

int Day::getHuangdiYear()
{
    return this->huangdi_year_;
}

int Day::getLunarYear(bool chineseNewYearBoundary)
{
    SSQ &SSQPtr = SSQ::getInstance();
    if (this->lunar_lichun_year_ == 0 || this->lunar_jun_year_ == 0)
    {
        this->checkSSQ();
    }
    std::vector<long double> vecZQ = SSQPtr.getZhongQi();
    std::vector<int> vecHS = SSQPtr.getHS();
    std::vector<int> vecYm = SSQPtr.getYm();
    int iLeap = SSQPtr.getLeap();
    long double D = vecZQ[3] + (this->d0 < vecZQ[3] ? -365 : 0) + 365.25 * 16 - 35; // 以立春为界定纪年
    this->lunar_lichun_year_ = int2(D / 365.2422 + 0.5);

    D = vecHS[2]; // 一般第3个月为春节
    for (int j = 0; j < 14; j++)
    { // 找春节
        if (vecYm[j] != 2 || iLeap == j && j)
        {
            continue;
        }
        D = vecHS[j];
        if (this->d0 < D)
        {
            D -= 365;
            break;
        } // 无需再找下一个正月
    }
    D = D + 5810; // 计算该年春节与1984年平均春节(立春附近)相差天数估计
    this->lunar_jun_year_ = int2(D / 365.2422 + 0.5);
    this->huangdi_year_ = this->lunar_jun_year_ + 1984 + 2698; // 黄帝纪年

    // 以春节为界
    if (chineseNewYearBoundary)
    {
        return this->lunar_jun_year_ + 1984;
    }

    //立春为界
    return this->lunar_lichun_year_ + 1984;
}

GZ Day::getYearGZ(bool chineseNewYearBoundary)
{
    // 以春节为界
    if (chineseNewYearBoundary)
    {
        if (this->gz_jan_year_ == nullptr)
        {
            int year = this->getLunarYear(chineseNewYearBoundary) - 1984;
            int D = year + 12000;
            this->gz_jan_year_ = new GZ(D % 10, D % 12);
        }
        return *(this->gz_jan_year_);
    }

    // 以立春为界
    if (this->gz_lichun_year_ == nullptr)
    {
        int year = this->getLunarYear(false) - 1984;
        int D = year + 12000;
        this->gz_lichun_year_ = new GZ(D % 10, D % 12);
    }
    return *(this->gz_lichun_year_);
}

GZ Day::getMonthGZ()
{
    SSQ &SSQPtr = SSQ::getInstance();
    if (nullptr == this->gz_month_)
    {
        this->checkSSQ();
        std::vector<long double> vecZQ = SSQPtr.getZhongQi();
        int mk = int2((this->d0 - vecZQ[0]) / 30.43685);
        // 相对大雪的月数计算,mk的取值范围0-12
        if (mk < 12 && this->d0 >= vecZQ[2 * mk + 1])
            mk++;
        // 相对于1998年12月7(大雪)的月数,900000为正数基数
        int D = mk + int2((vecZQ[12] + 390) / 365.2422) * 12 + 900000;
        this->gz_month_ = new GZ(D % 10, D % 12);
    }
    return *(this->gz_month_);
}

GZ Day::getDayGZ()
{
    if (nullptr == this->gz_day_)
    {
        int D = this->d0 - 6 + 9000000;
        this->gz_day_ = new GZ(D % 10, D % 12);
    }
    return *(this->gz_day_);
}

GZ Day::getHourGZ(uint8_t hour, bool isZaoWanZiShi)
{
    GZ dayGZ = this->getDayGZ();
    return sxtwl::getShiGz(dayGZ.tg, hour, isZaoWanZiShi);
}

bool Day::isLunarLeap()
{
    this->checkLunarData();
    return this->is_lunar_leap_month_;
}

bool Day::isSpecNextMonth()
{
    this->checkLunarData();
    return is_lunar_spec_month_;
}

int Day::getSolarYear()
{
    this->checkSolarData();
    return this->solar_year_;
}

uint8_t Day::getSolarMonth()
{
    this->checkSolarData();
    return this->solar_month_;
}

int Day::getSolarDay()
{
    this->checkSolarData();
    return this->solar_day_;
}

int Day::getMoslemYear()
{
    this->checkMoslemData();
    return this->moslem_year_;
}

uint8_t Day::getMoslemMonth()
{
    this->checkMoslemData();
    return this->moslem_month_;
}

int Day::getMoslemDay()
{
    this->checkMoslemData();
    return this->moslem_day_;
}

uint8_t Day::getWeek()
{
    if (this->week_ == 0xFF)
    {
        this->week_ = (this->d0 + J2000 + 1 + 7000000) % 7; // 当前日的星期
    }
    return this->week_;
}

uint8_t Day::getFirstWeekDayOfMonth()
{
    if (this->week0_ == 0xFF)
    {
        checkExtSolarData();
        this->week0_ = (this->bd0 + J2000 + 1 + 7000000) % 7; // 本月第一天星期
    }

    return this->week0_;
}

uint8_t Day::getTotalWeekNumsOfMonth()
{
    if (this->weekn_ == 0xFF)
    {
        checkExtSolarData();
        this->weekn_ = int2((this->week0_ + this->bdn - 1) / 7) + 1; // 本月的总周数
    }
    return this->weekn_;
}

// 处于该月的第几周
uint8_t Day::getWeekIndex()
{
    if (this->weeki_ == 0xFF)
    {
        int i = (this->getSolarDay() - 1) % 7;

        int w0 = 0;
        if (this->getWeek() >= i)
        {
            w0 = this->getWeek() - i;
        }
        else
        {
            w0 = this->getWeek() + 7 - i;
        }

        this->weeki_ = int2((w0 + this->getSolarDay() - 1) / 7) + 1;
    }
    return this->weeki_;
}

bool Day::hasYueXiang()
{
    this->checkYXData();
    return this->yx_idx_ != -1;
}

uint8_t Day::getYueXiang()
{
    this->checkYXData();
    return this->yx_idx_;
}

double Day::getYueXiangJD()
{
    if (this->yxjd_ != 0)
    {
        return this->yxjd_ + J2000;
    }

    checkExtSolarData();

    // 以下是月相与节气的处理
    long double d, jd2 = this->bd0 + dt_T(this->bd0) - 8.0 / 24.0;
    // 月相查找
    long double w = XL::MS_aLon(jd2 / 36525, 10, 3);
    w = int2((w - 0.78) / M_PI * 2) * M_PI / 2;
    int D = 0;
    do
    {
        d = so_accurate(w);
        D = int2(d + 0.5);
        auto xn = int2(w / pi2 * 4 + 4000000.01) % 4;
        w += pi2 / 4;
        if (D >= this->bd0 + this->bdn)
        {
            break;
        }
        if (D < this->bd0)
        {
            continue;
        }
        if (D == this->d0)
        {
            this->yxjd_ = d;
            this->yx_idx_ = xn;
            break;
        }
    } while (D + 5 < this->bd0 + this->bdn);
    return this->yxjd_ + J2000;
}
// 是否有节气
bool Day::hasJieQi()
{
    this->checkJQData();
    return this->jieling_ != -1;
}
// 获取节气
uint8_t Day::getJieQi()
{
    this->checkJQData();
    return this->jieling_;
}

double Day::getJieQiJD()
{
    if (this->jieqi_jd_ != 0)
    {
        return this->jieqi_jd_ + J2000;
    }

    checkExtSolarData();

    long double d, xn, jd2 = this->bd0 + dt_T(this->bd0) - 8.0 / 24.0;
    long double w = XL::S_aLon(jd2 / 36525, 3);
    w = int2((w - 0.13) / pi2 * 24) * pi2 / 24;
    int D = 0;

    do
    {
        d = qi_accurate(w);
        D = int2(d + 0.5);
        // 计算出的节令值
        xn = int2(w / pi2 * 24 + 24000006.01) % 24;
        w += pi2 / 24;
        if (D > this->bd0 + this->bdn)
        {
            break;
        }
        if (D < this->bd0)
        {
            continue;
        }

        if (D == this->d0)
        {
            this->jieqi_jd_ = d;
            this->jieling_ = xn;
            break;
        }
    } while (D + 12 < this->bd0 + this->bdn);

    return this->jieqi_jd_ + J2000;
}

// 获取星座
uint8_t Day::getConstellation()
{
    SSQ &SSQPtr = SSQ::getInstance();
    if (this->xingzuo_ == 0xFF)
    {
        this->checkSSQ();
        std::vector<long double> vecZQ = SSQPtr.getZhongQi();
        int mk = int2((this->d0 - vecZQ[0] - 15) / 30.43685);
        // 星座所在月的序数,(如果j=13,ob.d0不会超过第14号中气)
        if (mk < 11 && this->d0 >= vecZQ[2 * mk + 2])
        {
            mk++;
        }
        this->xingzuo_ = (mk + 12) % 12;
    }
    return this->xingzuo_;
}

uint8_t Day::getLunarDay12Jian() {
    GZ gzMonth = getMonthGZ();
    GZ gzDay = getDayGZ();
    uint8_t posYueJian = gzMonth.dz;
    uint8_t posRiZhi = gzDay.dz;

    this->lunar_day_12jian_ = (12 + posRiZhi - posYueJian) % 12;

    return this->lunar_day_12jian_;
}

Day *Day::fromDate(const Time &t)
{
    Time t0{};
    t0.h = 12, t0.m = 0, t0.s = 0.1;
    t0.Y = t.getYear();
    t0.M = t.getMonth();
    t0.D = t.getDay();
    auto d = int2(JD::toJD(t0)) - J2000;
    auto jd = JD::toJD(t) - J2000;
    return new Day(d,jd - d);
}

Day *Day::fromJulianDay(double jd)
{
    auto t = JD::JD2DD(jd);
    t.h = 12,t.m = 0, t.s = 0.1;
    auto d = int2(JD::toJD(t)) - J2000;
    auto jf = jd - d - J2000;
    return new Day(d,jf);
}

Day *Day::fromSolar(int _year, uint8_t _month, int _day)
{
    Time t{};
    t.h = 12, t.m = 0, t.s = 0.1;
    t.Y = _year;
    t.M = _month;
    t.D = _day;
    auto d0 = int2(JD::toJD(t)) - J2000;
    return new Day(d0);
}

// 19年7闰法计算逻辑
Day *Day::getLunarDateSpecial(int year, uint8_t month, int day, bool isRun, bool isSpec)
{
    //-721年12月17日至-479年11月春秋战国历采用19年7闰,闰年的末月置闰并取名闰十三。年首为正月。
    //-221年10月31日至-104年11月秦汉历同样采用19年7闰，闰年的末月置闰并取名后九月。年首为十月。
    //本函数只处理这段时间的农历转阳历
    Time t{};
    t.h = 12, t.m = 0, t.s = 0.1;
    t.Y = year;
    t.M = 2;
    t.D = 4;

    int yearHeadMonth = 1;

    int Bd0 = int2(JD::toJD(t)) - J2000;
    if (Bd0 + J2000 >= 1640641) //-221年10月31日
    {
        yearHeadMonth = 10;
    }

    //月首为十月，当是12月时处于交界情况，往下一年计算月序，但查找时不要往后查
    bool findNextFlag = true;
    if (10 == yearHeadMonth && 12 == month)
    {
        t.Y++;
        Bd0 = int2(JD::toJD(t)) - J2000;
        findNextFlag = false;
    }

    SSQ& SSQPtr = SSQ::getInstance();
    std::vector<long double> vecZQ = SSQPtr.getZhongQi();
    if (vecZQ.empty() || Bd0 < vecZQ[0] || Bd0 >= vecZQ[24])
    {
        SSQPtr.calcY(Bd0);
    }

    int leap = SSQPtr.getLeap();
    std::vector<int> vecHS = SSQPtr.getHS();
    std::vector<int> vecYm = SSQPtr.getYm();
    std::vector<int> vecDx = SSQPtr.getDx();
    std::vector<int> vecYueName = SSQPtr.getYueName();
    std::vector<int> vecNs = SSQPtr.getBDNS();

    int jd = 0;
    if(isRun)
    {
        if (leap > 0 && month == vecYueName[leap])
        {
            if (0 < day && day <= vecDx[leap])
            {
                jd = vecHS[leap] + day - 1;
            } else {
                return nullptr;
            }
        }
    }
    else
    {
        long startIdx = 0;
        bool bWithYearHeadFlag = false;
        auto it = std::find(vecYueName.begin(), vecYueName.end(), yearHeadMonth);
        if (it != vecYueName.end())
        {
            startIdx = std::distance(vecYueName.begin(), it);
            if (startIdx <= 2)
            {
                bWithYearHeadFlag = true;
            }
        }

        //春秋战国历年首临界情况特殊处理
        if ((month == 1 && year < -220) && !bWithYearHeadFlag )
        {
            vecHS.insert(vecHS.begin(),vecNs[0]);
            vecDx.insert(vecDx.begin(),vecHS[1] - vecHS[0]);
            vecYueName.insert(vecYueName.begin(),yearHeadMonth);
        }

        //秦汉历年首十月特殊处理: 往下一年推算月序时，十冬腊这些月先查前面的，否则往后查
        auto curStartIdx = vecYueName.begin();
        if (vecYueName[0] >= 10 && month >= 10 && findNextFlag)
        {
            //查据正月的距离
            it = std::find(curStartIdx , vecYueName.end(), 1);
            if (it != vecYueName.end())
            {
                startIdx = std::distance(vecYueName.begin(), it);
            } else
            {
                return nullptr;
            }
            curStartIdx += startIdx;
        }

        //从start开始找月序中有没有当前输入的月
        it = std::find(curStartIdx , vecYueName.end(), month);
        if (it != vecYueName.end())
        {
            startIdx = std::distance(vecYueName.begin(), it);
        }
        else
        {
            //表明输入的月还是有问题
            return nullptr;
        }

        if(isSpec)
        {
            it = std::find(vecYueName.begin() + startIdx + 1, vecYueName.end(), month);
            if (it != vecYueName.end())
            {
                startIdx = std::distance(vecYueName.begin(), it);
            }
            else
            {
                //表明输入的月还是有问题
                return nullptr;
            }
        }
        if (0 < day && day <= vecDx[startIdx] ){
            jd = vecHS[startIdx] + day - 1;
        } else {
            return nullptr;
        }
    }

    return new Day(jd);
}

// 农历转公历核心处理
Day *Day::getLunarDate(int year, uint8_t month, int day, bool isRun, bool isSpec)
{
    if (year < -722 )
    {
        return nullptr;
    }

    if (year < -104 && year > -722)
    {
        return getLunarDateSpecial(year,month,day,isRun,isSpec);
    }

    Time t{};
    t.h = 12, t.m = 0, t.s = 0.1;
    t.Y = year;
    t.M = 2;
    t.D = 4;

    // 正常来说一个公历年只有10个农历月，但有些特殊的年份会含有多个农历月
    uint8_t stdMonth = 0;
    if ((9 <= year && year <= 23) || (237 <= year && year <= 239))
    {
        stdMonth = 11;
    }
    else if (690 <= year && year <= 699)
    {
        stdMonth = 12;
    }
    else
    {
        stdMonth = 10;
    }

    bool isHead = false;
    if ((0 < month && month <= stdMonth) || (year == 700 && month == 12 && isSpec == false))
    {
        isHead = true;
    }
    else if (stdMonth < month && month < 13)
    {
        isHead = false;
        t.Y = year + 1;
    }
    else
    {
        return nullptr;
    }

    int Bd0 = int2(JD::toJD(t)) - J2000;
    SSQ &SSQPtr = SSQ::getInstance();
    std::vector<long double> vecZQ = SSQPtr.getZhongQi();
    if (vecZQ.empty() || Bd0 < vecZQ[0] || Bd0 >= vecZQ[24])
    {
        SSQPtr.calcY(Bd0);
    }

    int leap = SSQPtr.getLeap();
    std::vector<int> vecHS = SSQPtr.getHS();
    std::vector<int> vecYm = SSQPtr.getYm();
    std::vector<int> vecDx = SSQPtr.getDx();
    std::vector<int> vecYueName = SSQPtr.getYueName();

    int jd = 0;
    if (isRun)
    {
        if (leap > 0 && month == vecYueName[leap])
        {
            if (0 < day && day <= vecDx[leap])
            {
                jd = vecHS[leap] + day - 1;
            }
            else
            {
                return nullptr;
            }
        }
    }
    else
    {
        long startIdx = 0;
        if (isHead)
        {
            // 月序第一个月不是正月,说明月序有误
            auto it = std::find(vecYueName.begin(), vecYueName.end(), 1);
            if (it != vecYueName.end())
            {
                startIdx = std::distance(vecYueName.begin(), it);
            }
            else
            {
                // 没有正月，表明月序有误
                return nullptr;
            }
            // 从start开始找月序中有没有当前输入的月
            it = std::find(vecYueName.begin() + startIdx, vecYueName.end(), month);
            if (it != vecYueName.end())
            {
                startIdx = std::distance(vecYueName.begin(), it);
            }
            else
            {
                // 表明输入的月还是有问题
                return nullptr;
            }
            if (isSpec)
            {
                it = std::find(vecYueName.begin() + startIdx + 1, vecYueName.end(), month);
                if (it != vecYueName.end())
                {
                    startIdx = std::distance(vecYueName.begin(), it);
                }
                else
                {
                    // 表明输入的月还是有问题
                    return nullptr;
                }
            }
            if (0 < day && day <= vecDx[startIdx])
            {
                jd = vecHS[startIdx] + day - 1;
            }
            else
            {
                return nullptr;
            }
        }
        else
        {
            // 查正月
            long pos = 0;
            auto it = std::find(vecYueName.begin(), vecYueName.end(), 1);
            if (it != vecYueName.end())
            {
                pos = std::distance(vecYueName.begin(), it);
            }
            else
            {
                // 没有正月，表明月序有误
                return nullptr;
            }
            // 农历700年特殊，头尾各有一个腊月
            it = std::find(vecYueName.begin(), vecYueName.begin() + pos, month);
            if (it != vecYueName.end())
            {
                startIdx = std::distance(vecYueName.begin(), it);
            }
            else
            {
                // 表明输入的月还是有问题
                return nullptr;
            }
            if (isSpec && year != 700)
            {
                it = std::find(vecYueName.begin() + startIdx + 1, vecYueName.end(), month);
                if (it != vecYueName.end())
                {
                    startIdx = std::distance(vecYueName.begin(), it);
                }
                else
                {
                    // 表明输入的月还是有问题
                    return nullptr;
                }
            }
            if (0 < day && day <= vecDx[startIdx])
            {
                jd = vecHS[startIdx] + day - 1;
            }
            else
            {
                return nullptr;
            }
        }
    }

    return new Day(jd);
}

Day *Day::fromLunar(int year, uint8_t month, int day, bool isRun, bool isSpec)
{
    return getLunarDate(year,month,day,isRun,isSpec);
}

Day *Day::fromMoslem(int _year, uint8_t _month, int _day)
{
    Time t{};
    t.h = 12, t.m = 0, t.s = 0.1;
    t.Y = _year;
    t.M = _month;
    t.D = _day;
    auto d0 = int2(JD::toJD(t)) - J2000;

    return new Day(d0);
}
