#include <cmath>
#include <sstream>
#include <iomanip>
#include "bazi.h"
#include "SSQ.h"
#include "eph.h"
#include "lunar_ob.h"
#include "sx_lang_zh.h"

BaziBase::BaziBase(const SBaziInputPara &input)
{
    mGender_ = input.gender;
    mUserName_ = input.name;
    mSUserJw_ = input.jw;
    mLifa_ = input.lifa;
    if (mLifa_ == YuWuWeiZiPingLifa_DingDongZhi || mLifa_ == YuWuWeiZiPingLifa_DingXiaZhi)
    {
        mBzCalcMode_ = PingQiPaiPanMode;
    }
    else
    {
        mBzCalcMode_ = DingQiPaiPanMode;
    }

    mIsAst_ = input.isAst;
    Time bzTime{};
    if (CalendarLunar == input.calendar)
    {
        // 农历转公历
        auto lyYear = input.birthDayTime.getYear();
        auto lyMonth = u_int8_t(input.birthDayTime.getMonth());
        auto lyDay = input.birthDayTime.getDay();
        auto hour = input.birthDayTime.getHour();
        auto min = input.birthDayTime.getMin();
        auto sec = input.birthDayTime.getSec();

        Day *dayPtr = sxtwl::fromLunar(lyYear, lyMonth, lyDay, input.isRun, input.isSpec);
        mSolarYear_ = dayPtr->getSolarYear();
        mSolarMonth_ = dayPtr->getSolarMonth();
        mSolarDay_ = dayPtr->getSolarDay();
        delete dayPtr;
        bzTime = {mSolarYear_, mSolarMonth_, mSolarDay_, hour, min, sec};
    }
    else
    {
        bzTime = input.birthDayTime;
        mSolarYear_ = bzTime.getYear();
        mSolarMonth_ = bzTime.getMonth();
        mSolarDay_ = bzTime.getDay();
    }

    mBdJd_ = JD::toJD(bzTime);
    Time bzAstTime = bzTime;
    if (input.isAst)
    {
        auto res = JD::calcAST(bzTime, mSUserJw_.J);
        mAstJd_ = std::get<0>(res);
        bzAstTime = std::get<1>(res);
    }
    else
    {
        mAstJd_ = mBdJd_;
    }

    mSLunarDay_ = sxtwl::solar2Lunar2(bzAstTime);

    mHeadJieQiJd_ = J2000;
    mTailJieQiJd_ = J2000;
    mHeadJieQiAstJd_ = J2000;
    mTailJieQiAstJd_ = J2000;
    mJyJd_ = J2000;
    arrayQiYunDelta_ = {0, 0, 0};
    arraySiZhu_ = {-1, -1, -1, -1, -1, -1, -1, -1};
    std::array<int, 10> arrayShiShen_ = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    for (auto i = 0; i < 10; i++)
    {
        vecShiShen_.push_back(-1);
    }
    mJyYear_ = 0;
}

BaziBase::~BaziBase()
{
}

void BaziBase::calcBaziPaiPan()
{
    switch (mBzCalcMode_)
    {
    case PingQiPaiPanMode:
        calcPingQiPaiPan();
        break;
    case DingQiPaiPanMode:
        calcDingQiPaiPan();
        break;
    default:
        calcPingQiPaiPan();
        break;
    }
}

std::array<long double, 4> BaziBase::bkCalc(int year, int month, LiFaType lifa)
{
	auto monthTemp = month, yearTemp = year;
	//定冬至
	if (lifa == YuWuWeiZiPingLifa_DingDongZhi)
	{
		monthTemp = 12;
		yearTemp = year - 1;
	}
	//定夏至
	else if (lifa == YuWuWeiZiPingLifa_DingXiaZhi && month < 7)
	{
		yearTemp = year - 1;
		monthTemp = 6;
	}
	else if (lifa == YuWuWeiZiPingLifa_DingXiaZhi && month >= 7)
	{
		yearTemp = year;
		monthTemp = 6;
	}

    Time t1 = {yearTemp,monthTemp,21,12,0,0};
	//auto jdt1 = JD::gcal2jd(yearTemp, monthTemp, 21);
	//auto stjd = std::get<0>(jdt1) + std::get<1>(jdt1) + 0.5 - J2000;
    double stjd = JD::toJD(t1) - J2000;
	//auto jdt2 = JD::gcal2jd(yearTemp+1,monthTemp, 21);
	//auto mdjd = std::get<0>(jdt2) + std::get<1>(jdt2) + 0.5 - J2000;
    Time t2 = {yearTemp+1,monthTemp,21,12,0,0};
    double mdjd = JD::toJD(t2) - J2000;
	//auto jdt3 = JD::gcal2jd(yearTemp+2,monthTemp, 21);
	//auto edjd = std::get<0>(jdt3) + std::get<1>(jdt3) + 0.5 - J2000;
    Time t3 = {yearTemp+2,monthTemp, 21,12,0,0};
    double edjd = JD::toJD(t3) - J2000;

	//前一年中气
	auto pzq = qi_accurate2(stjd, false, 120) + J2000;
	auto zq = qi_accurate2(mdjd, false, 120) + J2000;
	//下一年中气
	auto nzq = qi_accurate2(edjd, false, 120) + J2000;
	std::array<long double,4> z{};
	z[0] = (zq - pzq) / 12.0; //k，当年平气
	z[1] = pzq + z[0] / 2.0;  //b，当年第一个节气小寒
	z[2] = (nzq - zq) / 12.0; //k1，下年平气
	z[3] = zq + z[2]/2.0;     //b1，下一年第一个节气小寒

	return z;
}

// 计算立春，确保此函数中没有修改类的成员函数
long double BaziBase::calcLichun(int year) const
{
    Time t0{};
    t0.Y = year;
    t0.M = 2;
    t0.D = 2;
    t0.h = 12;
    t0.m = 0;
    t0.s = 0;
    auto jd = JD::toJD(t0) - J2000;
    auto lichun = qi_accurate2(jd, false, 120) + J2000;
    if (jd < lichun)
    {
        lichun -= 365.2422;
    }

    return lichun;
}

std::tuple<int, double, double> BaziBase::calcJieQiTermsByPingQi()
{
    int index = 0;
    long double D = 0.0;
    auto year = mSolarYear_;
    auto month = mSolarMonth_;

    auto array = bkCalc(year, month, mLifa_);
    auto k = array[0], b = array[1], k1 = array[2], b1 = array[3];
    auto daxue = b - k;
    auto lichun = b + k;
    if (mBdJd_ < b)
    {
        index = 0;
        mHeadJieQiJd_ = b - k; // 大雪
        mTailJieQiJd_ = b;     // 小寒
    }
    else if (b <= mBdJd_ && mBdJd_ < b1)
    {
        auto n = std::floor((mBdJd_ - b) / k);
        D = k * n + b;
        index = int(n) + 1;
        mHeadJieQiJd_ = D;
        mTailJieQiJd_ = D + k;
    }
    else
    {
        auto n = std::floor((mBdJd_ - b1) / k1);
        D = k1 * n + b1;
        mHeadJieQiJd_ = b1;
        mTailJieQiJd_ = b1 + k1;
        index = 13;
    }

    if (mIsAst_)
    {
        mHeadJieQiAstJd_ = JD::calcAST(mHeadJieQiJd_, mSUserJw_.J);
        mTailJieQiAstJd_ = JD::calcAST(mTailJieQiJd_, mSUserJw_.J);
    }
    else
    {
        mHeadJieQiAstJd_ = mHeadJieQiJd_;
        mTailJieQiAstJd_ = mTailJieQiJd_;
    }

    return std::make_tuple(index, daxue, lichun);
}

void BaziBase::calcJiaoYunDate(bool flag)
{
    long double bd, hJq, tJq, dt1, dt2;
    if (mIsAst_)
    {
        bd = mAstJd_;
        hJq = mHeadJieQiAstJd_;
        tJq = mTailJieQiAstJd_;
    }
    else
    {
        bd = mBdJd_;
        hJq = mHeadJieQiJd_;
        tJq = mTailJieQiJd_;
    }

    if (flag)
    {
        dt1 = hJq;
        dt2 = bd;
    }
    else
    {
        dt1 = bd;
        dt2 = tJq;
    }

    auto offset = dt2 - dt1;
    auto deltaY = std::floor(offset / (365.2422 / 120.0));
    offset = (offset / (365.2422 / 120.0) - deltaY) * 365.2422;
    mJyJd_ = bd + 365.2422 * deltaY + offset;
    Time t = JD::JD2DD(mJyJd_);

    auto sLunar = sxtwl::solar2Lunar2(t);
    int mJyYeaer = sLunar.year;

    auto deltaM = std::floor(offset / (365.2422 / 12.0));
    offset -= deltaM * (365.2422 / 12.0);
    auto deltaD = std::floor(offset);
    arrayQiYunDelta_ = {int(deltaY), int(deltaM), int(deltaD)};
}

void BaziBase::calcJieQiTermsByDingQi()
{
    auto year = mSolarYear_;
    auto month = mSolarMonth_;
    auto day = mSolarDay_;
    auto jd = mBdJd_;
    auto jqjd = getJieByDingQi(year, month);
    if (jqjd > jd)
    {
        long double jqjdPre = jd;
        while (jd <= jqjdPre)
        {
            month--;
            if (month < 1)
            {
                year--;
                month = 12; 
            }
            jqjdPre = getJieByDingQi(year, month);
        }
        mHeadJieQiJd_ = jqjdPre;
        mTailJieQiJd_ = jqjd;
    }
    else
    {
        long double jqjdNext = jd;
        while (jd >= jqjdNext)
        {
            month++;
            if (month > 12)
            {
                year++;
                month = 1;
            }
            jqjdNext = getJieByDingQi(year, month);
        }

        mHeadJieQiJd_ = jqjd;
        mTailJieQiJd_ = jqjdNext; 
    }

    if (mIsAst_)
    {
        mHeadJieQiAstJd_ = JD::calcAST(mHeadJieQiJd_, mSUserJw_.J);
        mTailJieQiAstJd_ = JD::calcAST(mTailJieQiJd_, mSUserJw_.J);
    }
    else
    {
        mHeadJieQiAstJd_ = mHeadJieQiJd_;
        mTailJieQiAstJd_ = mTailJieQiJd_;
    }
}

long double BaziBase::getJieByDingQi(int year, int month)
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
    bool isJie = false;
    for (int i = 0; i < Bdn; i++)
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

        if (d0 == vecZQ[qk] && (qk & 1) == 1)
        {
            isJie = true;
            break;
        }
    }
    if (isJie)
    {
        return qi_accurate2(d0, false, 120) + J2000;
    }
    else
    {
        month++;
        if (month > 12)
        {
            year++;
            month = 1;
        }
        return getJieByDingQi(year, month);
    }
}

void BaziBase::calcPingQiPaiPan()
{
    auto jd = mBdJd_;
    if (mIsAst_)
    {
        jd = mAstJd_;
    }

    auto res = calcJieQiTermsByPingQi();
    auto index = std::get<0>(res);
    auto daxue = std::get<1>(res), lichun = std::get<2>(res);
    if (mBdJd_ < lichun)
    {
        lichun -= 365;
    }
    // 1984年立春作为起点
    lichun -= 2445735;
    auto nianZhu = int2(lichun / 365.2422 + 0.5) + 9000;

    // #1998 12 7日 大雪作为起点
    daxue -= 2451155;
    auto yueZhu = int2(daxue / 365.2422 + 0.4) * 12 + index + 90000;

    // #2000年1月7日作为起点
    jd += 13 / 24.0 - J2000; // 转为前一日23点起算(原jd为本日中午12点起算)
    auto D = int2(jd);
    auto SC = int2((jd - D) * 12); // 日数与时辰
    auto riZhu = D - 6 + 9000000;
    auto shiZhu = (D - 1) * 12 + 90000000 + SC;

    arraySiZhu_[0] = nianZhu % 10;
    arraySiZhu_[1] = nianZhu % 12;
    arraySiZhu_[2] = yueZhu % 10;
    arraySiZhu_[3] = yueZhu % 12;
    arraySiZhu_[4] = riZhu % 10;
    arraySiZhu_[5] = riZhu % 12;
    arraySiZhu_[6] = shiZhu % 10;
    arraySiZhu_[7] = shiZhu % 12;

    auto j = arraySiZhu_[4];
    if (j % 2 == 0)
    {
        for (auto i = 0; i < 10; i++)
        {
            vecShiShen_[j] = i;
            j += 1;
            j %= 10;
        }
    }
    else
    {
        for (auto k = 0; k < 9; k += 2)
        {
            vecShiShen_[j] = k;
            vecShiShen_[j - 1] = k + 1;
            j += 2;
            j %= 10;
        }
    }
    vecShiShen_.push_back(10);
    auto isFemale = 0;
    if (mGender_)
    {
        isFemale = 1; // 女性
    }

    bool flag = false;
    if ((arraySiZhu_[0] % 2) ^ ((isFemale) == 1))
    {
        flag = true;
    }
    calcJiaoYunDate(flag);
}

void BaziBase::calcDingQiPaiPan()
{
    auto jd = mBdJd_;
    if (mIsAst_)
    {
        jd = mAstJd_;
    }

    auto jd1 = mBdJd_ + (-8.0) / 24 - J2000;
    auto jd2 = jd1 + dt_T(jd1);                         // 力学时
    auto w = XL::S_aLon(jd2 / 36525.0, -1);             // 此刻太阳视黄经
    int k = int2((w / pi2 * 360 + 45 + 15 * 360) / 30); // 1984年立春起算的节气数(不含中气)
    // 年
    int v = int2(k / 12.0 + 6000000);
    arraySiZhu_[0] = v % 10;
    arraySiZhu_[1] = v % 12;
    // 月
    v = k + 2 + 60000000;
    arraySiZhu_[2] = v % 10;
    arraySiZhu_[3] = v % 12;

    // 子时转日
    jd += 13.0 / 24 - J2000 ;                            // 转为前一日23点起算(原jd为本日中午12点起算)
    int D = int2(jd), SC = int2((jd - D) * 12); // 日数与时辰
    // 日
    v = D - 6 + 9000000;
    arraySiZhu_[4] = v % 10;
    arraySiZhu_[5] = v % 12;
    // 时
    v = (D - 1) * 12 + 90000000 + SC;
    arraySiZhu_[6] = v % 10;
    arraySiZhu_[7] = v % 12;
    // 十神
    auto j = arraySiZhu_[4];
    if (j % 2 == 0)
    {
        for (auto i = 0; i < 10; i++)
        {
            vecShiShen_[j] = i;
            j += 1;
            j %= 10;
        }
    }
    else
    {
        for (auto p = 0; p < 9; p += 2)
        {
            vecShiShen_[j] = p;
            vecShiShen_[j - 1] = p + 1;
            j += 2;
            j %= 10;
        }
    }
    vecShiShen_.push_back(10);
    auto isFemale = 0;
    if (mGender_)
    {
        isFemale = 1; // 女性
    }

    bool flag = false;
    if ((arraySiZhu_[0] % 2) ^ ((isFemale) == 1))
    {
        flag = true;
    }
    // 计算节气
    calcJieQiTermsByDingQi();
    // 计算交运
    calcJiaoYunDate(flag);
}

std::string BaziBase::getSolarBirth() const
{
    std::string s;
    // 1582-10-15 12:00:00
    if (int(mBdJd_) >= 2299161)
    {
        s = "『公历生日』 ";
    }
    else
    {
        s = "『儒略历生日』 ";
    }
    s += JD::JD2str(mBdJd_);
    return s;
}

std::string BaziBase::getUserName() const
{
    return "『姓名』  " + mUserName_;
}
std::string BaziBase::getUserGender() const
{
    if (mGender_)
    {
        return "『性别』 女";
    }
    else
    {
        return "『性别』 男";
    }
}
std::tuple<std::string, std::string, std::string> BaziBase::getLunarInfo() const
{
    int solarYear = mSolarYear_;

    std::string lySx = "『生肖』 " + std::string(ShengXiao[(solarYear - 1984 + 12000) % 12]);

    // 年号
    LunarYear ob(solarYear);
    std::string lyNh = ob.getNianHao();

    // 干支，以立春为界
    std::string lyGan = Gan[mSLunarDay_.lyGanIdx];
    std::string lyZhi = Zhi[mSLunarDay_.lyZhiIdx];

    // 月信息
    bool bdLeapYear = false;
    std::string strTemp1, strTemp2, lmMc;
    if (mSLunarDay_.isLeap)
    {
        if (solarYear >= -721 && solarYear < -220)
        {
            strTemp1 = BDLeapYueName[0];
            bdLeapYear = true;
        }
        else if (solarYear >= -220 && solarYear <= -104)
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
        lmMc = strTemp1;
    }
    else
    {
        if (mSLunarDay_.isNext)
        {
            strTemp2 = SYmc[mSLunarDay_.monthIdx];
        }
        else
        {
            strTemp2 = Ymc[mSLunarDay_.monthIdx];
        }
        lmMc = strTemp1 + strTemp2 + "月";
    }

    // 农历日信息
    std::string ldMc = std::string(Rmc[mSLunarDay_.dayIdx]) + "日";
    std::string retStr = "农历" + lyGan + lyZhi + "年" + lmMc + ldMc;
    return std::make_tuple(lyNh, lySx, retStr);
}

std::string BaziBase::getLifa() const
{
    std::string s;
    if (mLifa_ > LifaUnknown && mLifa_ < XianDaiNongLifa_DingQiFa)
    {
        s = "依据" + std::string(LifaList[mLifa_ - 1]) + "数据拟合";
        return "『定气方式』 " + s;
    }
    else if (XianDaiNongLifa_DingQiFa == mLifa_)
    {
        s = "现代农历定气";
        return "『定气方式』 " + s;
    }
    else
    {
        return "『定气方式』 依据尤氏子平历计算节气交接日期";
    }
}

std::string BaziBase::getAst() const
{
    std::string s = "『出生地』" + mSUserJw_.s + mSUserJw_.x;
    s += ",『经度』" + toFixed(mSUserJw_.J, 6) + ",『纬度』" + toFixed(mSUserJw_.W, 6) + "\n";
    s += "『出生地真太阳时』" + JD::formatStr(mAstJd_) + "\n";
    return s;
}

std::string BaziBase::getShenXiao() const
{
    auto &Bazi = arraySiZhu_;
    auto lySx = std::string(ShengXiao[Bazi[1]]);
    return "『生肖』 " + lySx;
}

std::string BaziBase::getAge() const
{
    std::string age;
    Time nowTime = JD::getNowTime();
    auto nowJd = JD::toJD(nowTime);
    if (nowJd < mBdJd_)
    {
        age = "『年龄』未出生";
        return age;
    }
    else
    {
        auto solarYear = mSolarYear_;
        auto startLichun = calcLichun(solarYear);
        auto endLichun = calcLichun(nowTime.getYear());
        auto iAge = int2((endLichun - startLichun) / 365.2422 + 0.4) + 1;
        if (iAge > 120.0)
        {
            return "『年龄』 历史人物";
        }
        return "『年龄』 虚岁" + std::to_string(iAge);
    }
}

std::string BaziBase::getJieQiIterms() const
{
    long double hJq, tJq;
    std::string astStr;
    if (mIsAst_)
    {
        astStr = "出生地真太阳时";
        hJq = mHeadJieQiAstJd_;
        tJq = mTailJieQiAstJd_;
    }
    else
    {
        astStr = "";
        hJq = mHeadJieQiJd_;
        tJq = mTailJieQiJd_;
    }

    std::string s = "『" + std::string(JieLing[arraySiZhu_[3]]) + "』 " + astStr + JD::formatStr(hJq) + "\n";
    s += "『" + std::string(JieLing[arraySiZhu_[3] + 1]) + "』 " + astStr + JD::formatStr(tJq);
    return s;
}

std::string BaziBase::getJieLing() const
{
    std::string s;
    if (mIsAst_)
    {
        s += getAst();
    }
    s += getJieQiIterms();
    s += "\n\n";
    return s;
}

std::vector<std::vector<std::string>> BaziBase::getSiZhu() const
{
    auto &bazi = arraySiZhu_;
    auto &ssgx = vecShiShen_;
    std::string sGender;
    auto isFemale = 0;
    if (mGender_)
    {
        isFemale = 1; // 女性
        sGender = Gender[isFemale];
    }
    else
    {
        isFemale = 0;
        sGender = Gender[isFemale];
    }
    std::vector<std::string> shiShenArray, tianGanArray, diZhiArray;
    shiShenArray.push_back("");
    shiShenArray.push_back(ShiShen[ssgx[bazi[0]]]);
    shiShenArray.push_back(ShiShen[ssgx[bazi[2]]]);
    shiShenArray.push_back("日元");
    shiShenArray.push_back(ShiShen[ssgx[bazi[6]]]);

    // tiangan
    tianGanArray.push_back(sGender);
    tianGanArray.push_back(Gan[bazi[0]]);
    tianGanArray.push_back(Gan[bazi[2]]);
    tianGanArray.push_back(Gan[bazi[4]]);
    tianGanArray.push_back(Gan[bazi[6]]);

    // dizhi
    diZhiArray.push_back("");
    diZhiArray.push_back(Zhi[bazi[1]]);
    diZhiArray.push_back(Zhi[bazi[3]]);
    diZhiArray.push_back(Zhi[bazi[5]]);
    diZhiArray.push_back(Zhi[bazi[7]]);

    std::vector<std::vector<std::string>> outPut;
    outPut.push_back(shiShenArray);
    outPut.push_back(tianGanArray);
    outPut.push_back(diZhiArray);

    for (size_t i = 0; i < 3; i++)
    {
        std::vector<std::string> cangGan;
        if (0 == i)
        {
            cangGan.push_back("藏干");
        }
        else
        {
            cangGan.push_back("");
        }

        for (size_t j = 1; j < 8; j += 2)
        {
            auto k = gCangGan[bazi[j]][i];
            if (k < 10)
            {
                std::string opt = std::string(Gan[k]) + " " + ShiShen[ssgx[k]];
                cangGan.push_back(opt);
            }
            else
            {
                cangGan.push_back("");
            }
        }

        outPut.push_back(cangGan);
    }

    return outPut;
}

std::string BaziBase::getQiYun() const
{
    std::string s;
    auto &qy = arrayQiYunDelta_;
    s += "起运 命主于出生后" + std::to_string(qy[0]) + "年" + std::to_string(qy[1]) + "个月" + std::to_string(qy[2]) + "天后起运\n";
    return s;
}

std::string BaziBase::getJiaoYun() const
{
    auto &jydt = mJyJd_;
    std::string res = "交运 命主于" + JD::formatStr(jydt) + "交运\n";
    return res;
}

std::vector<std::string> BaziBase::getDaYunList() const
{
    int isFemale = 0;
    if (mGender_)
    {
        isFemale = 1; // 女性
    }
    else
    {
        isFemale = 0;
    }

    auto &Bazi = arraySiZhu_;
    bool flag = (Bazi[0] % 2) ^ (isFemale);
    std::array<int, 2> offsets = {0, 0};
    if (flag)
    {
        offsets = {9, 11};
    }
    else
    {
        offsets = {1, 1};
    }

    auto j = Bazi[2];
    auto k = Bazi[3];
    std::string opt;
    std::vector<std::string> output;
    for (int i = 0; i < 8; i++)
    {
        j += offsets[0];
        j %= 10;
        k += offsets[1];
        k %= 12;

        opt = std::string(Gan[j]) + Zhi[k];
        output.push_back(opt);
    }
    return output;
}

std::vector<int> BaziBase::getStartYearList() const
{
    int jyYear = mJyYear_;

    std::vector<int> output;
    for (int i = 0; i < 8; i++)
    {
        output.push_back(jyYear);
        jyYear += 10;
    }
    return output;
}

std::vector<std::string> BaziBase::getFleetingYearList() const
{

    int jyYear = mJyYear_;

    auto p = (jyYear + 6) % 10;
    auto q = (jyYear + 8) % 12;
    std::vector<std::string> output;
    std::string opt;
    for (int i = 0; i < 10; i++)
    {
        int m = q;
        for (int j = 0; j < 8; j++)
        {
            opt = std::string(Gan[p]) + Zhi[m];
            m += 10;
            m %= 12;
            output.push_back(opt);
        }
        p += 1;
        p %= 10;
        q += 1;
        q %= 12;
    }
    return output;
}

std::vector<int> BaziBase::getEndYearList() const
{
    int jyYear = mJyYear_;
    auto endYear = jyYear + 9;
    std::vector<int> output;
    for (int i = 0; i < 8; i++)
    {
        output.push_back(endYear);
        endYear += 10;
    }
    return output;
}

std::string BaziBase::printBazi() const
{
    // 姓名，性别，生肖，年龄
    std::string s;
    auto res = getLunarInfo();
    s += getUserName() + getUserGender() + getShenXiao() + getAge() + "\n";
    // 公历生日
    s += getSolarBirth() + "\n";
    // 出生年代，农历生日

    s += "『出生年代』 " + std::get<0>(res) + "\n";
    s += "『农历生日』 " + std::get<2>(res) + "\n";
    // 定气方式
    s += getLifa() + "\n";

    // 节气
    s += getJieLing();

    // 四柱藏干盘面排列
    auto &bazi = arraySiZhu_;
    auto &ssgx = vecShiShen_;

    int isFemale = 0;
    std::string gender;
    if (mGender_)
    {
        isFemale = 1; // 女性
        gender = Gender[isFemale];
    }
    else
    {
        isFemale = 0;
        gender = Gender[isFemale];
    }

    std::ostringstream oss;
    oss << std::left << std::setw(11) << " " << std::setw(16) << ShiShen[ssgx[bazi[0]]] << std::setw(16) << ShiShen[ssgx[bazi[2]]] << std::setw(16) << "日元" << std::setw(17) << ShiShen[ssgx[bazi[6]]] << std::endl;
    s += "\033[31m" + oss.str() + "\033[0m";
    oss.str("");

    oss << std::left << std::setw(14) << gender << std::setw(15) << Gan[bazi[0]] << std::setw(15) << Gan[bazi[2]] << std::setw(15) << Gan[bazi[4]] << std::setw(15) << Gan[bazi[6]] << std::endl;
    s += oss.str();
    oss.str("");

    oss << std::left << std::setw(12) << " " << std::setw(15) << Zhi[bazi[1]] << std::setw(15) << Zhi[bazi[3]] << std::setw(15) << Zhi[bazi[5]] << std::setw(15) << Zhi[bazi[7]] << std::endl;
    s += oss.str() + "\n";

    for (int i = 0; i < 3; i++)
    {
        std::string opt;
        for (int j = 1; j < 8; j += 2)
        {
            auto k = gCangGan[bazi[j]][i];
            if (k < 10)
            {
                opt += std::string(Gan[k]) + " " + ShiShen[ssgx[k]] + std::string(7, ' ');
            }
            else
            {
                opt += std::string(14, ' ');
            }
        }
        if (0 == i)
        {
            s += "\033[34m" + std::string("藏干") + std::string(8, ' ') + opt + "\n" + "\033[0m";
        }
        else
        {
            s += "\033[34m" + std::string(12, ' ') + opt + "\n" + "\033[0m";
        }
    }
    s += "\n";
    s += getQiYun();
    s += getJiaoYun();
    s += "\n";
    return s;
}
