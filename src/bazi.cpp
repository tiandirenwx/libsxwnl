#include <cmath>
#include <sstream>
#include <iomanip>
#include "bazi.h"
#include "bazi_analysis.h"
#include "SSQ.h"
#include "eph.h"
#include "lunar_ob.h"
#include "sx_lang_zh.h"
#include "lunar_month_name.h"

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
        auto lyMonth = uint8_t(input.birthDayTime.getMonth());
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
    mMidQiJd_ = J2000;
    mMidQiAstJd_ = J2000;
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
    mPaiPanDone_ = true;
}

BaziBase::PingQiFrame BaziBase::bkCalc(int year, int month, LiFaType lifa)
{
	int yearTemp = year;
	// 定冬至: 以冬至为气历起点, 公历 1–11 月归属上一年气历
	if (lifa == YuWuWeiZiPingLifa_DingDongZhi)
	{
		yearTemp = year - 1;
	}
	// 定夏至: 以夏至为气历起点, 公历 1–6 月归属上一年气历
	else if (lifa == YuWuWeiZiPingLifa_DingXiaZhi && month < 7)
	{
		yearTemp = year - 1;
	}

	// 使用SSQ类精确计算中气。
	// 注意: calcY 会把 SSQ 单例的年份窗口挪走, 本函数不负责还原 —— 读 SSQ 的一方
	// (Day::checkSSQ 等) 都会自行校验窗口是否覆盖目标日, 不覆盖就重算一次。
	SSQ &ssq = SSQ::getInstance();

	// ppzq: 上上个至点(仅用于求上一气年的间隔 kPrev)
	long double ppzq, pzq, zq, nzq;

	if (lifa == YuWuWeiZiPingLifa_DingXiaZhi)
	{
		// 定夏至: 夏至位于 vecZQ 索引 12, 一次 calcY 只能取到一个夏至,
		// 故 yearTemp-1 / yearTemp / yearTemp+1 / yearTemp+2 各算一次。
		auto xiaZhiOf = [&ssq](int y) -> long double {
			ssq.calcY(int2((y - 2000.0) * 365.2422 + 180));
			return qi_accurate2(ssq.getZhongQi()[12], false, 120) + J2000;
		};
		ppzq = xiaZhiOf(yearTemp - 1);
		pzq  = xiaZhiOf(yearTemp);
		zq   = xiaZhiOf(yearTemp + 1);
		nzq  = xiaZhiOf(yearTemp + 2);
	}
	else
	{
		// 定冬至: 一次 calcY 的 vecZQ[0]/vecZQ[24] 就是相邻两个冬至, 三次即可覆盖四个。
		int baseJD0 = int2((yearTemp - 2000.0) * 365.2422 + 180);
		ssq.calcY(baseJD0);
		std::vector<long double> vecZQ0 = ssq.getZhongQi();
		ppzq = qi_accurate2(vecZQ0[0], false, 120) + J2000;  // yearTemp-1 年冬至

		int baseJD1 = int2((yearTemp + 1 - 2000.0) * 365.2422 + 180);
		ssq.calcY(baseJD1);
		std::vector<long double> vecZQ1 = ssq.getZhongQi();
		pzq = qi_accurate2(vecZQ1[0], false, 120) + J2000;  // yearTemp年冬至
		zq = qi_accurate2(vecZQ1[24], false, 120) + J2000;  // yearTemp+1年冬至

		// 计算yearTemp+2年的节气数据，vecZQ[24]是yearTemp+2年冬至
		int baseJD2 = int2((yearTemp + 2 - 2000.0) * 365.2422 + 180);
		ssq.calcY(baseJD2);
		std::vector<long double> vecZQ2 = ssq.getZhongQi();
		nzq = qi_accurate2(vecZQ2[24], false, 120) + J2000; // yearTemp+2年冬至
	}

	PingQiFrame f{};
	f.pzq   = pzq;
	f.kPrev = (pzq - ppzq) / 12.0; // 上一气年的节间隔
	f.k     = (zq - pzq) / 12.0;   // 本气年的节间隔
	f.b     = pzq + f.k / 2.0;     // 本气年首节: 定冬至=小寒, 定夏至=小暑
	f.k1    = (nzq - zq) / 12.0;   // 下一气年的节间隔
	f.b1    = zq + f.k1 / 2.0;     // 下一气年首节
	return f;
}

// 计算立春，确保此函数中没有修改类的成员函数
long double BaziBase::calcLichun(int year) const
{
	// 使用SSQ类精确计算立春
	SSQ &ssq = SSQ::getInstance();
	
	// 计算指定年份的节气信息，以该年年中为基准
	int baseJD = int2((year - 2000.0) * 365.2422 + 180);
	ssq.calcY(baseJD);
	std::vector<long double> vecZQ = ssq.getZhongQi();
	
	// 立春在vecZQ数组中的索引是3
	auto lichun = qi_accurate2(vecZQ[3], false, 120) + J2000;
	
	// 验证计算结果的年份，如果不是目标年份则使用上一年的立春
	Time t = JD::JD2DD(lichun);
	if (t.Y != year)
	{
		// 重新计算上一年的立春
		baseJD = int2((year - 1 - 2000.0) * 365.2422 + 180);
		ssq.calcY(baseJD);
		vecZQ = ssq.getZhongQi();
		lichun = qi_accurate2(vecZQ[3], false, 120) + J2000;
	}
	
	return lichun;
}

std::tuple<int, double, double> BaziBase::calcJieQiTermsByPingQi()
{
    int index = 0;
    auto year = mSolarYear_;
    auto month = mSolarMonth_;

    auto f = bkCalc(year, month, mLifa_);
    const auto kPrev = f.kPrev, k = f.k, b = f.b, k1 = f.k1, b1 = f.b1, pzq = f.pzq;

    // index 的零点: 上一气年的最后一个节 —— 定冬至=大雪, 定夏至=芒种。
    // 它位于 [ppzq, pzq] 这个气年里, 因此必须用 kPrev; 若借用本气年的 k,
    // 会因相邻至点间隔逐年起伏(十几分钟量级)而使该节偏出几十秒,
    // 造成"同一个节从 12 月侧和 1 月侧查到两个值"。
    const long double baseNode = pzq - kPrev / 2.0;

    if (mBdJd_ < b)
    {
        // 出生早于本气年首节 => 所在节属于上一气年, 整段用 kPrev 定位。
        // 实际只会落在 baseNode 这一个节上 (定冬至=公历 1 月初小寒前,
        // 定夏至=公历 7 月初小暑前), 这里仍按通式处理, 不对历法漂移做假设。
        auto steps = std::ceil((baseNode - mBdJd_) / kPrev);
        if (steps < 0) steps = 0;
        index = -int(steps);
        mHeadJieQiJd_ = baseNode - steps * kPrev;
        mTailJieQiJd_ = (steps == 0) ? b : (mHeadJieQiJd_ + kPrev);
        mMidQiJd_ = mHeadJieQiJd_ + kPrev / 2.0;
    }
    else if (mBdJd_ < b1)
    {
        auto n = std::floor((mBdJd_ - b) / k);
        // 由于 k(当年) 与 k1(下年) 因日心运动不匀而有微差,
        // 会产生 [b+12k, b1) 的一段"隙缝"(通常十几秒到几分钟),
        // 若 mBdJd 落在隙缝内, floor 可能得到 n=12, 使 index=13,
        // 与下面 ELSE 分支的 index=13 冲突, 导致月柱无法切换.
        // 此处把 n 限制在 [0,11]: 位于隙缝内视作 "第 12 个节 (index=12)
        // 已开始, 而月末应对齐下年首节 b1".
        if (n >= 12) n = 11;
        index = int(n) + 1;
        mHeadJieQiJd_ = b + k * n;
        // 若当前为最后一个节 (index=12), 月末节即下年首节 b1
        // (而非 b + 12k, 后者不精确, 会与 ELSE 分支的 head 存在数秒空隙).
        mTailJieQiJd_ = (index == 12) ? b1 : (mHeadJieQiJd_ + k);
        mMidQiJd_ = mHeadJieQiJd_ + k / 2.0;
    }
    else
    {
        auto n = std::floor((mBdJd_ - b1) / k1);
        index = 13 + int(n);
        mHeadJieQiJd_ = b1 + k1 * n;
        mTailJieQiJd_ = mHeadJieQiJd_ + k1;
        mMidQiJd_ = mHeadJieQiJd_ + k1 / 2.0;
    }

    // 立春在气年内的相位(以至点为 0):
    //   定冬至 冬至→小寒→大寒→立春, 立春 = pzq + 1.5k
    //   定夏至 夏至起第 15 个气即立春,  立春 = pzq + 7.5k
    // 用 k 而非 365.2422/12 等常数, 保证年柱与月柱在立春同一瞬间切换。
    const long double lcPhase = (mLifa_ == YuWuWeiZiPingLifa_DingXiaZhi) ? 7.5L : 1.5L;
    long double lichun = pzq + lcPhase * k;
    if (mBdJd_ < lichun)
    {
        // 出生在本气年立春之前, 年柱归上一气年; 该立春同样落在上一气年, 用 kPrev。
        lichun = pzq - (12.0L - lcPhase) * kPrev;
    }

    if (mIsAst_)
    {
        mHeadJieQiAstJd_ = JD::calcAST(mHeadJieQiJd_, mSUserJw_.J);
        mTailJieQiAstJd_ = JD::calcAST(mTailJieQiJd_, mSUserJw_.J);
        mMidQiAstJd_ = JD::calcAST(mMidQiJd_, mSUserJw_.J);
    }
    else
    {
        mHeadJieQiAstJd_ = mHeadJieQiJd_;
        mTailJieQiAstJd_ = mTailJieQiJd_;
        mMidQiAstJd_ = mMidQiJd_;
    }

    // 定夏至以芒种(午月)为 index 零点, 定冬至以大雪(子月)为零点, 两者相差 6 个月。
    // 下游月柱公式按子月基准写死, 故在此补偿。
    if (mLifa_ == YuWuWeiZiPingLifa_DingXiaZhi)
    {
        index += 6;
    }

    return std::make_tuple(index, (double)baseNode, (double)lichun);
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
    mJyYear_ = sLunar.year;

    auto deltaM = std::floor(offset / (365.2422 / 12.0));
    offset -= deltaM * (365.2422 / 12.0);
    auto deltaD = std::floor(offset);
    arrayQiYunDelta_ = {int(deltaY), int(deltaM), int(deltaD)};
}

void BaziBase::calcJieQiTermsByDingQi()
{
    // 直接使用 SSQ vecZQ (25 项, 从当期 冬至 到 下期 冬至),
    // 从中找出与 mBdJd_ 相邻的两个 "节" (奇数索引).
    //
    // 原实现用 getJieByDingQi(year, month) 逐月扫描, 只能返回该
    // 公历月中"第一个"节. 当 Julian/公历 月里恰好含有两个节时
    // (远古或 Julian 年份, 如 900-01 Julian 同时含 小寒 与 立春),
    // 该实现会错过第二个节, 造成 head 与真正的月柱不一致.
    long double jdOff = mBdJd_ - J2000;   // vecZQ 存储的是 J2000 偏移量

    SSQ &ssq = SSQ::getInstance();
    std::vector<long double> vecZQ = ssq.getZhongQi();
    // 若当前 vecZQ 未覆盖 jdOff, 重新计算
    if (vecZQ.empty() || jdOff < vecZQ[0] || jdOff >= vecZQ[24])
    {
        ssq.calcY(int2(jdOff));
        vecZQ = ssq.getZhongQi();
    }

    // vecZQ 内奇数索引 (1,3,5,...,23) 为 12 节, 偶数索引为 12 中气.
    // 对每个节做 qi_accurate2 精化, 得到绝对 JD (加 J2000).
    // 月首节取到索引 i 后, 紧随其后的 vecZQ[i+1] 即本月中气 (命宫论月用).
    long double head = 0, tail = 0, mid = 0;
    bool haveHead = false, haveTail = false;
    for (int i = 1; i <= 23; i += 2)
    {
        long double jNode = qi_accurate2(vecZQ[i], false, 120) + J2000;
        if (jNode <= mBdJd_)
        {
            head = jNode;
            mid = qi_accurate2(vecZQ[i + 1], false, 120) + J2000;
            haveHead = true;
        }
        else
        {
            tail = jNode;
            haveTail = true;
            break;
        }
    }
    // mBdJd_ 早于本期首个节: 计算上一期
    if (!haveHead)
    {
        // 位移 300 天回退, 覆盖前一个 冬至→冬至 周期
        ssq.calcY(int2(jdOff) - 300);
        auto vecZQPrev = ssq.getZhongQi();
        for (int i = 23; i >= 1; i -= 2)
        {
            long double jNode = qi_accurate2(vecZQPrev[i], false, 120) + J2000;
            if (jNode <= mBdJd_)
            {
                head = jNode;
                mid = qi_accurate2(vecZQPrev[i + 1], false, 120) + J2000;
                haveHead = true;
                break;
            }
        }
    }
    // mBdJd_ 晚于本期末尾节: 计算下一期
    if (!haveTail)
    {
        ssq.calcY(int2(jdOff) + 300);
        auto vecZQNext = ssq.getZhongQi();
        for (int i = 1; i <= 23; i += 2)
        {
            long double jNode = qi_accurate2(vecZQNext[i], false, 120) + J2000;
            if (jNode > mBdJd_)
            {
                tail = jNode;
                haveTail = true;
                break;
            }
        }
    }

    mHeadJieQiJd_ = head;
    mTailJieQiJd_ = tail;
    mMidQiJd_ = mid;

    if (mIsAst_)
    {
        mHeadJieQiAstJd_ = JD::calcAST(mHeadJieQiJd_, mSUserJw_.J);
        mTailJieQiAstJd_ = JD::calcAST(mTailJieQiJd_, mSUserJw_.J);
        mMidQiAstJd_ = JD::calcAST(mMidQiJd_, mSUserJw_.J);
    }
    else
    {
        mHeadJieQiAstJd_ = mHeadJieQiJd_;
        mTailJieQiAstJd_ = mTailJieQiJd_;
        mMidQiAstJd_ = mMidQiJd_;
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
    // daxue: index 的零点(定冬至=大雪, 定夏至=芒种); lichun: 已保证 <= 出生时刻
    auto daxue = std::get<1>(res), lichun = std::get<2>(res);
    // 1984年立春作为起点
    lichun -= 2445735;
    auto nianZhu = int2(lichun / 365.2422 + 0.5) + 9000;

    // #1998 12 7日 大雪作为起点
    daxue -= 2451155;
    auto yueZhu = int2(daxue / 365.2422 + 0.4) * 12 + index + 90000;

    // #2000年1月7日作为起点
    jd += 13 / 24.0 - J2000; // 转为前一日23点起算(原jd为本日中午12点起算)
    auto D = int2(jd);
    // 加微小 epsilon 抵消浮点误差: 时辰边界(整点)处 (jd-D)*12 可能因浮点下溢
    // 略小于整数(如 3.99999999627), 导致 floor 后归入前一个时辰(戊辰误判为丁卯)。
    // jd 量级约 2.4e6, long double 量化使此处误差可达 ~3e-9; 取 1e-7 既覆盖该误差,
    // 又远小于 1 秒(SC 空间约 1.39e-4), 不会误伤真实非边界时间。
    auto SC = int2((jd - D) * 12 + 1e-7); // 日数与时辰
    auto riZhu = D - 6 + 9000000;
    auto shiZhu = (D - 1) * 12 + 90000000 + SC;

    arraySiZhu_[0] = nianZhu % 10;
    if (arraySiZhu_[0] < 0) arraySiZhu_[0] += 10;
    arraySiZhu_[1] = nianZhu % 12;
    if (arraySiZhu_[1] < 0) arraySiZhu_[1] += 12;
    arraySiZhu_[2] = yueZhu % 10;
    if (arraySiZhu_[2] < 0) arraySiZhu_[2] += 10;
    arraySiZhu_[3] = yueZhu % 12;
    if (arraySiZhu_[3] < 0) arraySiZhu_[3] += 12;
    arraySiZhu_[4] = riZhu % 10;
    if (arraySiZhu_[4] < 0) arraySiZhu_[4] += 10;
    arraySiZhu_[5] = riZhu % 12;
    if (arraySiZhu_[5] < 0) arraySiZhu_[5] += 12;
    arraySiZhu_[6] = shiZhu % 10;
    if (arraySiZhu_[6] < 0) arraySiZhu_[6] += 10;
    arraySiZhu_[7] = shiZhu % 12;
    if (arraySiZhu_[7] < 0) arraySiZhu_[7] += 12;

    int ssMap[10];
    bazi::buildShiShenMap(arraySiZhu_[4], ssMap);
    for (int i = 0; i < 10; i++)
        vecShiShen_[i] = ssMap[i];
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
    if (arraySiZhu_[0] < 0) arraySiZhu_[0] += 10;
    arraySiZhu_[1] = v % 12;
    if (arraySiZhu_[1] < 0) arraySiZhu_[1] += 12;
    // 月
    v = k + 2 + 60000000;
    arraySiZhu_[2] = v % 10;
    if (arraySiZhu_[2] < 0) arraySiZhu_[2] += 10;
    arraySiZhu_[3] = v % 12;
    if (arraySiZhu_[3] < 0) arraySiZhu_[3] += 12;

    // 子时转日
    jd += 13.0 / 24 - J2000 ;                            // 转为前一日23点起算(原jd为本日中午12点起算)
    int D = int2(jd);
    // 加微小 epsilon 抵消浮点误差: 时辰边界(整点)处 (jd-D)*12 可能因浮点下溢
    // 略小于整数(如 3.99999999627), 导致 floor 后归入前一个时辰(戊辰误判为丁卯)。
    // jd 量级约 2.4e6, long double 量化使此处误差可达 ~3e-9; 取 1e-7 既覆盖该误差,
    // 又远小于 1 秒(SC 空间约 1.39e-4), 不会误伤真实非边界时间。
    int SC = int2((jd - D) * 12 + 1e-7); // 日数与时辰
    // 日
    v = D - 6 + 9000000;
    arraySiZhu_[4] = v % 10;
    if (arraySiZhu_[4] < 0) arraySiZhu_[4] += 10;
    arraySiZhu_[5] = v % 12;
    if (arraySiZhu_[5] < 0) arraySiZhu_[5] += 12;
    // 时
    v = (D - 1) * 12 + 90000000 + SC;
    arraySiZhu_[6] = v % 10;
    if (arraySiZhu_[6] < 0) arraySiZhu_[6] += 10;
    arraySiZhu_[7] = v % 12;
    if (arraySiZhu_[7] < 0) arraySiZhu_[7] += 12;
    int ssMap[10];
    bazi::buildShiShenMap(arraySiZhu_[4], ssMap);
    for (int i = 0; i < 10; i++)
        vecShiShen_[i] = ssMap[i];
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

int BaziBase::getDaysAfterJie() const
{
    long double bd = mIsAst_ ? mAstJd_ : mBdJd_;
    long double jie = mIsAst_ ? mHeadJieQiAstJd_ : mHeadJieQiJd_;
    int days = int(std::floor(bd - jie));
    return days < 0 ? 0 : days;
}

std::string BaziBase::getTaiYuan() const
{
    if (!mPaiPanDone_) return "";
    return bazi::taiYuanStr(arraySiZhu_[2], arraySiZhu_[3]);
}

std::string BaziBase::getMingGong() const
{
    if (!mPaiPanDone_) return "";
    return bazi::mingGongStr(arraySiZhu_[0], arraySiZhu_[3], arraySiZhu_[7],
                             isAfterZhongQi());
}

bool BaziBase::isAfterZhongQi() const
{
    long double bd = mIsAst_ ? mAstJd_ : mBdJd_;
    long double mid = mIsAst_ ? mMidQiAstJd_ : mMidQiJd_;
    long double jie = mIsAst_ ? mHeadJieQiAstJd_ : mHeadJieQiJd_;
    // 中气未落在本月区间内(未排盘/异常)时按"未过中气"处理, 避免命宫误进一月。
    if (!(mid > jie)) return false;
    return bd >= mid;
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

    // 年号
    LunarYear ob(solarYear);
    std::string lyNh = ob.getNianHao();

    // 干支纪年, 以立春为界。
    // 必须与本盘年柱取同一个立春: mSLunarDay_ 来自 solar2Lunar2, 走的是"现代定气 +
    // 历谱整日"的立春, 而年柱走的是所选历法(平气定冬至/定夏至/定气)的立春**时刻**。
    // 两者在立春前后会差 1~2 天 —— 平气定冬至在 1900-2100 的 2/1~2/12 区间内有 16%
    // 的日子对不上 —— 报告里同时出现"农历庚子年"和"年柱己亥"会自相矛盾。
    // 排盘尚未执行时无年柱可用, 退回 solar2Lunar2 的结果。
    int yGanIdx = mPaiPanDone_ ? arraySiZhu_[0] : mSLunarDay_.lyGanIdx;
    int yZhiIdx = mPaiPanDone_ ? arraySiZhu_[1] : mSLunarDay_.lyZhiIdx;

    // 生肖同样以年柱地支为准, 与 getShenXiao() 保持一致
    std::string lySx = "『生肖』 " + std::string(ShengXiao[yZhiIdx]);

    std::string lyGan = Gan[yGanIdx];
    std::string lyZhi = Zhi[yZhiIdx];

    // 月信息 —— 统一走 sxwnl::lunarMonthName(见 lunar_month_name.h), 与年历/月历/日详情一致
    std::string lmMc = sxwnl::lunarMonthName(solarYear, mSLunarDay_.monthIdx + 1,
                                             mSLunarDay_.monthStyle, mSLunarDay_.isLeap);

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

    if (XianDaiNongLifa_DingQiFa == mLifa_)
    {
        s = "现代农历定气";
        return "『定气方式』 " + s;
    }

    if (YuWuWeiZiPingLifa_DingDongZhi == mLifa_)
    {
        return "『定气方式』 依据尤氏子平历定冬至计算节气交接日期";
    }

    if (YuWuWeiZiPingLifa_DingXiaZhi == mLifa_)
    {
        return "『定气方式』 依据尤氏子平历定夏至计算节气交接日期";
    }
    return "『定气方式』 未知";
}

std::string BaziBase::getDingQiType() const
{
    // 返回不含标签/"依据"前缀的定气方式简洁标签, 供经典版"依据XXX。"文案直接拼接。
    if (mLifa_ > LifaUnknown && mLifa_ < XianDaiNongLifa_DingQiFa)
    {
        return std::string(LifaList[mLifa_ - 1]) + "数据拟合";
    }
    if (XianDaiNongLifa_DingQiFa == mLifa_)
    {
        return "现代农历定气法";
    }
    if (YuWuWeiZiPingLifa_DingDongZhi == mLifa_)
    {
        return "尤氏子平历定冬至";
    }
    if (YuWuWeiZiPingLifa_DingXiaZhi == mLifa_)
    {
        return "尤氏子平历定夏至";
    }
    return "未知定气方式";
}

std::string BaziBase::getJieQiTerms() const
{
    return getJieQiIterms();
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

std::vector<std::string> BaziBase::getFleetingYearList(bool vertical) const
{
    int jyYear = mJyYear_;
    auto p = (jyYear + 6) % 10;
    if (p < 0) p += 10;
    auto q = (jyYear + 8) % 12;
    if (q < 0) q += 12;
    std::vector<std::string> output;
    std::string opt;
    
    if (vertical)
    {
        // 竖排：外层循环10次（每个大运10年），内层循环8次（8个大运）
        for (int i = 0; i < 10; i++)
        {
            int m = q;
            for (int j = 0; j < 8; j++)
            {
                opt = std::string(Gan[p]) + Zhi[m];
                output.push_back(opt);
                m = (m + 10) % 12;
            }
            p = (p + 1) % 10;
            q = (q + 1) % 12;
        }
    }
    else
    {
        // 横排：外层循环8次（8个大运），内层循环10次（每个大运10年）
        for (int i = 0; i < 8; i++)
        {
            // 内层循环10次（每个大运10年）
            for (int j = 0; j < 10; j++)
            {
                opt = std::string(Gan[p]) + Zhi[q];
                output.push_back(opt);
                
                // 天干和地支都递增1
                p = (p + 1) % 10;
                q = (q + 1) % 12;
            }
        }
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
    
    // 打印大运和流年（竖排）
    s += "=== 大运流年 ===\n\n";
    
    // 获取当前年份
    Time nowTime = JD::getNowTime();
    int currentYear = nowTime.getYear();
    
    // 获取大运列表和年份范围
    auto daYunList = getDaYunList();
    auto startYearList = getStartYearList();
    auto endYearList = getEndYearList();
    auto fleetingYearList = getFleetingYearList(true); // 使用竖排方式
    
    // 定义列宽（视觉宽度）
    const int colVisWidth = 12;
    const int labelVisWidth = 8;
    
    // 辅助lambda：生成指定数量的空格
    auto getPad = [](int count) { return std::string(count > 0 ? count : 0, ' '); };
    
    // 打印大运标题行
    s += "大运" + getPad(labelVisWidth - 4); // "大运"视觉宽度为4
    for (size_t i = 0; i < daYunList.size(); i++)
    {
        s += daYunList[i] + getPad(colVisWidth - 4); // 干支视觉宽度为4
    }
    s += "\n";
    
    // 打印年份范围行
    s += "起止年" + getPad(labelVisWidth - 6); // "起止年"视觉宽度为6
    for (size_t i = 0; i < startYearList.size(); i++)
    {
        std::ostringstream yearRange;
        yearRange << startYearList[i] << "-" << endYearList[i];
        s += yearRange.str() + getPad(colVisWidth - 9); // "YYYY-YYYY"视觉宽度为9
    }
    s += "\n";
    s += std::string(labelVisWidth + 8 * colVisWidth, '=') + "\n";
    
    // 竖排打印流年（每行代表大运内的第N年，横向展示8个大运的该年份流年）
    for (int yearOffset = 0; yearOffset < 10; yearOffset++)
    {
        // 先输出流年干支行
        s += getPad(labelVisWidth);
        
        for (int dayunIdx = 0; dayunIdx < 8; dayunIdx++)
        {
            // 竖排索引：yearOffset * 8 + dayunIdx
            int fleetingYearIndex = yearOffset * 8 + dayunIdx;
            if (fleetingYearIndex < static_cast<int>(fleetingYearList.size()))
            {
                std::string fleetingYearGanZhi = fleetingYearList[fleetingYearIndex];
                int actualYear = startYearList[dayunIdx] + yearOffset;
                
                // 判断是否为当前流年
                bool isCurrentFleetingYear = (actualYear == currentYear);
                
                if (isCurrentFleetingYear)
                {
                    // 当前流年用红色标注
                    s += "\033[1;31m" + fleetingYearGanZhi + "\033[0m" + getPad(colVisWidth - 4);
                }
                else
                {
                    s += fleetingYearGanZhi + getPad(colVisWidth - 4);
                }
            }
            else
            {
                s += getPad(colVisWidth);
            }
        }
        s += "\n";
        
        // 输出年份行
        s += getPad(labelVisWidth);
        
        for (int dayunIdx = 0; dayunIdx < 8; dayunIdx++)
        {
            int actualYear = startYearList[dayunIdx] + yearOffset;
            bool isCurrentFleetingYear = (actualYear == currentYear);
            
            std::ostringstream tempYear;
            tempYear << "(" << actualYear << ")";
            std::string yearStr = tempYear.str();
            
            if (isCurrentFleetingYear)
            {
                // 当前流年用红色标注
                s += "\033[1;31m" + yearStr + "\033[0m" + getPad(colVisWidth - 6); // "(YYYY)"视觉宽度为6
            }
            else
            {
                s += yearStr + getPad(colVisWidth - 6);
            }
        }
        s += "\n";
    }
    
    s += "\n";
    
    return s;
}
