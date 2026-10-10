#pragma once
#include <string>
#include <vector>
#include "geo.h"
#include "SSQ.h"
#include "JD.h"
#include "const.h"
#include "sxtwl.h"


typedef enum
{
    LifaUnknown,
	TaiChuLifa,
    SiFenLifa,
    DaMingLifa,
    WuYinYuanLifa,
    LinDeLifa,
    ZhengYuanLifa,
    YingTianLifa,
    ChongTianLifa,
    ChunYouLifa,
    ShouShiLifa,
    XianDaiNongLifa_DingQiFa,
    YuWuWeiZiPingLifa_DingDongZhi,
    YuWuWeiZiPingLifa_DingXiaZhi,
    
}LiFaType;

typedef enum
{
    CalendarUnknown,
    CalendarSolar,
    CalendarLunar,
} CalendarType;

typedef enum {
    UnknownPaiPanMode,
    PingQiPaiPanMode,
    DingQiPaiPanMode,
}PaiPanMode;

typedef struct
{
    std::string name;
    bool gender; // 男：false, 女：true
    bool isAst;  // 是否真太阳时
    JINGWEI jw;
    LiFaType lifa;
    CalendarType calendar;
    Time birthDayTime;
    bool isRun;  // 是否农历闰年
    bool isSpec; // 是否特殊的下一月份
} SBaziInputPara;

typedef struct
{
    std::string shishen;
    std::string tiangan;
    std::string dizhi;
    std::string canggan1;
    std::string canggan2;
    std::string canggan3;
} SBaziSiZhu;

typedef struct
{
    std::string testTime;                     // 起盘时间
    std::string userName;                     // 姓名
    std::string gender;                       // 性别
    std::string shengxiao;                    // 生肖
    std::string solarBirth;                   // 公历生日或儒略历生日
    std::string lunarBirth;                   // 农历生日
    std::string dateOfBirth;                  // 出生年代
    std::string dingQiType;                   // 定气方式
    std::string jieQi;                        // 节气
    SBaziSiZhu szInfo;                        // 八字盘面信息
    std::string qiYun;                        // 起运
    std::string jiaoYun;                      // 交运
    std::vector<std::string> vecDaYun;        // 大运
    std::vector<std::string> vecFleetingYear; // 流年
    int age;                                  // 年龄
    std::vector<int> vecStartYear;            // 起于
    std::vector<int> vecEndYear;              // 止于
} SPaiPanInfo;

class BaziBase
{
public:
    explicit BaziBase(const SBaziInputPara &input);
    ~BaziBase();
    void calcBaziPaiPan();
    
    std::string getLifa() const;
    std::string getDingQiType() const; // 定气方式(历法)简洁标签, 供排盘经典版"依据..."文案使用
    std::string getJieQiTerms() const; // 出生月首/末节交接(不含经纬度/真太阳时块)
    std::string getAst() const;
    std::string getUserName() const;
    std::string getUserGender() const;
    std::string getSolarBirth() const;
    std::tuple<std::string, std::string, std::string> getLunarInfo() const;
    
    std::string getShenXiao() const;
    std::string getAge() const;
    std::string getJieLing() const;
    std::vector<std::vector<std::string>> getSiZhu() const;
    std::string getQiYun() const;
    std::string getJiaoYun() const;
    std::vector<std::string> getDaYunList() const;
    std::vector<int> getStartYearList() const;
    std::vector<std::string> getFleetingYearList(bool vertical = false) const;
    std::vector<int> getEndYearList() const;
    std::string printBazi() const;

    // 胎元: 月柱干进一位、支进三位, 例 "丁巳"
    std::string getTaiYuan() const;
    // 命宫: 月时相加取支、年干五虎遁配干, 论月以中气为界, 例 "丙辰"
    std::string getMingGong() const;

    // ── 供派生命理分析(bazi_analysis)使用的原始数据 ──
    // 四柱索引: [0]年干[1]年支[2]月干[3]月支[4]日干[5]日支[6]时干[7]时支
    std::array<int, 8> getSiZhuIndex() const { return arraySiZhu_; }
    int  getStartYear() const { return mJyYear_; }     // 起运/首个流年的公历年
    int  getSolarYearOfBirth() const { return mSolarYear_; } // 出生公历年
    bool getGenderIsFemale() const { return mGender_; } // true=女
    int  getDaysAfterJie() const;                       // 出生距本月节令的天数

    // 出生月首/末节的儒略日(JD): 平太阳时口径 (mBdJd 时基)
    long double getHeadJieQiJd() const { return mHeadJieQiJd_; }
    long double getTailJieQiJd() const { return mTailJieQiJd_; }
    // 出生月首/末节的儒略日(JD): 真太阳时口径 (mAstJd 时基, isAst=false 时同上)
    long double getHeadJieQiAstJd() const { return mHeadJieQiAstJd_; }
    long double getTailJieQiAstJd() const { return mTailJieQiAstJd_; }
    // 出生时刻 JD (平太阳时/真太阳时)
    long double getBirthJd() const { return mBdJd_; }
    long double getBirthAstJd() const { return mAstJd_; }

    // 出生月"中气"(月首节与月末节之间那个气)的儒略日。
    // 八字以"节"划月, 但命宫论月以"中气"为界, 故单独给出。
    long double getMidQiJd() const { return mMidQiJd_; }
    long double getMidQiAstJd() const { return mMidQiAstJd_; }
    // 出生时刻是否已过本月中气 (命宫论月用: 过中气者作次月计)
    bool isAfterZhongQi() const;

private:
    // 平气排盘所需的"气年框架"。
    // 一个气年 = 相邻两个至点之间, 被 24 等分; 其中 12 个"节"位于奇数分点。
    // 由于相邻至点间隔逐年有十几分钟的起伏, 每个气年必须用自己的间隔 k, 跨气年
    // 取节时不能借用邻年的 k, 否则节的时刻会偏出几十秒。
    struct PingQiFrame
    {
        long double pzq{};   // 本气年起点至点(定冬至=冬至, 定夏至=夏至)
        long double kPrev{}; // 上一气年的节间隔 = (pzq - ppzq)/12
        long double k{};     // 本气年的节间隔   = (zq  - pzq )/12
        long double b{};     // 本气年首节 = pzq + k/2  (小寒 / 小暑)
        long double k1{};    // 下一气年的节间隔 = (nzq - zq)/12
        long double b1{};    // 下一气年首节 = zq + k1/2
    };

    void calcPingQiPaiPan();
    void calcDingQiPaiPan();
    PingQiFrame bkCalc(int year, int month, LiFaType lifa);
    long double calcLichun(int year) const;
    std::string getJieQiIterms() const;
    std::tuple<int, double, double> calcJieQiTermsByPingQi();
    long double getJieByDingQi(int year,int month);
    void calcJieQiTermsByDingQi();
    void calcJiaoYunDate(bool flag);
    long double mBdJd_;  // 平太阳时对应儒略世纪数
    long double mAstJd_; // 真太阳对应儒略世经数
    LiFaType mLifa_;
    PaiPanMode mBzCalcMode_;
    bool mGender_;
    bool mIsAst_;
    std::string mUserName_;
    JINGWEI mSUserJw_;
    SLunarDay mSLunarDay_;
    int mSolarYear_;
    int mSolarMonth_;
    int mSolarDay_;
    
    long double mHeadJieQiJd_{};
    long double mTailJieQiJd_{};
    long double mHeadJieQiAstJd_{};
    long double mTailJieQiAstJd_{};
    long double mMidQiJd_{};    // 月首节与月末节之间的中气
    long double mMidQiAstJd_{};
    long double mJyJd_{}; // 起运对应儒略世经数
    int mJyYear_;
    std::array<int, 3> arrayQiYunDelta_{};
    std::array<int, 8> arraySiZhu_{};
    // arraySiZhu_ 零初始化后恰好等于"甲子", 与真实结果无法区分, 故单列一个完成标志
    bool mPaiPanDone_{false};
    std::vector<int> vecShiShen_;
};

