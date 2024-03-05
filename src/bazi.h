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
    std::vector<std::string> getFleetingYearList() const;
    std::vector<int> getEndYearList() const;
    std::string printBazi() const;

private:
    void calcPingQiPaiPan();
    void calcDingQiPaiPan();
    std::array<long double, 4> bkCalc(int year, int month, LiFaType lifa);
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
    long double mJyJd_{}; // 起运对应儒略世经数
    int mJyYear_;
    std::array<int, 3> arrayQiYunDelta_{};
    std::array<int, 8> arraySiZhu_{};
    std::vector<int> vecShiShen_;
};

