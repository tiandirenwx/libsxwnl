#pragma once
#include <vector>
#include <string>
#include <array>

typedef enum
{
	QType,   //气
	SType	 //朔
} QSType;

// 阴历月的"显示风格"。用于把 SSQ 排出的数字月序翻译成正确的显示名, 供
// 八字/年历/月历/日详情各处统一使用(见 lunar_month_name.h)。
//
// 背景: SSQ 计算出的月序是数字(建子=11/冬, 建丑=12/腊, 建寅=正...), 但历史上
// 帝王多次更改月建别名, 会出现"同一数字月名需要两种不同写法"的情形:
//   * 连续同名月(如 CE 23/239 年底出现两个十二月): 后一个用 SYmc(拾贰)区分;
//   * 689-701 武周历: 建子改称"正", 建寅改称"一"(而非"正"), 需要单独的"一"写法。
// 这些只是"显示"层面的区别, 与农历↔公历互转所用的 is_spec(重月回环标记)彻底解耦,
// 因此单独用本枚举承载, 不复用 is_spec, 避免显示误用 SYmc(如 690 年误显示"拾壹月")。
enum LunarMonthNameStyle
{
	LUNAR_MONTH_NORMAL       = 0, // 正常: 用 Ymc[]  (正/二/…/十/冬/腊)
	LUNAR_MONTH_SYMC         = 1, // 连续同名月: 用 SYmc[] (壹/贰/…/拾/拾壹/拾贰)
	LUNAR_MONTH_YI           = 2, // 689-701 建寅: 显示"一"(区别于建子改称的"正")
	// 古历年终置闰月(春秋/战国"十三月", 秦汉"后九月")。
	// 单独用风格承载其身份, 使显示不依赖 leap 下标——因为 leap_month_==0 兼作
	// "本年无闰"的哨兵值, 当置闰月恰好落在年历首月(下标 0, 如公元前 255 年)时,
	// 旧的 (leap!=0 && i==leap) 判定会漏判, 导致"十三月"被误显示成"腊月"。
	LUNAR_MONTH_ANCIENT_LEAP = 3,
};


class SSQ
{
public:
    // 说明：SSQ 的常量修正表（suoKB/qiKB/SB/QB）在进程内只初始化一次并全局共享（只读，线程安全）；
    // calcY 产生的可变计算状态则保存在各实例中。getInstance 返回 thread_local 实例，
    // 使每个线程拥有独立的计算状态，从而在 Android/iOS 等多线程环境下安全使用，
    // 同时避免因每线程重建而复制体积较大的常量表。
    static SSQ& getInstance() {
        thread_local SSQ instance;
        return instance;
    }

public:
	int calc(long double jd, QSType type);
	//较高精度气;
	long double qi_high(long double);
	//较高精度朔
	long double so_high(long double);
	long double so_low(long double W);
	long double qi_low(long double W);

	void calcY(int jd);
	//气朔表解压缩（不依赖实例状态）
	static std::string jieya(std::string s);
    std::vector<long double> getZhongQi() const;
    std::vector<int> getHS() const;
    std::vector<int> getYm() const;
	std::vector<int> getYueName() const;
    std::vector<int> getDx() const;
    std::vector<bool> getSpecificLunarMonth() const;
    //各月的显示风格(LunarMonthNameStyle), 供显示层统一翻译月名
    std::vector<int> getMonthDisplayStyle() const;
	std::vector<int> getBDNS() const;
    std::vector<int> getZQPe() const;
    int getLeap() const;

private:
    SSQ();
    ~SSQ();
    SSQ(const SSQ&) = delete;
    SSQ& operator =(const SSQ&) = delete;

private:
	// 以下为每次 calcY 计算得到的可变状态，随实例（thread_local）独立保存
	std::vector<long double> zhongqi_array_; //中气表,其中.liqiu是节气立秋的儒略日,计算三伏时用到
    std::vector<int> zhongqi_pe_array_; //补算二气
	//long double ZQ_pe1, ZQ_pe2;
	std::vector<int> sun_moon_hesuo_array_; //合朔表
	std::vector<int> month_order_array_;  //月名索引
	std::vector<int> month_name_array_; //月名
	std::vector<int> month_daxiao_array_; //各月大小
    std::vector<bool> specific_next_month_array_; //特殊月序是下一个月
    std::vector<int> month_display_style_; //各月显示风格(LunarMonthNameStyle)
	std::vector<int> bd_ns_array_; //公元前19年7闰法年首，年尾，闰月情况
	int leap_month_ = 0;

};


