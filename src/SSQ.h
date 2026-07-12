#pragma once
#include <vector>
#include <string>
#include <array>

typedef enum
{
	QType,   //气
	SType	 //朔
} QSType;


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
	std::vector<int> bd_ns_array_; //公元前19年7闰法年首，年尾，闰月情况
	int leap_month_ = 0;

};


