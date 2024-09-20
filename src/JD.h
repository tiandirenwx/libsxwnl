#pragma once
#include <cstdint>
#include <ctime>
#include <string>
#include <tuple>

struct Time
{
	int Y, M, D;
	double h, m, s;
	Time() = default;
	Time(int year, int month, int day, double hour, double min, double sec)
	{
		this->Y = year;
		this->M = month;
		this->D = day;
		this->h = hour;
		this->m = min;
		this->s = sec;
	}

	int getYear() const
	{
		return Y;
	}

	void setYear(int year)
	{
		Y = year;
	}

	void setMonth(int month)
	{
		M = month;
	}

	int getMonth() const
	{
		return M;
	}

	int getDay() const
	{
		return D;
	}

	void setDay(int day)
	{
		D = day;
	}

	double getHour() const
	{
		return h;
	}

	void setHour(double hour)
	{
		h = hour;
	}

	double getMin() const
	{
		return m;
	}

	void setMour(double min)
	{
		m = min;
	}

	double getSec() const
	{
		return s;
	}

	void setSec(double sec)
	{
		s = sec;
	}
};

const double MJD_0 = 2400000.5;
const double MJD_JD2000 = 51544.5;

class JD
{
public:
	// 公历转儒略日
	static long double DD2JD(int y, uint8_t m, long double d);
	static Time JD2DD(long double jd);

	static long double toJD(const Time &time);
	static std::string timeStr(long double jd);
	static std::string timeStr(const Time &t);
	static std::string formatStr(long double jd);
	static std::string JD2str(long double jd);
	static Time getNowTime();
	static std::tuple<long double, Time,std::string> calcAST(Time &dt,long double lon);
	static long double calcAST(long double jd, long double lon);
    
	static std::tuple<double, double> gcal2jd(int year, int month, int day);
	static std::tuple<int, int, int, double> jd2gcal(double jd1, double jd2);
	static std::tuple<int, int, int, double> jd2jcal(double jd1, double jd2);
	static std::tuple<double, double> jcal2jd(int year, int month, int day);
};
