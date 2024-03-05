#pragma once
#include <cctype>
#include <cstdio>
#include <vector>
#include <map>
#include <cstring>
#include <string>

struct JINGWEI
{
	double J;	// 经度
	double W;	// 纬度
	std::string s; // 省市
	std::string x; // 区县
};

class GeoPostion
{
public:
	static GeoPostion &getInstance()
	{
		static GeoPostion instance;
		return instance;
	}
	void init();
	void JWdecode(const std::string &v,double &jin,double &wei);
	void SQdecode();
	JINGWEI getCityGeoPos(const std::string &pName,const std::string &aName) const;
	JINGWEI getCityGeoPos() const;
	JINGWEI getDefaultGeoPos() const;

private:
	GeoPostion();
	~GeoPostion();
	GeoPostion(const GeoPostion &) = delete;
	GeoPostion &operator=(const GeoPostion &) = delete;

private:
	std::map<std::string,JINGWEI> jwMap;
};
