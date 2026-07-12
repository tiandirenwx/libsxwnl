#pragma once
#include <cctype>
#include <cstdio>
#include <vector>
#include <map>
#include <cstring>
#include <string>

struct JINGWEI
{
	double J;	// 经度 (度, 东+ 西-)
	double W;	// 纬度 (度, 北+ 南-)
	std::string s; // 省市
	std::string x; // 区县
	double tz = 8.0; // 时区 (小时, 东+ 西-); 国内统一 8, 国际按 SQv
};

// 一个国家/地区的时区记录 (解析自 SQv)
struct TimezoneGroup
{
	std::string continent;        // "亚洲"/"欧洲"/"北美洲"/...
	std::string country;          // "日本"/"加拿大东部时区"...
	double timezone;              // 8 / -5 / 5.5 ...
	std::vector<std::string> cities; // 该国常用城市名(展示用; 不含坐标)
};

class GeoPostion
{
public:
	// GeoPostion 是只读的经纬度查表：jwMap 仅在构造时由 init() 建立，之后全部为 const 只读访问。
	// 因此作为共享单例本身即线程安全（C++11 保证单例静态初始化线程安全，且并发只读 std::map 安全）。
	// 注意：init()/SQdecode() 会修改内部状态，已置为 private，禁止在构造完成后再调用，
	// 以保证“构造后不可变”这一线程安全前提。
	static GeoPostion &getInstance()
	{
		static GeoPostion instance;
		return instance;
	}
	void JWdecode(const std::string &v,double &jin,double &wei) const;
	JINGWEI getCityGeoPos(const std::string &pName,const std::string &aName) const;
	JINGWEI getCityGeoPos() const;
	JINGWEI getDefaultGeoPos() const;

	// ── 目录查询 (供 GUI 使用) ──────────────────────────
	// 列出所有省/直辖市/自治区 (按 JWv 原始顺序).
	std::vector<std::string> listProvinces() const;
	// 列出某省内所有城市/区县 (按 JWv 原始顺序). 缺省返回空表.
	std::vector<JINGWEI> listCitiesIn(const std::string &province) const;
	// 按关键词模糊匹配区县名 (子串). 最多返回 limit 条.
	std::vector<JINGWEI> search(const std::string &keyword, int limit = 64) const;
	// 解析后的国际时区分组 (大洲 → 国家 → 时区 + 城市名).
	const std::vector<TimezoneGroup>& timezoneGroups() const { return tzGroups; }

private:
	GeoPostion();
	~GeoPostion();
	GeoPostion(const GeoPostion &) = delete;
	GeoPostion &operator=(const GeoPostion &) = delete;

	// 仅供构造时初始化查表使用，不可在构造后调用（保证“构造后不可变”的线程安全前提）
	void init();
	void SQdecode();

	// SQv 解码内部实现 (init() 调用一次).
	void buildTimezoneGroups();

private:
	std::map<std::string,JINGWEI> jwMap;
	// 按 JWv 原序保存的省→城市列表, 用于 GUI 顺序枚举.
	std::vector<std::pair<std::string, std::vector<JINGWEI>>> orderedProvinces;
	// 解析后的国际时区表.
	std::vector<TimezoneGroup> tzGroups;
};
