#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <map>


struct Area
{
	std::string n;
};

struct City
{
	std::string n;
	std::vector<Area> a;
};

struct Province
{
	std::string p;
	std::vector<City> c;
};

class WorldArea
{
public:
	static WorldArea &getInstance()
	{
		static WorldArea instance;
		return instance;
	}

	std::vector<Province> getProvinces() const;
	std::map<std::string, std::map<std::string, std::vector<Area>>> getProvincesMap() const;
	std::vector<std::string> getProvinceLists() const;
	std::map<std::string, std::vector<Area>> findProvince(const std::string &name) const;
	std::vector<Area> findCity(const std::string &cityName) const;
	void printWorldArea(const std::vector<Province> &p);
	void printWorldArea();

private:
	WorldArea();
	~WorldArea() = default;
	// 禁止复制和赋值
	WorldArea(const WorldArea &) = delete;
	WorldArea &operator=(const WorldArea &) = delete;

	std::vector<Province> provinces_array_;
	std::vector<std::string> province_list_array_;
	std::map<std::string, std::vector<Area>> city_map_;
	std::map<std::string, std::map<std::string, std::vector<Area>>> provinces_map_;
};
