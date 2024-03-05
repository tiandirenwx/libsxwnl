#pragma once
#include <vector>
#include <string>
#include "SSQ.h"
#include "const.h"
#include "JD.h"
#include "eph.h"
#include "sx_lang_zh.h"

class LunarYear
{
private:
    int year_;
    std::vector<std::string> gJNB_;
    void init();
public:
    LunarYear(int year);
    ~LunarYear();
    std::string getNianHao();
    std::string getNianLiStr();
};



