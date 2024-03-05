#pragma once
#include <string>
#include <vector>
#include <sstream>
#include <locale>
#include <algorithm>

class Util
{
public:
    static std::vector<std::string> split(const std::string &str, char delimiter);
    static std::vector<std::vector<std::string>> splitList(const std::vector<std::string> &JWv,char delimiter);
    static void extractParts(const std::string &s, std::string &alphanumeric, std::string &chinese);
    static std::string extractAlphanumeric(const std::string &s);
    static bool isChinese(char32_t c);
};

