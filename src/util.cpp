#include "util.h"

std::vector<std::string> Util::split(const std::string &str, char delimiter)
{
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}

std::vector<std::vector<std::string>> Util::splitList(const std::vector<std::string> &JWv,char delimiter) {
    std::vector<std::vector<std::string>> result;
    for (const auto &str: JWv) {
        std::vector<std::string> temp;
        temp = split(str,delimiter);
        result.push_back(temp);
    }
    return result;
}

std::string Util::extractAlphanumeric(const std::string &s)
{
    std::string result;
    for (char i : s) {
        if (std::isalnum(static_cast<unsigned char>(i))) {
            result += i;
        } else {
            break;
        }
    }
    return result;
}

void Util::extractParts(const std::string &s, std::string &alphanumeric, std::string &chinese) {
    alphanumeric.clear();
    chinese.clear();

    for (char c : s) {
        if (std::isalnum(c)) {
            alphanumeric += c;
        } else if (isChinese(c)) {
            chinese += c;
        }
    }
}

bool Util::isChinese(char32_t c)
{
    return (c >= 0x4E00 && c <= 0x9FFF) || (c >= 0x3400 && c <= 0x4DBF) || (c >= 0x20000 && c <= 0x2A6DF) || (c >= 0x2A700 && c <= 0x2B73F) || (c >= 0x2B740 && c <= 0x2B81F) || (c >= 0x2B820 && c <= 0x2CEAF) || (c >= 0xF900 && c <= 0xFAFF) || (c >= 0x2F800 && c <= 0x2FA1F);
}

