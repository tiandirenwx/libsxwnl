#include "util.h"
#include <regex>

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

    // 注意: JWv 字符串是 UTF-8 编码的, 不能逐 char 判断 (高位字节会被符号扩展
    // 成负数, 导致 isChinese 永远返回 false). 这里按 UTF-8 首字节确定字符长度,
    // 把整段非 ASCII 字节直接当作"中文"附加到 chinese.
    size_t i = 0;
    while (i < s.size()) {
        unsigned char c = static_cast<unsigned char>(s[i]);
        if (c < 0x80) {
            if (std::isalnum(c)) alphanumeric += static_cast<char>(c);
            ++i;
            continue;
        }
        int len = 1;
        if      ((c & 0xE0u) == 0xC0u) len = 2;
        else if ((c & 0xF0u) == 0xE0u) len = 3;
        else if ((c & 0xF8u) == 0xF0u) len = 4;
        if (i + static_cast<size_t>(len) > s.size()) break; // 防越界
        chinese.append(s, i, len);
        i += static_cast<size_t>(len);
    }
}

bool Util::isChinese(char32_t c)
{
    return (c >= 0x4E00 && c <= 0x9FFF) || (c >= 0x3400 && c <= 0x4DBF) || (c >= 0x20000 && c <= 0x2A6DF) || (c >= 0x2A700 && c <= 0x2B73F) || (c >= 0x2B740 && c <= 0x2B81F) || (c >= 0x2B820 && c <= 0x2CEAF) || (c >= 0xF900 && c <= 0xFAFF) || (c >= 0x2F800 && c <= 0x2FA1F);
}

std::vector<std::string> Util::replaceAndSplit(const std::string &s, char delimiter,const std::string &pattern)
{
    std::vector<std::string> result;
    std::regex regexPattern(pattern); // 正则表达式匹配
    std::string cleaned = std::regex_replace(s, regexPattern, ""); // 移除正则表达

    std::stringstream ss(cleaned);
    std::string item;
    while (std::getline(ss, item, delimiter)) {
        if (!item.empty()) { // 忽略空字符串
            result.push_back(item);
        }
    }
    return result;
}

