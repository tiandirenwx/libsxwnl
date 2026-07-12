#include "year_utils.h"

#include <cctype>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <string>

namespace year_utils
{

// ── 内部: 去除 UTF-8 "年" / 空白 ───────────────────────────────
//   "年" (U+5E74) 的 UTF-8 编码是 0xE5 0xB9 0xB4
//   "公元" -> 0xE5 0x85 0xAC 0xE5 0x85 0x83
//   "公元前" -> 公元 + 0xE5 0x89 0x8D
//   "前" -> 0xE5 0x89 0x8D
static bool startsWith(const std::string &s, const std::string &p)
{
    return s.size() >= p.size() && s.compare(0, p.size(), p) == 0;
}

static void stripTrailingNian(std::string &s)
{
    const std::string nian = "\xE5\xB9\xB4";
    while (s.size() >= nian.size() &&
           s.compare(s.size() - nian.size(), nian.size(), nian) == 0)
    {
        s.erase(s.size() - nian.size());
    }
}

static void stripWhitespace(std::string &s)
{
    std::string out;
    out.reserve(s.size());
    for (char c : s)
    {
        if (!std::isspace(static_cast<unsigned char>(c))) out.push_back(c);
    }
    s.swap(out);
}

int year2Ayear(const std::string &raw)
{
    if (raw.empty()) return INT_MIN;

    std::string s = raw;
    stripWhitespace(s);
    stripTrailingNian(s);

    // 与上游 tools.js year2Ayear 及两端 YearUtil 对齐:
    //   * B / b / * / 公元前 / 前  => 公元前记法, 天文年 = 1 - n
    //   * 前导 '-'                 => 直接天文负年, 天文年 = -n  (如 "-220" = 天文 -220)
    //   * 公元 / 纯数字             => 直接天文年, 天文年 = n
    // 注意: 历史上 C++ 曾把 '-' 也当"公元前", 导致 "-220" 被解析成天文 -219,
    //       与 JS/HarmonyOS 不一致(生肖差一位), 此处修正。
    const std::string kGongYuan    = "\xE5\x85\xAC\xE5\x85\x83";          // 公元
    const std::string kGongYuanQian= "\xE5\x85\xAC\xE5\x85\x83\xE5\x89\x8D";// 公元前
    const std::string kQian        = "\xE5\x89\x8D";                       // 前

    bool bc = false;       // 公元前记法: 天文 = 1 - n
    bool negAstro = false; // 前导 '-': 天文 = -n

    if (startsWith(s, kGongYuanQian))
    {
        s.erase(0, kGongYuanQian.size());
        bc = true;
    }
    else if (startsWith(s, kQian))
    {
        s.erase(0, kQian.size());
        bc = true;
    }
    else if (startsWith(s, kGongYuan))
    {
        s.erase(0, kGongYuan.size());
    }
    else if (!s.empty() && (s[0] == 'B' || s[0] == 'b' || s[0] == '*'))
    {
        s.erase(0, 1);
        bc = true;
    }
    else if (!s.empty() && s[0] == '-')
    {
        s.erase(0, 1);
        negAstro = true;
    }

    if (s.empty()) return INT_MIN;

    // 必须是纯数字
    for (char c : s)
    {
        if (!std::isdigit(static_cast<unsigned char>(c))) return INT_MIN;
    }

    long long n = 0;
    try
    {
        n = std::stoll(s);
    }
    catch (...)
    {
        return INT_MIN;
    }

    if (bc)
    {
        // 公元前从 B.C.1 起, 没有公元前 0 年
        if (n < 1) return INT_MIN;
        long long a = 1 - n;
        if (a < INT_MIN) return INT_MIN;
        return static_cast<int>(a);
    }
    if (negAstro)
    {
        return static_cast<int>(-n);
    }
    if (n > INT_MAX) return INT_MAX;
    return static_cast<int>(n);
}

std::string Ayear2year(int aYear, bool fullStyle)
{
    if (aYear > 0)
    {
        if (fullStyle)
            return std::string("\xE5\x85\xAC\xE5\x85\x83") + std::to_string(aYear) + "\xE5\xB9\xB4";
        return std::to_string(aYear);
    }
    // aYear <= 0 -> 公元前 (1 - aYear)
    long long hist = 1 - static_cast<long long>(aYear);
    std::string num = std::to_string(hist);
    if (fullStyle)
        return std::string("\xE5\x85\xAC\xE5\x85\x83\xE5\x89\x8D") + num + "\xE5\xB9\xB4";
    return std::string("\xE5\x89\x8D") + num;
}

double timeStr2hour(const std::string &raw)
{
    if (raw.empty()) return std::nan("");
    std::string s = raw;
    stripWhitespace(s);

    // 拆分 ':'
    int parts[3] = {0, 0, 0};
    double secf = 0.0;
    int colon = 0;
    int idx = 0;
    std::string cur;

    auto tryParseInt = [](const std::string &x, int &out) -> bool {
        if (x.empty()) return false;
        for (char c : x)
            if (!std::isdigit(static_cast<unsigned char>(c)) && c != '-') return false;
        try { out = std::stoi(x); } catch (...) { return false; }
        return true;
    };

    // 计数冒号
    for (char c : s) if (c == ':') ++colon;

    if (colon == 0)
    {
        // 纯数字: 当作小数小时
        try
        {
            size_t pos = 0;
            double v = std::stod(s, &pos);
            if (pos != s.size()) return std::nan("");
            return v;
        }
        catch (...)
        {
            return std::nan("");
        }
    }

    for (size_t i = 0; i <= s.size(); ++i)
    {
        if (i == s.size() || s[i] == ':')
        {
            if (idx >= 3) return std::nan("");
            if (idx < 2)
            {
                if (!tryParseInt(cur, parts[idx])) return std::nan("");
            }
            else
            {
                try
                {
                    size_t pos = 0;
                    secf = std::stod(cur, &pos);
                    if (pos != cur.size()) return std::nan("");
                }
                catch (...) { return std::nan(""); }
            }
            cur.clear();
            ++idx;
        }
        else
        {
            cur.push_back(s[i]);
        }
    }

    double h = parts[0];
    double m = parts[1];
    double sec = (colon >= 2) ? secf : 0.0;

    return h + m / 60.0 + sec / 3600.0;
}

} // namespace year_utils
