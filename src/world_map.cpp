#include "world_map.h"
#include "world_map_data_full.h"
#include "const.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <unordered_map>

namespace world_map
{

// ─── dituJM 解码器 ──────────────────────────────────────────────
//
//  规则与 JS 完全一致:
//   单字符 'a'..'z' 编码 0..25 (差分量)
//   单字符 'A'..'Z' 编码 0..-25
//   '#' 表示新线段开始 (输出 1e7, 累加经纬度归零)
//   多字符: 直到遇到 ',' 表示正数, 遇到 ':' 表示负数
//
std::vector<long double> dituJM(const std::string &p, long double jb, long double wb)
{
    std::vector<long double> a;
    a.reserve(p.size() * 2);

    long double J = 0.0L, W = 0.0L;
    int k2 = 0;

    for (std::size_t i = 0; i < p.size(); )
    {
        char c = p[i];
        long b = 0;

        if (c == '#')
        {
            J = 0.0L;
            W = 0.0L;
            k2 = 0;
            a.push_back(1e7L);
            ++i;
            continue;
        }

        if (c >= 'a' && c <= 'z')
        {
            b = c - 'a';
            ++i;
        }
        else if (c >= 'A' && c <= 'Z')
        {
            b = 'A' - c;  // 与 JS "65 - charCode" 等价
            ++i;
        }
        else
        {
            // 多字符数字: 累积直到 ',' 或 ':'
            std::string num;
            num.push_back(c);
            ++i;
            bool negative = false;
            while (i < p.size())
            {
                char d = p[i];
                if (d == ',') { ++i; break; }
                if (d == ':') { ++i; negative = true; break; }
                num.push_back(d);
                ++i;
                if (num.size() > 20) break;
            }
            b = std::atol(num.c_str());
            if (negative) b = -b;
        }

        ++k2;
        if (k2 % 2 == 1)
        {
            J += b;
            a.push_back(J * jb);
        }
        else
        {
            W += b;
            a.push_back(W * wb);
        }
    }
    return a;
}

// ═══════════════════════════════════════════════════════════════
//  ditu0 (小图, CSV 像素坐标, 2009 x 970)
// ═══════════════════════════════════════════════════════════════
//
//  数据格式: 经度像素, 纬度像素 交替, "1e7" 是分隔符。
//  像素 ↔ 弧度换算:
//      lon = (px / 2009 - 0.5) * 2π
//      lat = (0.5 - py / 970)  * π
//
//  此处仅嵌入小图; 大图 ditu1/ditu2 由调用方通过 setMapData 提供。
static const char *kDitu0CSV =
"1e7,2,212,58,180,121,180,128,143,150,129,130,129,99,147,92,171,74,170,42,150,122,104,175,102,230,118,223,124,192,125,215,138,405,92,413,99,411,126,445,117,422,100,590,65,639,76,606,88,724,88,741,101,818,93,942,111,1010,111,1064,126,1045,134,1003,130,1003,146,918,159,896,189,881,189,879,170,926,146,922,142,861,159,793,160,767,176,776,181,793,180,793,213,758,253,740,255,720,270,731,299,714,297,713,282,702,273,689,273,688,263,667,274,673,283,689,283,676,296,689,325,678,353,642,373,604,373,596,387,619,422,598,435,569,417,564,440,594,475,586,479,551,441,554,415,531,400,519,371,459,401,453,433,437,441,390,359,335,353,278,328,273,336,299,361,317,354,340,375,259,420,206,339,187,334,240,427,291,431,269,470,220,511,228,566,163,663,108,672,72,575,85,524,41,453,2,454,2,290,65,286,67,308,193,320,208,292,163,291,158,273,203,260,231,267,234,253,208,229,175,229,161,260,129,282,86,232,74,238,99,265,94,275,54,246,2,272,2,212,"
"1e7,67,54,115,46,154,50,155,54,119,58,100,71,67,54,"
"1e7,300,96,321,77,385,68,388,71,315,91,327,101,317,105,300,96,"
"1e7,522,45,547,44,592,58,573,62,550,60,522,45,"
"1e7,775,75,780,82,808,83,818,74,806,72,775,75,"
"1e7,809,180,799,188,800,232,811,230,809,180,"
"1e7,802,238,790,256,788,284,748,295,738,310,744,317,754,302,791,301,797,269,816,244,802,238,"
"1e7,686,354,681,357,677,365,683,371,686,354,"
"1e7,544,460,583,487,602,515,640,527,597,528,544,460,"
"1e7,662,452,619,482,638,507,646,507,672,484,665,469,673,460,662,452,"
"1e7,680,480,671,513,684,514,703,479,680,480,"
"1e7,747,492,750,504,808,532,846,537,820,507,776,494,747,492,"
"1e7,771,548,708,560,641,604,659,675,746,653,790,690,842,689,865,666,868,628,802,541,796,576,767,564,771,548,"
"1e7,973,672,981,693,941,736,952,742,1003,693,973,672,"
"1e7,282,550,257,571,254,619,268,619,289,564,282,550,"
"1e7,1092,112,1117,120,1119,126,1091,124,1082,130,1118,134,1091,146,1135,163,1105,176,1107,178,1140,169,1172,151,1255,160,1316,204,1322,278,1425,388,1517,419,1550,449,1579,450,1560,516,1587,566,1620,582,1593,744,1614,786,1644,785,1621,773,1643,745,1634,739,1771,605,1815,536,1815,519,1629,429,1557,442,1543,405,1522,404,1526,377,1486,392,1464,371,1477,332,1546,326,1559,351,1559,320,1618,256,1672,238,1639,204,1685,203,1699,192,1657,162,1582,146,1574,200,1502,169,1483,153,1573,108,1607,115,1607,126,1582,134,1644,145,1667,119,1559,96,1558,75,1666,31,1593,27,1334,68,1329,109,1259,112,1142,101,1092,112,"
"1e7,1680,39,1617,62,1634,77,1684,75,1734,110,1721,121,1744,151,1774,156,1788,131,1890,101,1864,88,1897,88,1906,54,1947,36,1927,36,1862,49,1892,37,1831,41,1836,33,1900,28,1833,22,1750,39,1680,39,"
"1e7,1889,122,1884,136,1914,138,1933,130,1889,122,"
"1e7,1969,175,1959,182,1958,198,1979,192,1980,179,1969,175,"
"1e7,1987,162,1980,172,1995,182,1976,203,2010,201,2010,180,1996,173,2002,166,1987,162,"
"1e7,2010,208,1989,215,2001,222,2001,247,1961,249,1962,287,1999,285,2010,273,2010,208,"
"1e7,2010,297,1980,297,1918,361,1916,415,1943,452,1970,466,2010,456,2010,297";

static constexpr long double kDitu0Width  = 2009.0L;
static constexpr long double kDitu0Height =  970.0L;

// 像素 (px, py) → 弧度 (lon, lat)
//   px∈[0,2009] → lon∈[-π, π]    (左 = 180°W, 右 = 180°E)
//   py∈[0, 970] → lat∈[π/2,-π/2]  (上 = 90°N, 下 = 90°S)
static const std::vector<long double>& buildDitu0()
{
    static std::vector<long double> data;
    static std::once_flag flag;
    std::call_once(flag, []() {
        const long double P    = static_cast<long double>(PI);
        const long double TWO  = static_cast<long double>(pi2);
        const char *p = kDitu0CSV;
        std::string tok;
        auto flushToken = [&](int field) {
            if (tok.empty()) return;
            if (tok == "1e7")
            {
                data.push_back(1e7L);
                tok.clear();
                return;
            }
            long double v = std::strtold(tok.c_str(), nullptr);
            if (field == 0)
            {
                // 经度像素
                long double lon = (v / kDitu0Width - 0.5L) * TWO;
                data.push_back(lon);
            }
            else
            {
                // 纬度像素
                long double lat = (0.5L - v / kDitu0Height) * P;
                data.push_back(lat);
            }
            tok.clear();
        };

        int field = 0;     // 0 = lon, 1 = lat
        for (; *p; ++p)
        {
            if (*p == ',')
            {
                bool wasMove = (tok == "1e7");
                flushToken(field);
                if (!wasMove) field = 1 - field;
                continue;
            }
            tok.push_back(*p);
        }
        flushToken(field);
    });
    return data;
}

const std::vector<long double>& ditu0() { return buildDitu0(); }

// ═══════════════════════════════════════════════════════════════
//  自定义地图(ditu1/ditu2)
// ═══════════════════════════════════════════════════════════════
struct MapStore
{
    std::vector<long double> data;
};

static std::unordered_map<int, MapStore> &storage()
{
    static std::unordered_map<int, MapStore> s;
    return s;
}

bool setMapData(int idx, const std::string &raw)
{
    if (idx != 1 && idx != 2) return false;
    const long double P   = static_cast<long double>(PI);
    const long double TWO = static_cast<long double>(pi2);
    long double jb = 0, wb = 0;
    if (idx == 1 || idx == 2)
    {
        jb = TWO / 4200.0L;
        wb = P   / 2100.0L;
    }
    storage()[idx].data = dituJM(raw, jb, wb);
    return true;
}

const std::vector<long double>& getMapData(int idx)
{
    // 首次访问时, 自动从内建 world_map_data_full.h 注入数据.
    static std::once_flag ditu1_flag, ditu2_flag;
    if (idx == 1) {
        std::call_once(ditu1_flag, []() {
            setMapData(1, std::string(world_map_builtin::DITU1));
        });
    } else if (idx == 2) {
        std::call_once(ditu2_flag, []() {
            setMapData(2, std::string(world_map_builtin::DITU2));
        });
    }

    auto it = storage().find(idx);
    if (it == storage().end())
    {
        static const std::vector<long double> empty;
        return empty;
    }
    return it->second.data;
}

void clearMapData(int idx)
{
    storage().erase(idx);
}

} // namespace world_map
