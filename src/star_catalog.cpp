#include "star_catalog.h"
#include "star_eph.h"
#include "const.h"
#include "eph.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <vector>

namespace star_catalog
{

// ═══════════════════════════════════════════════════════════════
//  88 星座数据 (与上游 xz88 完全一致)
// ═══════════════════════════════════════════════════════════════
namespace {
struct RawXZ { const char *nameAbbr, *area, *ra, *dec, *qfm, *nameEn; };

static const RawXZ kXZ88[] = {
    {"仙女座And", "722.278", " 0 48.46", " 37 25.91", "NQ1 英仙", "Andromeda"},
    {"唧筒座Ant", "238.901", "10 16.43", "-32 29.01", "SQ2 拉卡伊", "Antlia"},
    {"天燕座APS", "206.327", "16 08.65", "-75 18.00", "SQ3 拜耳", "Apus"},
    {"宝瓶座Aqr", "979.854", "22 17.38", "-10 47.35", "SQ4 黄道", "Aquarius"},
    {"天鹰座Aql", "652.473", "19 40.02", "  3 24.65", "NQ4 武仙", "Aquila"},
    {"天坛座Ara", "237.057", "17 22.49", "-56 35.30", "SQ3 武仙", "Ara"},
    {"白羊座Ari", "441.395", " 2 38.16", " 20 47.54", "NQ1 黄道", "Aries"},
    {"御夫座Aur", "657.438", " 6 04.42", " 42 01.68", "NQ2 英仙", "Auriga"},
    {"牧夫座Boo", "906.831", "14 42.64", " 31 12.16", "NQ3 大熊", "Bootes"},
    {"雕具座Cae", "124.865", " 4 42.27", "-37 52.90", "SQ1 拉卡伊", "Caelum"},
    {"鹿豹座Cam", "756.828", " 8 51.37", " 69 22.89", "NQ2 大熊", "Camelopardalis"},
    {"巨蟹座Cnc", "505.872", " 8 38.96", " 19 48.35", "NQ2 黄道", "Cancer"},
    {"猎犬座CVn", "465.194", "13 06.96", " 40 06.11", "NQ3 大熊", "Canes Venatici"},
    {"大犬座CMa", "380.118", " 6 49.74", "-22 08.42", "SQ2 猎户", "Canis Major"},
    {"小犬座CMi", "183.367", " 7 39.17", "  6 25.63", "NQ2 猎户", "Canis Minor"},
    {"摩羯座CAP", "413.947", "21 02.93", "-18 01.39", "SQ4 黄道", "Capricornus"},
    {"船底座Car", "494.184", " 8 41.70", "-63 13.16", "SQ2 幻之水", "Carina"},
    {"仙后座Cas", "598.407", " 1 19.16", " 62 11.04", "NQ1 英仙", "Cassiopeia"},
    {"半人马Cen", "1060.422","13 04.27", "-47 20.72", "SQ3 武仙", "Centaurus"},
    {"仙王座Cep", "587.787", "22 00.00", " 71 00.51", "NQ4 英仙", "Cepheus"},
    {"鲸鱼座Cet", "1231.411"," 1 40.10", " -7 10.76", "SQ1 英仙", "Cetus"},
    {"堰蜒座Cha", "131.592", "10 41.53", "-79 12.30", "SQ2 拜耳", "Chamaeleon"},
    {"圆规座Cir", " 93.353", "14 34.54", "-63 01.82", "SQ3 拉卡伊", "Circinus"},
    {"天鸽座Col", "270.184", " 5 51.76", "-35 05.67", "SQ1 幻之水", "Columba"},
    {"后发座Com", "386.475", "12 47.27", " 23 18.34", "NQ3 大熊", "Coma Berenices"},
    {"南冕座CrA", "127.696", "18 38.79", "-41 08.85", "SQ4 武仙", "Corona Australis"},
    {"北冕座CrB", "178.710", "15 50.59", " 32 37.49", "NQ3 大熊", "Corona Borealis"},
    {"乌鸦座Crv", "183.801", "12 26.52", "-18 26.20", "SQ3 武仙", "Corvus"},
    {"巨爵座Crt", "282.398", "11 23.75", "-15 55.74", "SQ2 武仙", "Crater"},
    {"南十字Cru", " 68.447", "12 26.99", "-60 11.19", "SQ3 武仙", "Crux"},
    {"天鹅座Cyg", "803.983", "20 35.28", " 44 32.70", "NQ4 武仙", "Cygnus"},
    {"海豚座Del", "188.549", "20 41.61", " 11 40.26", "NQ4 幻之水", "Delphinus"},
    {"剑鱼座Dor", "179.173", " 5 14.51", "-59 23.22", "SQ1 拜耳", "Dorado"},
    {"天龙座Dra", "1082.952","15 08.64", " 67 00.40", "NQ3 大熊", "Draco"},
    {"小马座Equ", " 71.641", "21 11.26", "  7 45.49", "NQ4 幻之水", "Equuleus"},
    {"波江座Eri", "1137.919"," 3 18.02", "-28 45.37", "SQ1 幻之水", "Eridanus"},
    {"天炉座For", "397.502", " 2 47.88", "-31 38.07", "SQ1 拉卡伊", "Fornax"},
    {"双子座Gem", "513.761", " 7 04.24", " 22 36.01", "NQ2 黄道", "Gemini"},
    {"天鹤座Gru", "365.513", "22 27.39", "-46 21.11", "SQ4 拜耳", "Grus"},
    {"武仙座Her", "1225.148","17 23.16", " 27 29.93", "NQ3 武仙", "Hercules"},
    {"时钟座Hor", "248.885", " 3 16.56", "-53 20.18", "SQ1 拉卡伊", "Horologium"},
    {"长蛇座Hya", "1302.844","11 36.73", "-14 31.91", "SQ2 武仙", "Hydra"},
    {"水蛇座Hyi", "243.035", " 2 20.65", "-69 57.39", "SQ1 拜耳", "Hydrus"},
    {"印第安Ind", "294.006", "21 58.33", "-59 42.40", "SQ4 拜耳", "Indus"},
    {"蝎虎座Lac", "200.688", "22 27.68", " 46 02.51", "NQ4 英仙", "Lacerta"},
    {"狮子座Leo", "946.964", "10 40.03", " 13 08.32", "NQ2 黄道", "Leo"},
    {"小狮座LMi", "231.956", "10 14.72", " 32 08.08", "NQ2 大熊", "Leo Minor"},
    {"天兔座Lep", "290.291", " 5 33.95", "-19 02.78", "SQ1 猎户", "Lepus"},
    {"天秤座Lib", "538.052", "15 11.96", "-15 14.08", "SQ3 黄道", "Libra"},
    {"豺狼座Lup", "333.683", "15 13.21", "-42 42.53", "SQ3 武仙", "Lupus"},
    {"天猫座Lyn", "545.386", " 7 59.53", " 47 28.00", "NQ2 大熊", "Lynx"},
    {"天琴座Lyr", "286.476", "18 51.17", " 36 41.36", "NQ4 武仙", "Lyra"},
    {"山案座Men", "153.484", " 5 24.90", "-77 30.24", "SQ1 拉卡伊", "Mensa"},
    {"显微镜Mic", "209.513", "20 57.88", "-36 16.49", "SQ4 拉卡伊", "Microscopium"},
    {"麒麟座Mon", "481.569", " 7 03.63", "  0 16.93", "NQ2 猎户", "Monoceros"},
    {"苍蝇座Mus", "138.355", "12 35.28", "-70 09.66", "SQ3 拜耳", "Musca"},
    {"矩尺座Nor", "165.290", "15 54.18", "-51 21.09", "SQ3 拉卡伊", "Norma"},
    {"南极座Oct", "291.045", "23 00.00", "-82 09.12", "SQ4 拉卡伊", "Octans"},
    {"蛇夫座Oph", "948.340", "17 23.69", " -7 54.74", "SQ3 武仙", "Ophiuchus"},
    {"猎户座Ori", "594.120", " 5 34.59", "  5 56.94", "NQ1 猎户", "Orion"},
    {"孔雀座Pav", "377.666", "19 36.71", "-65 46.89", "SQ4 拜耳", "Pavo"},
    {"飞马座Peg", "1120.794","22 41.84", " 19 27.98", "NQ4 英仙", "Pegasus"},
    {"英仙座Per", "614.997", " 3 10.50", " 45 00.79", "NQ1 英仙", "Perseus"},
    {"凤凰座Phe", "469.319", " 0 55.91", "-48 34.84", "SQ1 拜耳", "Phoenix"},
    {"绘架座Pic", "246.739", " 5 42.46", "-53 28.45", "SQ1 拉卡伊", "Pictor"},
    {"双鱼座Psc", "889.417", " 0 28.97", " 13 41.23", "NQ1 黄道", "Pisces"},
    {"南鱼座PsA", "245.375", "22 17.07", "-30 38.53", "SQ4 幻之水", "Piscis Austrinus"},
    {"船尾座Pup", "673.434", " 7 15.48", "-31 10.64", "SQ2 幻之水", "Puppis"},
    {"罗盘座Pyx", "220.833", " 8 57.16", "-27 21.10", "SQ2 幻之水", "Pyxis"},
    {"网罟座Ret", "113.936", " 3 55.27", "-59 59.85", "SQ1 拉卡伊", "Reticulum"},
    {"天箭座Sge", " 79.932", "19 39.05", " 18 51.68", "NQ4 武仙", "Sagitta"},
    {"人马座Sgr", "867.432", "19 05.94", "-28 28.61", "SQ4 黄道", "Sagittarius"},
    {"天蝎座Sco", "496.783", "16 53.24", "-27 01.89", "SQ3 黄道", "Scorpius"},
    {"玉夫座Scl", "474.764", " 0 26.28", "-32 05.30", "SQ1 拉卡伊", "Sculptor"},
    {"盾牌座Sct", "109.114", "18 40.39", " -9 53.32", "SQ4 武仙", "Scutum"},
    {"巨蛇座Ser", "636.928", "16 57.04", "  6 07.32", "NQ3 武仙", "Serpens"},
    {"六分仪Sex", "313.515", "10 16.29", " -2 36.88", "SQ2 武仙", "Sextans"},
    {"金牛座Tau", "797.249", " 4 42.13", " 14 52.63", "NQ1 黄道", "Taurus"},
    {"望远镜Tel", "251.512", "19 19.54", "-51 02.21", "SQ4 拉卡伊", "Telescopium"},
    {"三角座Tri", "131.847", " 2 11.07", " 31 28.56", "NQ1 英仙", "Triangulum"},
    {"南三角TrA", "109.978", "16 04.95", "-65 23.28", "SQ3 武仙", "Triangulum Australe"},
    {"杜鹃座Tuc", "294.557", "23 46.64", "-65 49.80", "SQ4 拜耳", "Tucana"},
    {"大熊座UMa", "1279.660","11 18.76", " 50 43.27", "NQ2 大熊", "Ursa Major"},
    {"小熊座UMi", "255.864", "15 00.00", " 77 41.99", "NQ3 大熊", "Ursa Minor"},
    {"船帆座Vel", "499.649", " 9 34.64", "-47 10.03", "SQ2 幻之水", "Vela"},
    {"室女座Vir", "1294.428","13 24.39", " -4 09.51", "SQ3 黄道", "Virgo"},
    {"飞鱼座Vol", "141.354", " 7 47.73", "-69 48.07", "SQ2 拜耳", "Volans"},
    {"狐狸座Vul", "268.165", "20 13.88", " 24 26.56", "NQ4 武仙", "Vulpecula"},
};
}  // namespace

const std::vector<Constellation>& list88()
{
    static std::vector<Constellation> v;
    static std::once_flag flag;
    std::call_once(flag, []() {
        for (const auto &r : kXZ88)
        {
            Constellation c;
            c.nameAbbr   = r.nameAbbr;
            c.areaSq     = std::strtold(r.area, nullptr);
            c.raStr      = r.ra;
            c.decStr     = r.dec;
            c.quadFamily = r.qfm;
            c.nameEn     = r.nameEn;
            v.push_back(c);
        }
    });
    return v;
}

const Constellation* findByAbbr(const std::string &abbr)
{
    const auto &v = list88();
    for (const auto &c : v)
    {
        // c.nameAbbr 最后 3 个 ASCII 字符就是缩写
        if (c.nameAbbr.size() >= 3 &&
            c.nameAbbr.substr(c.nameAbbr.size() - 3) == abbr)
        {
            return &c;
        }
    }
    return nullptr;
}

// ═══════════════════════════════════════════════════════════════
//  辅助: str2rad
//   "  0 48.46"     → 0h48.46m   (isHour=1) → 弧度
//   " 37 25.91"     → 37°25.91′  (isHour=0) → 弧度
//   "0 48 12.345"   → 0h48m12.345s 等
// ═══════════════════════════════════════════════════════════════
long double str2rad(const std::string &s, bool isHour)
{
    // 去首尾空白
    std::string t = s;
    size_t a = 0, b = t.size();
    while (a < b && std::isspace((unsigned char)t[a])) ++a;
    while (b > a && std::isspace((unsigned char)t[b - 1])) --b;
    t = t.substr(a, b - a);
    if (t.empty()) return 0.0L;

    // 处理符号
    long double sign = 1.0L;
    if (t[0] == '-') { sign = -1.0L; t.erase(0, 1); }
    else if (t[0] == '+') { t.erase(0, 1); }

    // 拆分按空格 / ':'
    std::vector<long double> parts;
    std::string cur;
    for (char c : t)
    {
        if (c == ' ' || c == ':' || c == '\t')
        {
            if (!cur.empty())
            {
                parts.push_back(std::strtold(cur.c_str(), nullptr));
                cur.clear();
            }
        }
        else cur.push_back(c);
    }
    if (!cur.empty()) parts.push_back(std::strtold(cur.c_str(), nullptr));

    long double v = 0.0L;
    if (parts.size() >= 1) v += parts[0];
    if (parts.size() >= 2) v += parts[1] / 60.0L;
    if (parts.size() >= 3) v += parts[2] / 3600.0L;
    v *= sign;
    // 时单位 → 度
    if (isHour) v *= 15.0L;
    // 度 → 弧度
    v = v / 180.0L * static_cast<long double>(PI);
    return v;
}

// ═══════════════════════════════════════════════════════════════
//  恒星库 (HXK)
// ═══════════════════════════════════════════════════════════════
namespace {
// 内置示例库 (来自上游 ephB.js 中的演示数据)
static const char *kBuiltinLib0 =
"库0#"
"* 0 01 57.620,- 6 00 50.68, 0.0031, -0.041, 0.008, 4.37 ,星1630 ,Psc 30 M3#"
"* 0 03 44.391,-17 20 09.58, 0.0020, -0.007, 0.014, 4.55 ,星905 ,Cet 2 B9#"
"* 0 05 20.142,- 5 42 27.45,-0.0009, 0.089, 0.025, 4.61 ,星1002 ,Psc 33 K1#"
"* 0 08 23.260, 29 05 25.54, 0.0104, -0.163, 0.034, 2.07 ,星1 ,And α B9#"
"* 0 09 10.686, 59 08 59.19, 0.0681, -0.180, 0.060, 2.28 ,星2 ,Cas β F2#"
"* 0 10 19.247, 46 04 20.17, 0.0005, 0.001, 0.003, 5.01 ,星4 ,And 22 F2#"
"* 0 11 34.421,-27 47 59.06, 0.0003, 0.016, 0.006, 5.41 ,星5 ,Scl κ2 K2#"
"* 0 11 44.010,-35 07 59.24, 0.0138, 0.115, 0.046, 5.24 ,星6 ,Scl θ F3#"
"* 0 13 14.154, 15 11 00.93, 0.0003, -0.008, 0.010, 2.83 ,星7 ,Peg γ B2#"
"* 0 14 36.165, 20 12 24.12, 0.0064, -0.001, 0.010, 4.79 ,星1004 ,Peg χ M2#"
"* 0 17 05.500, 38 40 53.87,-0.0046, -0.013, 0.013, 4.61 ,N30 ,And θ A2#"
"* 0 18 19.658, 36 47 06.79,-0.0055, -0.042, 0.023, 4.51 ,星1005 ,And σ A2#"
"* 0 18 38.258, 31 31 02.01, 0.0044, -0.004, 0.006, 5.88 ,星1006 ,Pi 0h38 A0#"
"* 0 19 25.676,- 8 49 26.14,-0.0010, -0.037, 0.011, 3.56 ,星9 ,Cet ι K2#"
"* 0 20 35.863, 8 11 24.96,-0.0003, 0.010, 0.008, 5.38 ,星1008 ,Psc 41 K3#"
"* 0 21 07.270, 37 58 06.95, 0.0049, -0.039, 0.020, 5.16 ,星1009 ,And ρ F5#"
"* 0 24 47.506, 61 49 51.80, 0.0018, -0.002, 0.004, 5.38 ,GC ,Cas 12 B9#"
"* 0 25 24.210, 1 56 22.87,-0.0010, -0.013, 0.006, 5.77 ,星1010 ,Psc 44 G5#"
"* 0 25 45.092,-77 15 15.30, 0.6689, 0.323, 0.134, 2.82 ,星11 ,Hyi β G2#"
"* 0 26 17.052,-42 18 21.55, 0.0210, -0.354, 0.042, 2.40 ,星12 ,Phe α K0#";

static const char *kBuiltinLib1 =
"库1#"
"* 0 48 22.978, 5 16 50.19, 0.0507,-1.141,0.134, 5.74,星1019,G.Psc 96 K2#"
"* 0 26 17.052,-42 18 21.55, 0.0210,-0.354,0.042, 2.40,星12 ,Phe α K0#"
"* 2 36 00.049,- 7 49 53.77, -0.0022,-0.060,0.006, 5.53,星1074,Cet 80 M0#"
"* 2 35 52.472, 5 35 35.67, -0.0019,-0.024,0.009, 4.87,星1072,Cet υ G8#"
"*18 36 27.834, 9 07 20.98, -0.0001,-0.132,0.026, 5.38,星1484,Oph 9 F5#"
"*18 36 56.338, 38 47 01.29, 0.0172, 0.288,0.129, 0.03,星699, Lyr α A0#"
"*18 37 54.426,-21 23 51.81, -0.0001,-0.066,0.010, 5.93,星1485,Sgr 83 A5#"
"*18 42 16.428,- 9 03 09.14, 0.0005,-0.002,0.017, 4.70,星1486,Sct δ F2#"
"*18 43 31.254,- 8 16 30.76, 0.0013, 0.008,0.006, 4.88,星702, Sct ε G8#";

struct LibStore
{
    std::string raw;
    std::vector<Star> all;        // 全部
    std::vector<Star> starOnly;   // 仅 *
};

static std::map<std::string, LibStore>& storage()
{
    static std::map<std::string, LibStore> s;
    return s;
}
} // namespace

// 解析"星名" "Phe α K0" 等字段, 并组装为 Star 列表
//   raw  : 原始字符串
//   onlyStar=true 时只取以 '*' 起始的记录
static std::vector<Star> parseLibrary(const std::string &raw, bool onlyStar)
{
    std::vector<Star> out;
    if (raw.empty()) return out;

    // 1. 去 \r, \n → '#'
    std::string s;
    s.reserve(raw.size());
    for (char c : raw)
    {
        if (c == '\r') continue;
        if (c == '\n') s.push_back('#');
        else s.push_back(c);
    }
    // 2. 去除"库?#" 第 1 行
    auto firstHash = s.find('#');
    if (firstHash == std::string::npos) return out;
    s = s.substr(firstHash + 1);
    // 3. # → ','
    for (auto &c : s) if (c == '#') c = ',';
    // 4. 拆分
    std::vector<std::string> fields;
    {
        std::string cur;
        for (char c : s)
        {
            if (c == ',')
            {
                // 去空白
                size_t a = 0, b = cur.size();
                while (a < b && std::isspace((unsigned char)cur[a])) ++a;
                while (b > a && std::isspace((unsigned char)cur[b - 1])) --b;
                fields.push_back(cur.substr(a, b - a));
                cur.clear();
            }
            else cur.push_back(c);
        }
        if (!cur.empty()) fields.push_back(cur);
    }
    // 5. 按 8 个一组组装
    const long double RAD = static_cast<long double>(rad);
    for (size_t i = 0; i + 7 < fields.size(); i += 8)
    {
        if (fields[i].empty() || fields[i].size() < 5) continue;
        bool starred = (fields[i][0] == '*');
        if (!starred && onlyStar) continue;
        std::string raStr = starred ? fields[i].substr(1) : fields[i];

        Star st;
        st.ra       = str2rad(raStr,        true);
        st.dec      = str2rad(fields[i + 1], false);
        st.pmRa     = std::strtold(fields[i + 2].c_str(), nullptr) / RAD * 15.0L;
        st.pmDec    = std::strtold(fields[i + 3].c_str(), nullptr) / RAD;
        st.parallax = std::strtold(fields[i + 4].c_str(), nullptr) / RAD;
        st.mag      = std::strtold(fields[i + 5].c_str(), nullptr);
        st.name     = fields[i + 6];
        st.info     = fields[i + 7];
        out.push_back(st);
    }
    return out;
}

void registerLibrary(const std::string &key, const std::string &raw)
{
    auto &slot   = storage()[key];
    slot.raw     = raw;
    slot.all     = parseLibrary(raw, false);
    slot.starOnly= parseLibrary(raw, true);
}

std::vector<std::string> libraryKeys()
{
    std::vector<std::string> v;
    for (const auto &kv : storage()) v.push_back(kv.first);
    return v;
}

std::vector<Star> getLibrary(const std::string &key, bool includeAll)
{
    static std::once_flag flag;
    std::call_once(flag, []() {
        registerLibrary("库0", kBuiltinLib0);
        registerLibrary("库1", kBuiltinLib1);
    });
    auto it = storage().find(key);
    if (it == storage().end()) return {};
    return includeAll ? it->second.all : it->second.starOnly;
}

std::vector<Star> search(const std::string &key)
{
    // 触发内置库的注册
    (void)getLibrary("库0");
    std::vector<Star> out;
    for (auto &kv : storage())
    {
        for (const auto &st : kv.second.starOnly)
        {
            if (st.name.find(key)  != std::string::npos ||
                st.info.find(key)  != std::string::npos)
            {
                out.push_back(st);
            }
        }
    }
    // 加上 88 星座中心位置
    if (key.size() == 3)
    {
        const Constellation *c = findByAbbr(key);
        if (c)
        {
            Star st{};
            // 转换分制到 0.6 (分秒) 与 JS 一致: "0 48.46" → "0 48 27.6"
            auto cvtPart = [](const std::string &str)->std::string{
                if (str.size() < 5) return str;
                std::string head = str.substr(0, str.size() >= 6 ? 6 : str.size());
                std::string tail = str.size() > 6 ? str.substr(7) : std::string();
                long double frac = std::strtold(tail.c_str(), nullptr);
                long double sec  = frac * 0.6L;
                std::ostringstream oss;
                oss << head << " " << std::fixed;
                oss.precision(1);
                oss << (double)sec;
                return oss.str();
            };
            st.ra       = str2rad(cvtPart(c->raStr),  true);
            st.dec      = str2rad(cvtPart(c->decStr), false);
            st.pmRa     = 0;
            st.pmDec    = 0;
            st.parallax = 0;
            st.mag      = 0;
            st.name     = std::string("中心") + c->nameAbbr.substr(0, c->nameAbbr.size() - 3) + "方";
            st.info     = c->nameEn;
            out.push_back(st);
        }
    }
    return out;
}

// ═══════════════════════════════════════════════════════════════
//  hxCalc
// ═══════════════════════════════════════════════════════════════
std::vector<StarResult> hxCalc(const std::vector<Star> &stars,
                               long double t, long double nutationQ,
                               int mode,
                               long double longitudeRad,
                               long double latitudeRad)
{
    std::vector<StarResult> out;

    Vector2 d{0, 0};
    long double E = 0;
    Vector3 v{0, 0, 0};
    Vector3 p{0, 0, 0};
    Vector3 sunA{0, 0, 0};
    long double gstP = 0, gst = 0;

    if (mode == 0 || mode == 1)
    {
        d    = nutation(t, nutationQ);   // 章动
        E    = hcjj(t);                  // 黄赤交角
        v    = star_eph::evSSB(t);
        p    = star_eph::epSSB(t);
        sunA = star_eph::sun2000(t, 20);
        sunA = llrConv(sunA, 84381.406L / static_cast<long double>(rad)); // 太阳赤道坐标
        gstP = pGST2(t * 36525.0L);      // 平恒星时
        gst  = gstP + d[0] * std::cos(E);// 真恒星时
    }

    const long double P_2 = static_cast<long double>(pi_2);

    for (const auto &st : stars)
    {
        Vector3 z;
        z[0] = st.ra  + st.pmRa  * t * 100.0L;
        z[1] = st.dec + st.pmDec * t * 100.0L;
        z[2] = (st.parallax != 0.0L) ? (1.0L / st.parallax) : 1e11L;
        z[0] = rad2mrad(z[0]);

        if (mode == 0 || mode == 1)
        {
            z = star_eph::ylpz(z, sunA);    // 引力偏转
            z = star_eph::scGxc(z, p, true);  // 周年视差
            z = star_eph::scGxc(z, v, false); // 光行差
            z = CDllr_J2D(t, z, std::string("P03")); // 岁差
            z = CDnutation(z, E, d[0], d[1]);        // 章动

            if (mode == 1)
            {
                // 站心坐标
                z[0] += P_2 - gst - longitudeRad;
                z = llrConv(z, P_2 - latitudeRad);
                z[0] = rad2mrad(-P_2 - z[0]);
                if (z[1] > 0) z[1] += MQC(z[1]); // 大气折射
            }
        }
        else if (mode == 2)
        {
            z = CDllr_J2D(t, z, std::string("P03"));
        }

        StarResult r;
        r.name = st.name + " " + st.info;
        r.mag  = st.mag;
        r.a    = z[0];
        r.b    = z[1];
        out.push_back(r);
    }
    return out;
}

} // namespace star_catalog
