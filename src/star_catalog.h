#pragma once

#include <array>
#include <string>
#include <vector>
#include "eph.h"

// ═══════════════════════════════════════════════════════════════
//  star_catalog —— 88 星座表 + 内置恒星库 + 检索/解析/hxCalc
//
//  对应上游 sxwnl/ephB.js 中的:
//    xz88        : 88 星座(中英文名 / 面积 / 中心位置 / 象限 / 族)
//    HXK         : 恒星库 (#分隔的多行记录, * 表示主星)
//    schHXK(key) : 按关键字检索星库 + 星座中心
//    getHXK(s)   : 解析星库字符串 → 数值数组
//    hxCalc(...) : 多颗恒星视位置/平位置/站心位置计算
// ═══════════════════════════════════════════════════════════════

namespace star_catalog
{

// ── 88 星座 ────────────────────────────────────────────────────
struct Constellation
{
    std::string nameAbbr;   // 例如 "仙女座And"
    long double areaSq;     // 面积(平方度)
    std::string raStr;      // "0 48.46" 等
    std::string decStr;     // "37 25.91" 等
    std::string quadFamily; // "NQ1 英仙"
    std::string nameEn;     // "Andromeda"
};

// 获取全部 88 个星座
const std::vector<Constellation>& list88();

// 按缩写检索 (大小写敏感, 例如 "And"); 找不到返回 nullptr
const Constellation* findByAbbr(const std::string &abbr);

// ── 恒星条目 ──────────────────────────────────────────────────
struct Star
{
    long double ra;       // 赤经(弧度, J2000)
    long double dec;      // 赤纬(弧度, J2000)
    long double pmRa;     // 赤经自行(弧度/世纪)
    long double pmDec;    // 赤纬自行(弧度/世纪)
    long double parallax; // 视差(弧度)
    long double mag;      // 星等
    std::string name;     // 星名
    std::string info;     // 星座/编号等信息
};

// 注册一个原始恒星库字符串 (上游 HXK 格式: # 分隔星, * 主星标志)
//   key : 库标识(如 "库0"), 必须唯一; 重复注册将覆盖
//   raw : 上游格式串
void registerLibrary(const std::string &key, const std::string &raw);

// 获取已注册库列表
std::vector<std::string> libraryKeys();

// 取出某库的全部恒星 (允许 includeAll=true 返回非 * 项)
std::vector<Star> getLibrary(const std::string &key, bool includeAll = false);

// 按关键字检索全部库. key 可匹配星名或附加信息 (如 "α" / "Lyr" / "织女")
std::vector<Star> search(const std::string &key);

// 把"度分秒"/"时分秒"字符串解析为弧度
//   isHour=true : "12 34 56.7" 视为时分秒(乘 15° 转为度)
//   isHour=false: "30 15 20" 视为度分秒
long double str2rad(const std::string &s, bool isHour);

// ── hxCalc 多颗恒星位置计算 ────────────────────────────────────
//
//   stars       : 待计算的恒星集合
//   t           : 儒略世纪TD (自 J2000 起算)
//   nutationQ   : 章动周期截断(天); 0 表示不限
//   mode        : 0=视位置(赤经/赤纬), 1=站心(方位/高度), 2=平位置
//   longitudeRad: 站心计算时的观测者地理经度(弧度); 其它模式忽略
//   latitudeRad : 站心计算时的观测者地理纬度(弧度); 其它模式忽略
struct StarResult
{
    std::string name;     // 星名 + 附加信息
    long double mag;      // 星等(原样传出)
    long double a;        // 视赤经 / 方位角 / 平赤经
    long double b;        // 视赤纬 / 高度角 / 平赤纬
};

std::vector<StarResult> hxCalc(const std::vector<Star> &stars,
                               long double t, long double nutationQ,
                               int mode,
                               long double longitudeRad = 0.0L,
                               long double latitudeRad  = 0.0L);

} // namespace star_catalog
