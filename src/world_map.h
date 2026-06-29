#pragma once

#include <vector>
#include <string>
#include <cstddef>

// ═══════════════════════════════════════════════════════════════
//  world_map —— 世界地图数据 (海岸/岛屿/国界)
//
//  数据来源: 上游 sxwnl/eph0.js 中的 ditu0/ditu1/ditu2 三套数据
//
//  ditu0 (小图, 像素 2009 x 970): 已内置, 直接可用
//  ditu1 (大图海岸/岛屿)
//  ditu2 (国界)
//    -> 数据较大(共约 50 KB), 调用 setMapData() 由用户从文件/资源载入,
//       内置静态库默认仅含解码器和 ditu0。
//
//  输出格式: std::vector<long double>, 经度/纬度交替, 1e7 作为 moveto 分隔符。
//  精度: 与 libsxwnl 其它天文计算保持一致, 全部使用 long double。
// ═══════════════════════════════════════════════════════════════

namespace world_map
{

// 解码 dituJM 字符串
//   p   原始数据串
//   jb  经度倍率(弧度) = 2π / 图宽
//   wb  纬度倍率(弧度) = π / 图高
//   返回经纬度交替向量, 1e7 表示 moveto
std::vector<long double> dituJM(const std::string &p, long double jb, long double wb);

// ditu0: 小图(2009 x 970 像素), 已内置
const std::vector<long double>& ditu0();

// 设置/获取自定义地图数据
//   idx = 1 → ditu1 (4200 x 2100 大图海岸)
//   idx = 2 → ditu2 (国界)
//   data: 上游原始字符串
//   返回 true 设置成功
bool setMapData(int idx, const std::string &data);

// 返回已设置的解码后数据; 未设置时返回空 vector
const std::vector<long double>& getMapData(int idx);

// 释放自定义数据
void clearMapData(int idx);

} // namespace world_map
