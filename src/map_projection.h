#pragma once

#include <vector>
#include <cstddef>

// ═══════════════════════════════════════════════════════════════
//  map_projection —— 9 种地图投影 + Mollweide 插值 + 折线变换
//
//  对应上游 sxwnl/eph0.js 中 touY 对象的实现, 用于把(经度,纬度)
//  序列投影到平面 (x, y), 单位归一化到 [-1, 1] 区间。
//
//  投影类型:
//    0  平面正投(斜轴)
//    1  斜轴等距方位
//    2  斜轴等积方位
//    3  斜轴等角方位
//    4  摩尔威特(Mollweide)
//    5  正轴等距圆柱
//    6  正轴等角圆柱(Mercator)
//    7  多圆锥
//    8  中国灯笼投影
//
//  精度: 内部计算全部使用 long double, 与 libsxwnl 其它模块一致。
// ═══════════════════════════════════════════════════════════════

namespace map_projection
{

struct Point
{
    long double x;
    long double y;
    long double z;   // z<0 表示该点不可见(背面 / 失真过大)
};

// 矩形可见窗口: [cx-rx, cx+rx] x [cy-ry, cy+ry]
struct Window
{
    long double cx;
    long double cy;
    long double rx;
    long double ry;
};

class Projector
{
public:
    Projector();

    // 设置投影参数
    //   tylx :  0..8 投影类型
    //   J0   :  基准经度(弧度)
    //   W0   :  基准纬度(弧度)
    //   win  :  可见窗口(默认 cx=cy=0, rx=ry=1)
    void setlx(int tylx, long double J0, long double W0,
               const Window &win = {0.0L, 0.0L, 1.0L, 1.0L});

    int          type() const { return tylx_; }
    long double  J0()   const { return J0_;   }
    long double  W0()   const { return W0_;   }

    // 把单个经纬度点投影到 (x, y, z) 平面坐标
    Point toxy(long double J, long double W) const;

    // 折线变换:
    //   输入是 [J, W, J, W, ..., 1e7, J, W, ...] 格式的"曲线序列",
    //   1e7 作为 moveto(line break) 标记,
    //   输出按相同格式, 且自动剔除窗口外/背面点并避免水平回扫。
    std::vector<long double> lineArr(const std::vector<long double> &d) const;

    // 同上, 但接受连续点(无 moveto), 输入是一条完整折线
    std::vector<long double> lineSingle(const long double *jw, std::size_t pairCount) const;

private:
    int         tylx_;
    long double J0_;
    long double W0_;
    long double cosE0_;
    long double sinE0_;
    Window      win_;

    // Mollweide 插值表, 缓存
    mutable std::vector<long double> mollY_;
    long double mollCZ(long double W) const;

    // 各投影实现
    Point toxy0(long double J, long double W) const;
    Point toxy1(long double J, long double W) const;
    Point toxy2(long double J, long double W) const;
    Point toxy3(long double J, long double W) const;
    Point toxy4(long double J, long double W) const;
    Point toxy5(long double J, long double W) const;
    Point toxy6(long double J, long double W) const;
    Point toxy7(long double J, long double W) const;
    Point toxy8(long double J, long double W) const;
};

// 常量: moveto 分隔符
constexpr long double kMoveTo = 1e7L;

} // namespace map_projection
