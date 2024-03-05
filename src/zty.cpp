//
// Created by fumingzhang on 2023/7/29.
//

#include "zty.h"
#include "eph.h"
Zty::Zty()
{

}
/*
 *
 * 真太阳和平太阳时差计算方法： [1]
https://baike.baidu.com/item/%E7%9C%9F%E5%A4%AA%E9%98%B3%E6%97%B6/692540?fr=ge_ala

 */
void Zty::calc(const long double inputJd , double L)
{
    this->T = inputJd;
    this->L = L;
    this->dt = dt_T(inputJd); // TD-UT
    this->jd = inputJd - this->dt; // UT

    auto tempJd = inputJd/36525;
    auto zd = nutation2(tempJd);
    this->dL = zd[0];
    this->dE = zd[1]; // 交角章动
    this->E = hcjj(tempJd) + this->dE; // 真黄赤交角

    auto z = e_coord(tempJd, -1, -1, -1); // 地球坐标
    z[0] = rad2mrad(z[0] + M_PI + gxc_sunLon(tempJd) + this->dL); // 补上太阳光行差及章动
    z[1] = -z[1] + gxc_sunLat(tempJd); // z数组为太阳地心黄道视坐标
    this->sHJ = z[0];
    this->sHW = z[1];
    this->sR = z[2]; // 太阳视黄经,视黄纬,日地质心距
    z = llrConv(z, this->E); // 转为赤道坐标
    this->sCJ = z[0];
    this->sCW = z[1]; // 太阳视赤经,视赤纬

    auto t = tempJd / 10, t2 = t * t, t3 = t2 * t, t4 = t3 * t, t5 = t4 * t;
    auto Lon = (1753470142 + 6283319653318 * t + 529674 * t2 + 432 * t3 - 1124 * t4 - 9 * t5) / 1000000000 + M_PI - 20.5 / rad; // 修正了光行差的太阳平黄经
    Lon = rad2mrad(Lon - (this->sCJ - this->dL * std::cos(this->E))); // (修正了光行差的平黄经)-(不含dL*cos(E)的视赤经)
    if (Lon > M_PI) Lon -= pi2; // 得到时差,单位是弧度
    this->sc = Lon / pi2; // 时差(单位:日)，经度时差
    this->pty = this->jd + L / pi2; // 平太阳时
    this->zty = this->jd + L / pi2 + this->sc; // 真太阳时

    this->sqsc = (L - 2 / 3 * M_PI) / pi2; //真平时时差
}

long double Zty::getPty() const {
    return pty;
}

long double Zty::getZty() const {
    return zty;
}

long double Zty::getSqsc() const {
    return sqsc;
}

long double Zty::getSc() const {
    return sc;
}

/*

//g++ zty.cpp eph.cpp -o zty_test -std=c++20
int main() {
    Zty zty;
    double T = 2460155.4;
    double L = 110.73;

    zty.calc(T, L);

    // ... 输出计算结果（如果需要）
    std::cout<<"平太阳时："<<zty.getPty()<<std::endl;
    std::cout<<"真太阳时："<<zty.getZty()<<std::endl;

    return 0;
}
*/