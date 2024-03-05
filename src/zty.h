//
// Created by fumingzhang on 2023/7/29.
//

#ifndef SXTWL_ZTY_H
#define SXTWL_ZTY_H

#include <iostream>
#include <string>
#include <vector>

class Zty {
public:
    Zty();
    ~Zty() = default;
    void calc(const long double inputJd , double L);
    long double getPty() const;
    long double getZty() const;
    long double getSqsc() const;
    long double getSc() const;

private:

    long double T, L, dt, jd, dL, dE, E, sHJ, sHW, sR, sCJ, sCW, sc, pty, zty, sqsc;
};




#endif //SXTWL_ZTY_H
