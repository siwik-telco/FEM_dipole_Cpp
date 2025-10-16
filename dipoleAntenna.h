//
// Created by Bartek Siwik on 16/10/2025.
//

#ifndef DIPOLEANTENNA_H
#define DIPOLEANTENNA_H
#include <math.h>
#include <iostream>
#include <algorithm>
#include <complex>
const double Et0 = 376.73;
const double PI = M_PI;


class dipoleAntenna {
public:
    dipoleAntenna(double lenght, double width);
    double radationResistance(double f) const;
    double reactance(double f) const;
    std::complex<double> inputZ(double f) const;
    double inductance(double f) const;
private:
    double lenght;
    double width;
    double radius;
};



#endif //DIPOLEANTENNA_H
