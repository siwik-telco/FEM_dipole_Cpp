//
// Created by Bartek Siwik on 16/10/2025.
//

#include "dipoleAntenna.h"
const double C = 3e8;

dipoleAntenna::dipoleAntenna(double lenght, double width) : width(width), lenght(lenght), radius(width/2) {}

double dipoleAntenna::radationResistance(double f) const {
    double lambda = C/f;
    double kl = 2 * PI * lenght/lambda;
    if (std::abs(kl-PI) < 1e-6) {
        return 73.1;
    }
    double kL_over_2 = kl / 2.0;
    double ci = std::cos(kL_over_2);
    double si = std::sin(kL_over_2);

    double t1 = (std::cos(kl) - 1.0) / (si*si);
    double t2 = (kl * (si - ci - 0.5 * std::sin(kl)) / (si*si));
    return (Et0 / (2 * PI)) * (t1+t2);
}

double dipoleAntenna::reactance(double f) const {
    double w = 2 * PI * f;
    // double X = (Et0/ (2*PI)) *  (2*std::log(w*lenght/(C* PI *radius))- 2.25);
    double X = (Et0/ (2*PI)) * (2*std::log10(w*lenght/(C* PI *radius))- 2.25);

    double ohm = 2 * radius / lenght;
    if (ohm> 0.001) {
         X *= (1.0-0.2 * std::log(1.0/ohm));
    }
    return X;
}

double dipoleAntenna::inductance(double f) const {
    double w = 2 * PI * f;
    double x = reactance(f);
    if (x>0) {
        return (x/w) * 1e9;
    }
    else {
        return 0.0;
    }
}
std::complex<double> dipoleAntenna::inputZ(double f) const {
    double R = radationResistance(f);
    double X = reactance(f);
    return std::complex<double>(R,X);
    }





