//
// Created by Bartek Siwik on 16/10/2025.
//

#ifndef IMPEDANCEMATCHER_H
#define IMPEDANCEMATCHER_H
#include <complex>
#include "dipoleAntenna.h"


class ImpedanceMatcher {
public:
    ImpedanceMatcher(double R_load);
    double RL(const std::complex<double>& Z_ant) const;
    double targetImpedance();

private:
    double R_load;
};



#endif //IMPEDANCEMATCHER_H
