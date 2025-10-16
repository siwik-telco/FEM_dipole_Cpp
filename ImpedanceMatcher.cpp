//
// Created by Bartek Siwik on 16/10/2025.
//

#include "ImpedanceMatcher.h"

ImpedanceMatcher::ImpedanceMatcher(double R_load) : R_load(R_load) {}

double ImpedanceMatcher::targetImpedance() {
    return R_load;
}

double ImpedanceMatcher::RL(const std::complex<double> &Z_ant) const {
    std::complex<double> Gamma = (Z_ant - R_load) / (Z_ant + R_load);
    return 20 * std::log(std::abs(Gamma));
}


