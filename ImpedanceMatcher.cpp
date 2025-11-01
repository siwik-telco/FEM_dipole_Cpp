#include "ImpedanceMatcher.h"
#include <cmath>

ImpedanceMatcher::ImpedanceMatcher(double R_load) : R_load(R_load) {}

double ImpedanceMatcher::targetImpedance() const {
    return R_load;
}

double ImpedanceMatcher::RL(const std::complex<double> &Z_ant) const {
    if (Z_ant.real() < 0) {
        return 0.0;
    }

    std::complex<double> Gamma = (Z_ant - R_load) / (Z_ant + R_load);
    if (std::abs(Gamma) >= 1.0) return 0.0;

    return -20.0 * std::log10(std::abs(Gamma));
}
