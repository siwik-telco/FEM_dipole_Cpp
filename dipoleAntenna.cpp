#include "dipoleAntenna.h"
#include <cmath>

const double C = 3e8;

dipoleAntenna::dipoleAntenna(double lenght, double width)
    : width(width), lenght(lenght), radius(width/2) {}

double dipoleAntenna::radationResistance(double f) const {
    double lambda = C / f;
    double L_over_lambda = lenght / lambda;
    double kl = 2.0 * PI * L_over_lambda;

    // Use empirical formulas based on electrical length
    // For near half-wave dipole (most common case)
    if (L_over_lambda >= 0.4 && L_over_lambda <= 0.6) {
        // Near half-wave: use interpolation around 73 Ohms
        double deviation = (L_over_lambda - 0.5) * 200.0;
        return 73.1 + deviation;
    }

    // For short dipoles (L < 0.4λ)
    if (L_over_lambda < 0.4) {
        return 20.0 * PI * PI * L_over_lambda * L_over_lambda;
    }

    // For longer dipoles (0.6λ < L < 1.2λ)
    // Radiation resistance formula with current maximum
    double sin_kl_half = std::sin(kl / 2.0);

    // Prevent division by zero
    if (std::abs(sin_kl_half) < 1e-6) {
        return 73.1;
    }

    // Base radiation resistance (at current maximum)
    double Rr_max = 73.1; // Approximate for half-wave

    // Input resistance with current distribution
    double Rin = Rr_max / (sin_kl_half * sin_kl_half);

    return std::max(1.0, Rin);
}

double dipoleAntenna::reactance(double f) const {
    double lambda = C / f;
    double kl = 2.0 * PI * lenght / lambda;
    double ka = 2.0 * PI * radius / lambda;

    // Thin wire approximation for reactance
    double Omega = std::log(lenght / radius) - 1.0;

    // Reactance formula for center-fed dipole
    double X = -120.0 * (std::cos(kl) / std::sin(kl / 2.0)) *
               std::log(std::sin(kl / 2.0) / (ka));

    // Alternative simplified formula
    double kl_half = kl / 2.0;
    double cot_kl_half = std::cos(kl_half) / std::sin(kl_half);

    X = -120.0 * cot_kl_half * std::log(2.0 * lenght / radius);

    return X;
}

double dipoleAntenna::inductance(double f) const {
    double w = 2.0 * PI * f;
    double x = reactance(f);
    if (x > 0) {
        return (x / w) * 1e9;
    }
    return 0.0;
}

std::complex<double> dipoleAntenna::inputZ(double f) const {
    double R = radationResistance(f);
    double X = reactance(f);
    return std::complex<double>(R, X);
}
