#include <iostream>
#include <iomanip>
#include <complex>
#include <vector>
#include "FEM.h"
#include "ImpedanceMatcher.h"

int main() {
    std::cout << std::fixed << std::setprecision(3);

    double fmin_mhz, fmax_mhz, length_m, radius_m;
    int n_points, nElemHalf;

    std::cout << "Min freq [MHz]: ";
    std::cin >> fmin_mhz;
    std::cout << "Max freq [MHz]: ";
    std::cin >> fmax_mhz;
    double fmin = fmin_mhz * 1e6, fmax = fmax_mhz * 1e6;

    std::cout << "Dipole total length [m]: ";
    std::cin >> length_m;
    std::cout << "Wire radius [m]: ";
    std::cin >> radius_m;

    std::cout << "Number of frequency points: ";
    std::cin >> n_points;
    if (n_points < 2) n_points = 10;

    std::cout << "Elements per half (e.g. 20): ";
    std::cin >> nElemHalf;
    if (nElemHalf < 4) nElemHalf = 10;

    double fstep = (fmax - fmin) / (n_points - 1);
    FEM fem(length_m, radius_m, 2 * nElemHalf);
    ImpedanceMatcher matcher(50.0); // np. 50 Ω system

    std::cout << "\nFreq[MHz]\tRe(Z)[Ω]\tIm(Z)[Ω]\tRL[dB]\n";
    for (int i = 0; i < n_points; i++) {
        double f = fmin + i * fstep;
        std::complex<double> Z = fem.solve_impedance(f);
        double RL = matcher.RL(Z);
        std::cout << f / 1e6 << "\t\t" << Z.real() << "\t\t" << Z.imag()
                  << "\t\t" << RL << "\n";
    }

    std::cout << "\n(Uwaga: to uproszczony 1D-FEM, nie modeluje pełnego promieniowania 3D.)\n";
    return 0;
}
