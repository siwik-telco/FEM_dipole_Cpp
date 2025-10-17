#include <iostream>
#include "dipoleAntenna.h"
#include "ImpedanceMatcher.h"
int main () {
    double f_min, f_max, length, width, R_load1;
    int n_points;

    std::cout << "Min freq: \n";
    std::cin >> f_min;

    std::cout << "Max freq: \n";
    std::cin >> f_max;

    std::cout << "Length of dipole: \n";
    std::cin >> length;

    std::cout << "Width of dipole: \n";
    std::cin >> width;

    std::cout << "Imput impendance of system: \n";
    std::cin >> R_load1;

    std::cout << "Number of points in simulation: \n";
    std::cin >> n_points;

    dipoleAntenna dipole(length, width);
    ImpedanceMatcher matcher(R_load1);

    double df =  (f_max-f_min) / (n_points - 1.0);
    std::cout << "test";

    double minRL = -50;
    double fbest = 0.0;
    std::complex<double> Z_best;

    for (int i=0; i<n_points;i++) {
        double f = f_min + i * df;
        double freq = f/1e6;
        std::complex<double> Z =dipole.inputZ(f);
        double R = std::real(Z);
        double X = std::imag(Z);
        double L = dipole.inductance(f);
        double gam = matcher.RL(Z);
        std::cout << "RL = " << gam << ".\n";
        // std::cout << "R: X: L: " << R << " " << X << " " << L << ".\n";
    }
    return 0;
}