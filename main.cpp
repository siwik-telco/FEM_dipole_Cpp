#include <iostream>
#include <vector>
#include <iomanip>
#include <complex>
#include "dipoleAntenna.h"
#include "ImpedanceMatcher.h"

int main() {
    double f_min_mhz, f_max_mhz, length, width, R_load;
    int n_points;

    std::cout << "Min freq [MHz]: ";
    std::cin >> f_min_mhz;
    std::cout << "Max freq [MHz]: ";
    std::cin >> f_max_mhz;

    double f_min = f_min_mhz * 1e6;
    double f_max = f_max_mhz * 1e6;

    std::cout << "Length of dipole [m]: ";
    std::cin >> length;
    std::cout << "Width of dipole [m]: ";
    std::cin >> width;


    std::cout << "Input impedance of system [Ohm]: ";
    std::cin >> R_load;


    std::cout << "Number of points in simulation: ";
    std::cin >> n_points;
    if (n_points < 2) n_points = 10;

    dipoleAntenna dipole(length, width);
    ImpedanceMatcher matcher(R_load);

    double df = (f_max - f_min) / (n_points - 1.0);

    std::cout << std::fixed << std::setprecision(2);
    std::cout << "\nFreq [MHz]\tR [Ω]\tX [Ω]\tL [nH]\tRL [dB]" << std::endl;

    double minRL = -50;
    double f_best = 0.0;
    std::complex<double> Z_best;

    for (int i = 0; i < n_points; ++i) {
        double f = f_min + i * df;
        double freq_mhz = f / 1e6;
        std::complex<double> Z = dipole.inputZ(f);
        double R = std::real(Z);
        double X = std::imag(Z);
        double L_ind = dipole.inductance(f);
        double rl_db = matcher.RL(Z);

        std::cout << freq_mhz << "\t\t" << R << "\t" << X << "\t" << L_ind << "\t" << rl_db << std::endl;

        if (rl_db < minRL) {
            minRL = rl_db;
            f_best = freq_mhz;
            Z_best = Z;
        }
    }

    std::cout << "\nBest impedance match: " << f_best << " MHz, RL = " << minRL << " dB" << std::endl;
//std::cout << "For RL < -10 dB, change the lenght for ~2-5% jeśli X >0." << std::endl;

    return 0;
}
