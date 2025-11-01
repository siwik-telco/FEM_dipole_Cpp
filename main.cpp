
#include <iostream>
#include <iomanip>
#include <complex>

#include "FEM.h"
//#include "FEM2D.h"

int main() {
    std::cout << std::fixed << std::setprecision(3);

    // Parametry anteny
    double fmin_mhz = 350.0;
    double fmax_mhz = 450.0;
    double length_m = 0.356;  // Poprawiona długość dla 400 MHz
    double radius_m = 0.0006045; // 0.6045 mm (jak w PDF dla 2400 MHz)
    int n_points = 11;

    // Promień domeny = 4λ (zalecenie z PDF, Tabela I)
    double f_center = 400e6;
    double lambda = 3e8 / f_center;
    double domain_radius = 4.0 * lambda; // 60 mm dla 2400 MHz w PDF

    std::cout << "=== FEM 2D Simulation - Dipole Antenna ===\n";
    std::cout << "Według: \"Modeling Cylindrical Dipole Antenna by FEM\" (PDF)\n\n";
    std::cout << "Parametry:\n";
    std::cout << "  Długość dipola: " << length_m * 1000 << " mm\n";
    std::cout << "  Promień drutu: " << radius_m * 1000 << " mm\n";
    std::cout << "  Promień domeny: " << domain_radius * 1000 << " mm\n";
    std::cout << "  Zakres częstotliwości: " << fmin_mhz << "-" << fmax_mhz << " MHz\n\n";

    // Inicjalizacja FEM
    FEM2D fem(length_m, radius_m, domain_radius);

    double fstep = (fmax_mhz - fmin_mhz) / (n_points - 1);

    std::cout << "Freq[MHz]\tRe(Z)[Ω]\tIm(Z)[Ω]\t|Z|[Ω]\t\tSWR\n";
    std::cout << std::string(70, '-') << "\n";

    for (int i = 0; i < n_points; i++) {
        double f_mhz = fmin_mhz + i * fstep;
        double f = f_mhz * 1e6;

        std::complex<double> Z = fem.solve_impedance(f);

        // VSWR względem 50 Ω
        double R0 = 50.0;
        std::complex<double> Gamma = (Z - R0) / (Z + R0);
        double VSWR = (1.0 + std::abs(Gamma)) / (1.0 - std::abs(Gamma));

        std::cout << f_mhz << "\t\t"
                  << Z.real() << "\t\t"
                  << Z.imag() << "\t\t"
                  << std::abs(Z) << "\t\t"
                  << VSWR << "\n";
    }

    std::cout << "\n(Symulacja 2D FEM - poziomy przekrój anteny)\n";
    std::cout << "Dla pełnej 3D symulacji potrzebna implementacja pionowego przekroju.\n";

    return 0;
}
