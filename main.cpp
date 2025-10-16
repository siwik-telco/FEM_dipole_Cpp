#include <iostream>
#include "dipoleAntenna.h"
#include "ImpedanceMatcher.h"
int main () {
    double f_min, f_max, length, width, R_load;
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
    std::cin >> R_load;

    std::cout << "Number of points in simulation: \n";
    std::cin >> n_points;

    dipoleAntenna(length, width);
    ImpedanceMatcher(R_load);



}