//
// Created by Bartosz on 17.10.2025.
//
#include <iostream>
#include "DisplayResults.h"
#include <fstream>
#include <cstdlib>
DisplayResults::DisplayResults() {

}

void DisplayResults::display() {
    std::ofstream gnu_file("plot.gnu");
    gnu_file << "set terminal wxt size 1200,800 enhanced font 'Arial,10' persist\n";
    gnu_file << "set multiplot layout 2,2 title 'Results of dipole antenna simulation'\n";
    gnu_file << "set xlabel 'Frequency [MHz]'\n";
    gnu_file << "set ylabel 'RL [dB]'\n";
    gnu_file << "set title 'Return Loss vs Frequency'\n";
    gnu_file << "set grid\n";
    gnu_file << "set yrange [*:0]  # Odwrócona skala Y dla RL (niższe = lepsze)\n";
    gnu_file << "plot 'data.dat' using 1:5 with linespoints lt 1 lw 2 pt 7 title 'RL'\n";  // Kolumna 1: freq, 5: RL

    gnu_file << "set ylabel 'R [Ω]'\n";
    gnu_file << "set title 'Rezystancja R vs Częstotliwość'\n";
    gnu_file << "unset yrange\n";
    gnu_file << "plot 'data.dat' using 1:2 with linespoints lt 2 lw 2 pt 5 title 'R'\n";  // Kolumna 2: R

    gnu_file << "set ylabel 'X [Ω]'\n";
    gnu_file << "set title 'Reaktancja X vs Częstotliwość'\n";
    gnu_file << "plot 'data.dat' using 1:3 with linespoints lt 3 lw 2 pt 6 title 'X'\n";  // Kolumna 3: X

    gnu_file << "set ylabel 'L [nH]'\n";
    gnu_file << "set title 'Indukcyjność L vs Częstotliwość'\n";
    gnu_file << "plot 'data.dat' using 1:4 with linespoints lt 4 lw 2 pt 7 title 'L' every :::0::  # Tylko gdzie L > 0, ale gnuplot filtruje ręcznie\n";  // Kolumna 4: L
    gnu_file << "unset multiplot\n";
    gnu_file << "pause -1  # Czekaj na zamknięcie okna\n";
    gnu_file.close();

    // Wywołaj gnuplot (działa na Linuxie/Windows z gnuplot)
    int result = system("gnuplot plot.gnu");
    if (result != 0) {
        std::cout << "Błąd: Gnuplot nie uruchomiony. Zainstaluj gnuplot i spróbuj ponownie." << std::endl;
    }
}