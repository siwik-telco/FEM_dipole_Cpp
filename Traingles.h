//
// Created by Bartek Siwik on 01/11/2025.
//

#ifndef TRIANGLEELEMENT_H
#define TRIANGLEELEMENT_H

#include <complex>
#include <array>
#include <vector>
#include <cmath>

using dcomp = std::complex<double>;

// Węzeł 2D (x, y)
struct Node {
    double x, y;
    int id;

    Node(double x_ = 0, double y_ = 0, int id_ = -1)
        : x(x_), y(y_), id(id_) {}
};

// Element trójkątny z 3 węzłami
class Triangles {
public:
    std::array<int, 3> nodes; // Indeksy węzłów (globalne)

    Triangles(int n1, int n2, int n3);

    // Obliczenie pola powierzchni
    double area(const std::vector<Node>& nodeList) const;

    // Funkcje kształtu N_i (shape functions) w punkcie (xi, eta)
    // Współrzędne naturalne: N1 = xi, N2 = eta, N3 = 1-xi-eta
    static std::array<double, 3> shapeFunc(double xi, double eta);

    // Pochodne funkcji kształtu: dN/dx, dN/dy
    std::array<std::array<double, 2>, 3> shapeFuncDerivatives(
        const std::vector<Node>& nodeList) const;

    // Macierze lokalne M_ij i T_ij według PDF (równania 10-11)
    void computeLocalMatrices(
        const std::vector<Node>& nodeList,
        double k0,
        std::array<std::array<dcomp, 3>, 3>& M_local,
        std::array<std::array<dcomp, 3>, 3>& T_local) const;

    // Wektor lokalny f_i dla wzbudzenia (równanie 12)
    std::array<dcomp, 3> computeLocalForce(
        const std::vector<Node>& nodeList,
        double Jz_magnitude) const;

    // Macierz Jacobiego dla transformacji współrzędnych
    std::array<std::array<double, 2>, 2> jacobian(
        const std::vector<Node>& nodeList) const;
};

#endif // TRIANGLEELEMENT_H
