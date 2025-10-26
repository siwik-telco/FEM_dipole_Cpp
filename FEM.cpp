
//
// Created by Bartosz on 19.10.2025.
//
#include "FEM.h"
#include <cmath>
#include <iostream>

const double c0 = 299792458.0;
const double mu0 = 4.0 * M_PI * 1e-7;
const double eps0 = 1.0 / (mu0 * c0 * c0);
const std::complex<double> Ij(0.0, 1.0);

std::vector<std::complex<double>> FEM::denseSolver(
        std::vector<std::vector<std::complex<double>>>& A,
        std::vector<std::complex<double>>& b)
{
    int n = (int)A.size();
    for (int k = 0; k < n; k++) {
        int piv = k;
        double best = std::abs(A[k][k]);
        for (int r = k + 1; r < n; r++) {
            double val = std::abs(A[r][k]);
            if (val > best) { best = val; piv = r; }
        }
        if (piv != k) {
            std::swap(A[k], A[piv]);
            std::swap(b[k], b[piv]);
        }
        if (std::abs(A[k][k]) < 1e-18) continue;
        for (int i = k + 1; i < n; i++) {
            std::complex<double> factor = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }
    std::vector<std::complex<double>> x(n, {0.0, 0.0});
    for (int i = n - 1; i >= 0; i--) {
        std::complex<double> s = b[i];
        for (int j = i + 1; j < n; j++) s -= A[i][j] * x[j];
        if (std::abs(A[i][i]) < 1e-18)
            x[i] = {0.0, 0.0};
        else
            x[i] = s / A[i][i];
    }
    return x;
}


void FEM::apply_dirichlet(std::vector<std::vector<std::complex<double>>>& A,
                          std::vector<std::complex<double>>& b,
                          int node,
                          std::complex<double> value)
{
    int n = (int)A.size();
    for (int j = 0; j < n; j++) A[node][j] = {0.0, 0.0};
    for (int i = 0; i < n; i++) A[i][node] = {0.0, 0.0};
    A[node][node] = {1.0, 0.0};
    b[node] = value;
}


// POPRAWKA: FEM zamiast FEMDipole
FEM::FEM(double length_m, double radius_m, int nElems)
    : L(length_m), r(radius_m), Nelem(nElems)
{
    if (Nelem < 2) Nelem = 10;
    Nnodes = Nelem + 1;
}


std::complex<double> FEM::solve_impedance(double f)
{
    double omega = 2.0 * M_PI * f;
    std::complex<double> k = omega / std::complex<double>(c0, 0.0);

    double he = L / Nelem;
    std::vector<std::vector<std::complex<double>>> A(
        Nnodes, std::vector<std::complex<double>>(Nnodes, {0.0, 0.0}));
    std::vector<std::complex<double>> rhs(Nnodes, {0.0, 0.0});

    std::complex<double> coeff_grad = 1.0 / he;
    std::complex<double> coeff_mass = he / 6.0;

    for (int e = 0; e < Nelem; e++) {
        int n1 = e;
        int n2 = e + 1;

        std::complex<double> ksq = k * k;
        std::complex<double> ke00 = coeff_grad * 1.0 - ksq * coeff_mass * 2.0;
        std::complex<double> ke01 = coeff_grad * -1.0 - ksq * coeff_mass * 1.0;
        std::complex<double> ke10 = ke01;
        std::complex<double> ke11 = coeff_grad * 1.0 - ksq * coeff_mass * 2.0;

        A[n1][n1] += ke00;
        A[n1][n2] += ke01;
        A[n2][n1] += ke10;
        A[n2][n2] += ke11;
    }

    int node_center = Nnodes / 2;
    rhs[node_center] = {1.0, 0.0};

    apply_dirichlet(A, rhs, 0, {0.0, 0.0});
    apply_dirichlet(A, rhs, Nnodes - 1, {0.0, 0.0});

    std::vector<std::complex<double>> I = denseSolver(A, rhs);
    std::complex<double> I_feed = I[node_center];
    return 1.0 / I_feed;
}
