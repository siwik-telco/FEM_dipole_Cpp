//
// Created by Bartosz on 19.10.2025.
//
#include <cmath>
#include <iostream>
#include "FEM.h"


const double c0 = 299792458.0;
const double mu0 = 4.0 * M_PI * 1e-7;
const double eps0 = 1.0 / (mu0 * c0 * c0);
const std::complex<double> Ij(0.0, 1.0);

std::vector<dcomp> FEM::denseSolver(std::vector<std::vector<dcomp> > &A, std::vector<dcomp> &B) {
    int n = (int)A.size();

    for (int k=0; k< n; k++) {
        int piv = k;
        double best = std::abs(A[k][k]);
        for (int r = k+1; r < n; r++) {
            double val = std::abs(A[k][r]);
            if (val >  best ) {
                best = val;
                piv = r;
            }

            if (piv != k) {
                std::swap(A[k], A[piv]);
                std::swap(B[k], B[piv]);
            }
            if ((std::abs(A[k][k])) < 1e-18) {continue;}
            for (int i= k+1; i < n; i++) {
                std::complex<double> factor = A[i][k] / A[k][k];
                for (int j= k+1; j < n; j++) {
                    A[i][j] -=  factor * A[k][j];
                }
                B[i] -= factor * B[k];
            }
        }
        std::vector<std::complex<double> > x(n, {0.0, 0.0});

        for (int i=n-1; i>=0; i--) {
            std::complex<double> s = B[i];
            for (int j=i+1; j < n; j++) {
                s -= A[i][j] * x[j];

                if (std::abs(A[i][j])< 1e-18) {
                    x[i] = {0.0, 0.0};
                }
                else {
                    x[i] = s/ A[i][j];

                }
                return x;
            }
        }
    }
}

void FEM::apply_dirichlet(std::vector<std::vector<dcomp> > &A, std::vector<dcomp> &B, int node, dcomp value) {

}

