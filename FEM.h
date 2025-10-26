// //
// // Created by Bartosz on 19.10.2025.
// //
//

#ifndef FEM_H
#define FEM_H
#include <vector>
#include <complex>

using dcomp = std::complex<double>;

class FEM {
public:
    FEM(double length_m, double radius_m, int nElems);
    dcomp solve_impedance(double f);

private:
    double L;
    double r;
    int Nelem;
    int Nnodes;

    void apply_dirichlet(std::vector<std::vector<dcomp>>& A, std::vector<dcomp>& B, int node, dcomp value);
    std::vector<dcomp> denseSolver(std::vector<std::vector<dcomp>>& A, std::vector<dcomp>& B);
};

#endif //FEM_H
