// //
// // Created by Bartosz on 19.10.2025.
// //
//

#ifndef FEM_H
#define FEM_H


#include "Traingles.h"
#include "MeshGen.h"
#include <vector>
#include <complex>

using dcomp = std::complex<double>;

class FEM2D {
public:
    FEM2D(double length_m, double radius_m, double domain_radius_m);

    dcomp solve_impedance(double freq);

    // Eksport wyników do pliku
    void exportField(const std::string& filename,
                     const std::vector<dcomp>& field) const;

private:
    MeshGen mesh;
    double L, r, r_domain;

    // Składanie macierzy globalnej (assembly)
    void assembleGlobalSystem(
        double k0,
        std::vector<std::vector<dcomp>>& K_global,
        std::vector<dcomp>& F_global) const;

    // Nałożenie warunków brzegowych
    void applyDirichletBC(
        std::vector<std::vector<dcomp>>& K,
        std::vector<dcomp>& F,
        const std::vector<int>& bc_nodes,
        dcomp value) const;

    void applyAbsorbingBC(
        std::vector<std::vector<dcomp>>& K,
        double k0,
        const std::vector<int>& abc_nodes) const;

    // Solver układu równań (Gauss)
    std::vector<dcomp> solveLinearSystem(
        std::vector<std::vector<dcomp>>& A,
        std::vector<dcomp>& b) const;
};

#endif // FEM_H
