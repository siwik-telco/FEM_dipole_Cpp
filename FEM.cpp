#include "FEM.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

// Stałe fizyczne
const double c0 = 299792458.0;
const double mu0 = 4.0 * M_PI * 1e-7;
const double eps0 = 1.0 / (mu0 * c0 * c0);
const std::complex<double> Ij(0.0, 1.0);

FEM2D::FEM2D(double length_m, double radius_m, double domain_radius_m)
    : L(length_m), r(radius_m), r_domain(domain_radius_m),
      mesh(length_m, radius_m, domain_radius_m)
{
    // Generuj siatkę z rozmiarem krawędzi = λ/10 dla 400 MHz
    double lambda_ref = c0 / 400e6;
    double edge_size = lambda_ref / 10.0;

    mesh.generateMesh(edge_size);
}

dcomp FEM2D::solve_impedance(double freq) {
    double k0 = 2.0 * M_PI * freq / c0;

    const auto& nodes = mesh.getNodes();
    int n_nodes = nodes.size();

    // Inicjalizacja macierzy globalnej i wektora
    std::vector<std::vector<dcomp>> K_global(n_nodes,
        std::vector<dcomp>(n_nodes, {0,0}));
    std::vector<dcomp> F_global(n_nodes, {0,0});

    // Składanie macierzy globalnej (assembly)
    assembleGlobalSystem(k0, K_global, F_global);

    // Warunki brzegowe Dirichleta (PEC na powierzchni anteny)
    auto dirichlet_nodes = mesh.getDirichletNodes();
    applyDirichletBC(K_global, F_global, dirichlet_nodes, {0,0});

    // Absorbing Boundary Conditions
    auto abc_nodes = mesh.getABCNodes();
    applyAbsorbingBC(K_global, k0, abc_nodes);

    // Wzbudzenie w punkcie zasilania (feed point)
    int feed = mesh.getFeedNode();
    F_global[feed] += 1.0; // 1V wzbudzenie

    // Rozwiązanie układu równań
    std::vector<dcomp> E_field = solveLinearSystem(K_global, F_global);

    // Obliczenie impedancji wejściowej
    // Z = V / I, gdzie I można przybliżyć z pola E
    dcomp I_feed = {0,0};

    // Prąd z prawa Ampère'a: ∇×H = J + jωεE
    // Uproszczenie: I ≈ sum(E * area) w otoczeniu feed point
    const auto& elements = mesh.getElements();
    for (const auto& elem : elements) {
        bool contains_feed = false;
        for (int n : elem.nodes) {
            if (n == feed) {
                contains_feed = true;
                break;
            }
        }

        if (contains_feed) {
            double area = elem.area(nodes);
            for (int n : elem.nodes) {
                I_feed += E_field[n] * area / 3.0;
            }
        }
    }

    I_feed *= Ij * 2.0 * M_PI * freq * eps0; // jωεE ≈ J

    if (std::abs(I_feed) < 1e-12) {
        return {1e9, 1e9};
    }

    return 1.0 / I_feed; // Z = V/I
}

void FEM2D::assembleGlobalSystem(
    double k0,
    std::vector<std::vector<dcomp>>& K_global,
    std::vector<dcomp>& F_global) const
{
    const auto& nodes = mesh.getNodes();
    const auto& elements = mesh.getElements();

    // Dla każdego elementu: oblicz macierze lokalne i dodaj do globalnej
    for (const auto& elem : elements) {
        std::array<std::array<dcomp, 3>, 3> M_local, T_local;
        elem.computeLocalMatrices(nodes, k0, M_local, T_local);

        // K^e = M^e - k0^2 * T^e (równanie 15 z PDF)
        std::array<std::array<dcomp, 3>, 3> K_local;
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                K_local[i][j] = M_local[i][j] - k0*k0 * T_local[i][j];
            }
        }

        // Wektor siły (wzbudzenie)
        double Jz = 1.0; // Amplituda gęstości prądu
        auto f_local = elem.computeLocalForce(nodes, Jz);

        // Assembly: dodaj do macierzy globalnej
        for (int i = 0; i < 3; i++) {
            int gi = elem.nodes[i]; // Globalny indeks węzła i

            for (int j = 0; j < 3; j++) {
                int gj = elem.nodes[j]; // Globalny indeks węzła j
                K_global[gi][gj] += K_local[i][j];
            }

            F_global[gi] += f_local[i];
        }
    }
}

void FEM2D::applyDirichletBC(
    std::vector<std::vector<dcomp>>& K,
    std::vector<dcomp>& F,
    const std::vector<int>& bc_nodes,
    dcomp value) const
{
    // Równanie (16) z PDF: E_z = 0 na powierzchni PEC
    int n = K.size();

    for (int node : bc_nodes) {
        // Wyzeruj wiersz i kolumnę
        for (int i = 0; i < n; i++) {
            K[node][i] = {0,0};
            K[i][node] = {0,0};
        }

        // Przekątna = 1, prawa strona = wartość BC
        K[node][node] = {1,0};
        F[node] = value;
    }
}

void FEM2D::applyAbsorbingBC(
    std::vector<std::vector<dcomp>>& K,
    double k0,
    const std::vector<int>& abc_nodes) const
{
    // Równanie (17) z PDF: ∂u/∂n + q*u = 0
    // Dla fali wychodzącej: q = -jk0

    const auto& nodes = mesh.getNodes();
    const auto& elements = mesh.getElements();

    dcomp q = -Ij * k0;

    // Dla każdego węzła ABC: dodaj składnik ABC do przekątnej
    for (int node : abc_nodes) {
        // Znajdź elementy zawierające ten węzeł
        double total_edge_length = 0.0;

        for (const auto& elem : elements) {
            bool has_node = false;
            int local_idx = -1;

            for (int i = 0; i < 3; i++) {
                if (elem.nodes[i] == node) {
                    has_node = true;
                    local_idx = i;
                    break;
                }
            }

            if (has_node) {
                // Przybliżenie: długość krawędzi
                int next = (local_idx + 1) % 3;
                const Node& n1 = nodes[elem.nodes[local_idx]];
                const Node& n2 = nodes[elem.nodes[next]];

                double dx = n2.x - n1.x;
                double dy = n2.y - n1.y;
                double edge_len = std::sqrt(dx*dx + dy*dy);

                total_edge_length += edge_len;
            }
        }

        // p_i = q * u * integral (równanie 13)
        K[node][node] += q * total_edge_length / 2.0;
    }
}

std::vector<dcomp> FEM2D::solveLinearSystem(
    std::vector<std::vector<dcomp>>& A,
    std::vector<dcomp>& b) const
{
    int n = A.size();

    // Eliminacja Gaussa z częściowym wyborem elementu głównego
    for (int k = 0; k < n; k++) {
        // Wybór elementu głównego
        int pivot = k;
        for (int i = k + 1; i < n; i++) {
            if (std::abs(A[i][k]) > std::abs(A[pivot][k])) {
                pivot = i;
            }
        }

        if (pivot != k) {
            std::swap(A[k], A[pivot]);
            std::swap(b[k], b[pivot]);
        }

        // Regularyzacja dla małych elementów
        if (std::abs(A[k][k]) < 1e-14) {
            A[k][k] = 1e-10;
        }

        // Eliminacja
        for (int i = k + 1; i < n; i++) {
            dcomp factor = A[i][k] / A[k][k];

            for (int j = k; j < n; j++) {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }

    // Podstawienie wsteczne
    std::vector<dcomp> x(n);
    for (int i = n - 1; i >= 0; i--) {
        dcomp sum = {0,0};
        for (int j = i + 1; j < n; j++) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }

    return x;
}

void FEM2D::exportField(const std::string& filename,
                        const std::vector<dcomp>& field) const
{
    std::ofstream file(filename);
    const auto& nodes = mesh.getNodes();

    file << "x,y,Re(E),Im(E),|E|\n";
    for (size_t i = 0; i < nodes.size(); i++) {
        file << std::setprecision(6)
             << nodes[i].x << ","
             << nodes[i].y << ","
             << field[i].real() << ","
             << field[i].imag() << ","
             << std::abs(field[i]) << "\n";
    }

    file.close();
    std::cout << "Pole wyeksportowane do " << filename << "\n";
}
