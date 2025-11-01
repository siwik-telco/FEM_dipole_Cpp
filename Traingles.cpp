//
// Created by Bartek Siwik on 01/11/2025.
//

#include "Traingles.h"
// #include "TriangleElement.h"

const double PI = M_PI;
const std::complex<double> Ij(0.0, 1.0);

Triangles::Triangles(int n1, int n2, int n3) {
    nodes[0] = n1;
    nodes[1] = n2;
    nodes[2] = n3;
}

double Triangles::area(const std::vector<Node>& nodeList) const {
    const Node& n1 = nodeList[nodes[0]];
    const Node& n2 = nodeList[nodes[1]];
    const Node& n3 = nodeList[nodes[2]];

    // Wzór: A = 0.5 * |x1(y2-y3) + x2(y3-y1) + x3(y1-y2)|
    return 0.5 * std::abs(
        n1.x * (n2.y - n3.y) +
        n2.x * (n3.y - n1.y) +
        n3.x * (n1.y - n2.y)
    );
}

std::array<double, 3> Triangles::shapeFunc(double xi, double eta) {
    // Funkcje kształtu dla elementu trójkątnego (liniowe)
    // N1 = xi, N2 = eta, N3 = 1 - xi - eta
    return {xi, eta, 1.0 - xi - eta};
}

std::array<std::array<double, 2>, 2> Triangles::jacobian(
    const std::vector<Node>& nodeList) const
{
    const Node& n1 = nodeList[nodes[0]];
    const Node& n2 = nodeList[nodes[1]];
    const Node& n3 = nodeList[nodes[2]];

    // Macierz Jacobiego: J = [dx/dxi  dx/deta]
    //                        [dy/dxi  dy/deta]
    std::array<std::array<double, 2>, 2> J;
    J[0][0] = n2.x - n1.x;  // dx/dxi
    J[0][1] = n3.x - n1.x;  // dx/deta
    J[1][0] = n2.y - n1.y;  // dy/dxi
    J[1][1] = n3.y - n1.y;  // dy/deta

    return J;
}

std::array<std::array<double, 2>, 3> Triangles::shapeFuncDerivatives(
    const std::vector<Node>& nodeList) const
{
    auto J = jacobian(nodeList);
    double detJ = J[0][0] * J[1][1] - J[0][1] * J[1][0];

    if (std::abs(detJ) < 1e-12) {
        throw std::runtime_error("Singular Jacobian - degenerate triangle");
    }

    // Odwrotność Jacobiego
    std::array<std::array<double, 2>, 2> Jinv;
    Jinv[0][0] =  J[1][1] / detJ;
    Jinv[0][1] = -J[0][1] / detJ;
    Jinv[1][0] = -J[1][0] / detJ;
    Jinv[1][1] =  J[0][0] / detJ;

    // Pochodne funkcji kształtu w współrzędnych naturalnych
    // dN/dxi = [1, 0, -1], dN/deta = [0, 1, -1]
    std::array<std::array<double, 2>, 3> dN_dxy;

    for (int i = 0; i < 3; i++) {
        double dN_dxi  = (i == 0) ? 1.0 : ((i == 2) ? -1.0 : 0.0);
        double dN_deta = (i == 1) ? 1.0 : ((i == 2) ? -1.0 : 0.0);

        // [dN/dx] = [Jinv] * [dN/dxi ]
        // [dN/dy]           [dN/deta]
        dN_dxy[i][0] = Jinv[0][0] * dN_dxi + Jinv[0][1] * dN_deta;
        dN_dxy[i][1] = Jinv[1][0] * dN_dxi + Jinv[1][1] * dN_deta;
    }

    return dN_dxy;
}

void Triangles::computeLocalMatrices(
    const std::vector<Node>& nodeList,
    double k0,
    std::array<std::array<dcomp, 3>, 3>& M_local,
    std::array<std::array<dcomp, 3>, 3>& T_local) const
{
    double A = area(nodeList);
    auto dN = shapeFuncDerivatives(nodeList);

    // Macierz M_ij według równania (10) z PDF:
    // M_ij = integral{ (dNi/dx * dNj/dx + dNi/dy * dNj/dy) dxdy }
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            M_local[i][j] = A * (
                dN[i][0] * dN[j][0] +
                dN[i][1] * dN[j][1]
            );
        }
    }

    // Macierz T_ij według równania (11) z PDF:
    // T_ij = integral{ Ni * Nj dxdy }
    // Dla elementów liniowych: integral po elemencie
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            if (i == j) {
                T_local[i][j] = A / 6.0;  // diagonal
            } else {
                T_local[i][j] = A / 12.0; // off-diagonal
            }
        }
    }
}

std::array<dcomp, 3> Triangles::computeLocalForce(
    const std::vector<Node>& nodeList,
    double Jz_magnitude) const
{
    double A = area(nodeList);

    // f_i = integral{ Jz * Ni dxdy } (równanie 12)
    // Dla stałego Jz i liniowych funkcji kształtu
    std::array<dcomp, 3> f_local;
    for (int i = 0; i < 3; i++) {
        f_local[i] = Jz_magnitude * A / 3.0; // równomierne rozłożenie
    }

    return f_local;
}
