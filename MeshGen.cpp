//
// Created by Bartek Siwik on 01/11/2025.
//

#include "MeshGen.h"
#include <cmath>
#include <iostream>

MeshGen::MeshGen(double dipole_length, double dipole_radius,
                             double domain_radius)
    : L(dipole_length), r_dipole(dipole_radius), r_domain(domain_radius)
{
}

void MeshGen::generateMesh(double edge_size) {
    // UWAGA: To uproszczona wersja!
    // Prawdziwa implementacja według PDF wymaga:
    // 1. Persson mesh generator (adaptive Delaunay)
    // 2. Biblioteki typu delaunator-cpp lub CGAL

    std::cout << "Generowanie uproszczonej siatki...\n";

    // Parametry siatki
    int n_radial = 20;    // Warstwy radialne
    int n_angular = 36;   // Punkty kątowe (co 10 stopni)

    // Generuj węzły w układzie cylindrycznym
    nodes.clear();
    int node_id = 0;

    // Węzły wokół dipola (poziomy przekrój)
    for (int ir = 0; ir < n_radial; ir++) {
        double r = r_dipole + (r_domain - r_dipole) * ir / (n_radial - 1);

        int n_angle = (ir == 0) ? 8 : n_angular; // Mniej punktów przy dipolu

        for (int ia = 0; ia < n_angle; ia++) {
            double theta = 2.0 * M_PI * ia / n_angle;
            double x = r * std::cos(theta);
            double y = r * std::sin(theta);

            nodes.emplace_back(x, y, node_id++);
        }
    }

    // Węzeł środkowy (feed point)
    nodes.emplace_back(0.0, 0.0, node_id);
    feed_node = node_id++;

    std::cout << "Wygenerowano " << nodes.size() << " węzłów\n";

    // Generuj elementy (uproszczona triangulacja)
    // W rzeczywistości użyj delaunator-cpp lub CGAL!
    elements.clear();

    // Przykładowa triangulacja koncentryczna
    int n0 = 0;
    for (int ir = 0; ir < n_radial - 1; ir++) {
        int n_curr = (ir == 0) ? 8 : n_angular;
        int n_next = n_angular;

        for (int ia = 0; ia < n_curr; ia++) {
            int i1 = n0 + ia;
            int i2 = n0 + (ia + 1) % n_curr;
            int i3 = n0 + n_curr + ia * n_next / n_curr;
            int i4 = n0 + n_curr + ((ia + 1) * n_next / n_curr) % n_next;

            // Dwa trójkąty na każdy quad
            elements.emplace_back(i1, i2, i3);
            elements.emplace_back(i2, i4, i3);
        }

        n0 += n_curr;
    }

    std::cout << "Wygenerowano " << elements.size() << " elementów\n";
}

std::vector<int> MeshGen::getDirichletNodes() const {
    // Węzły na powierzchni dipola (PEC): r ≈ r_dipole
    std::vector<int> dirichlet_nodes;

    for (const auto& node : nodes) {
        double r = std::sqrt(node.x * node.x + node.y * node.y);
        if (r < r_dipole * 1.1) { // Tolerancja 10%
            dirichlet_nodes.push_back(node.id);
        }
    }

    return dirichlet_nodes;
}

std::vector<int> MeshGen::getABCNodes() const {
    // Węzły na zewnętrznej granicy: r ≈ r_domain
    std::vector<int> abc_nodes;

    for (const auto& node : nodes) {
        double r = std::sqrt(node.x * node.x + node.y * node.y);
        if (r > r_domain * 0.95) { // Tolerancja 5%
            abc_nodes.push_back(node.id);
        }
    }

    return abc_nodes;
}
