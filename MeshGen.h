//
// Created by Bartek Siwik on 01/11/2025.
//

#ifndef MESHGEN_H
#define MESHGEN_H

#include "Traingles.h"
#include <vector>
#include <set>

class MeshGen {
public:
    MeshGen
    (double dipole_length, double dipole_radius,
                  double domain_radius);

    void generateMesh(double edge_size);

    const std::vector<Node>& getNodes() const { return nodes; }
    const std::vector<Triangles>& getElements() const { return elements; }

    // Identyfikacja węzłów z warunkami brzegowymi
    std::vector<int> getDirichletNodes() const; // PEC (na powierzchni anteny)
    std::vector<int> getABCNodes() const;       // ABC (zewnętrzna granica)

    int getFeedNode() const { return feed_node; } // Węzeł wzbudzenia

private:
    double L;           // Długość dipola
    double r_dipole;    // Promień drutu
    double r_domain;    // Promień domeny obliczeniowej

    std::vector<Node> nodes;
    std::vector<Triangles> elements;

    int feed_node;      // Węzeł środkowy (feed point)

    void generateCircularMesh(double edge_size);
    void identifyBoundaryNodes();
};

#endif
