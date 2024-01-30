#include "ExplicitScheme.hpp"

#include <iostream>

ExplicitScheme::ExplicitScheme(const std::shared_ptr<Mesh>& m) : mesh(m) {
}

void ExplicitScheme::doAdvance(const double& dt) const {
    diffuse(dt);
    reset();
    updateBoundaries();
}

void ExplicitScheme::updateBoundaries() const {
#pragma omp parallel for schedule(static) default(none) shared(mesh)
    for (int i = 0; i < 4; i++) {
        reflectBoundaries(i);
    }
}

void ExplicitScheme::init() const {
    updateBoundaries();
}

void ExplicitScheme::reset() const {
    auto& u0 = mesh->getU0();
    const auto& u1 = mesh->getU1();
    const int x_min = mesh->getMin()[0];
    const int x_max = mesh->getMax()[0];
    const int y_min = mesh->getMin()[1];
    const int y_max = mesh->getMax()[1];

    const int nx = mesh->getNx()[0] + 2;

#pragma omp parallel for collapse(2) schedule(static) default(none) shared(u0, u1, x_min, x_max, y_min, y_max, nx)
    for (int k = y_min - 1; k <= y_max + 1; k++) {
        for (int j = x_min - 1; j <= x_max + 1; j++) {
            const int i = Mesh::poly2(j, k, x_min - 1, y_min - 1, nx);
            u0[i] = u1[i];
        }
    }
}

void ExplicitScheme::diffuse(const double& dt) const {
    const auto& u0 = mesh->getU0();
    auto& u1 = mesh->getU1();
    const int x_min = mesh->getMin()[0];
    const int x_max = mesh->getMax()[0];
    const int y_min = mesh->getMin()[1];
    const int y_max = mesh->getMax()[1];
    const double dx = mesh->getDx()[0];
    const double dy = mesh->getDx()[1];

    const int nx = mesh->getNx()[0] + 2;

    const double rx = dt / (dx * dx);
    const double ry = dt / (dy * dy);

#pragma omp parallel for collapse(2) schedule(static) default(none) shared(u0, u1, x_min, x_max, y_min, y_max, nx, rx, ry)
    for (int k = y_min; k <= y_max; k++) {
        for (int j = x_min; j <= x_max; j++) {
            const int n1 = Mesh::poly2(j, k, x_min - 1, y_min - 1, nx);
            const int n2 = Mesh::poly2(j - 1, k, x_min - 1, y_min - 1, nx);
            const int n3 = Mesh::poly2(j + 1, k, x_min - 1, y_min - 1, nx);
            const int n4 = Mesh::poly2(j, k - 1, x_min - 1, y_min - 1, nx);
            const int n5 = Mesh::poly2(j, k + 1, x_min - 1, y_min - 1, nx);

            u1[n1] = (1.0 - 2.0 * rx - 2.0 * ry) * u0[n1] + rx * u0[n2] + rx * u0[n3] + ry * u0[n4] + ry * u0[n5];
        }
    }
}

void ExplicitScheme::reflectBoundaries(const int& boundary_id) const {
    auto& u0 = mesh->getU0();
    const int x_min = mesh->getMin()[0];
    const int x_max = mesh->getMax()[0];
    const int y_min = mesh->getMin()[1];
    const int y_max = mesh->getMax()[1];

    const int nx = mesh->getNx()[0] + 2;

    switch (boundary_id) {
        case 0: {
            /* top */
#pragma omp parallel for schedule(static) default(none) shared(u0, x_min, x_max, y_min, y_max, nx)
            for (int j = x_min; j <= x_max; j++) {
                const int n1 = Mesh::poly2(j, y_max, x_min - 1, y_min - 1, nx);
                const int n2 = Mesh::poly2(j, y_max + 1, x_min - 1, y_min - 1, nx);

                u0[n2] = u0[n1];
            }
            break;
        }

        case 1: {
            /* right */
#pragma omp parallel for schedule(static) default(none) shared(u0, x_min, x_max, y_min, y_max, nx)
            for (int k = y_min; k <= y_max; k++) {
                const int n1 = Mesh::poly2(x_max, k, x_min - 1, y_min - 1, nx);
                const int n2 = Mesh::poly2(x_max + 1, k, x_min - 1, y_min - 1, nx);

                u0[n2] = u0[n1];
            }
            break;
        }

        case 3: {
            /* left */
#pragma omp parallel for schedule(static) default(none) shared(u0, x_min, x_max, y_min, y_max, nx)
            for (int k = y_min; k <= y_max; k++) {
                const int n1 = Mesh::poly2(x_min, k, x_min - 1, y_min - 1, nx);
                const int n2 = Mesh::poly2(x_min - 1, k, x_min - 1, y_min - 1, nx);

                u0[n2] = u0[n1];
            }
            break;
        }

        case 2: {
            /* bottom */
#pragma omp parallel for schedule(static) default(none) shared(u0, x_min, x_max, y_min, y_max, nx)
            for (int j = x_min; j <= x_max; j++) {
                const int n1 = Mesh::poly2(j, y_min, x_min - 1, y_min - 1, nx);
                const int n2 = Mesh::poly2(j, y_min - 1, x_min - 1, y_min - 1, nx);

                u0[n2] = u0[n1];
            }
            break;
        }

        default:
            std::cerr << "Error in reflectBoundaries(): unknown boundary id (" << boundary_id << ")" << '\n';
    }
}
