// TODO: remove the following rule disable and fix the issue
// ReSharper disable CppDFAUnusedValue CppDFAUnreadVariable CppDFALoopConditionNotUpdated
#include "Diffusion.hpp"
#include "ExplicitScheme.hpp"

#include <iostream>

Diffusion::Diffusion(const InputFile& input, const std::shared_ptr<Mesh>& m) : mesh{m} {
    if (const std::string scheme_str = input.getString("scheme", "explicit"); scheme_str == "explicit") {
        scheme = std::make_unique<ExplicitScheme>(mesh);
    } else {
        std::cerr << "Error: unknown scheme \"" << scheme_str << "\"" << '\n';
        std::exit(1);
    }

    subregion = input.getDoubleList("subregion", std::vector<double>{});

    if (!subregion.empty() && subregion.size() != 4) {
        std::cerr << "Error:  region must have 4 entries (xmin, ymin, xmax, ymax)" << '\n';
        std::exit(1);
    }

    init();
}

void Diffusion::init() const {
    auto& u0 = mesh->getU0();

    const int x_max = mesh->getNx()[0];
    const int y_max = mesh->getNx()[1];

    const auto& cellx = mesh->getCellX();
    const auto& celly = mesh->getCellY();

    const int nx = x_max + 2;

#pragma omp parallel default(none) shared(u0, x_max, y_max, cellx, celly, subregion, nx)
    {
#pragma omp for collapse(2) schedule(static) nowait
        for (int j = 0; j < y_max + 2; j++) {
            for (int i = 0; i < x_max + 2; i++) {
#pragma omp critical
                {
                    if (
                        !subregion.empty() &&
                        celly[j] > subregion[1] && celly[j] <= subregion[3] &&
                        cellx[i] > subregion[0] && cellx[i] <= subregion[2]
                    ) {
                        u0[i + j * nx] = 10.0;
                    } else {
                        u0[i + j * nx] = 0.0;
                    }
                }
            }
        }
    }

    scheme->init();
}

void Diffusion::doCycle(const double& dt) const {
    scheme->doAdvance(dt);
}
