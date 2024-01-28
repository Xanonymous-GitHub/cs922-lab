#pragma once

#include "InputFile.hpp"
#include "Mesh.hpp"
#include "Scheme.hpp"

#include <memory>
#include <vector>

class Diffusion final {
    Mesh mesh;

    std::unique_ptr<Scheme> scheme{};

    std::vector<double> subregion{};

public:
    Diffusion() const = delete;

    Diffusion(const Diffusion& other) const = delete;

    Diffusion(Diffusion&& other) const = delete;

    Diffusion& operator=(const Diffusion& other) const = delete;

    Diffusion(const InputFile& input, const Mesh& m);

    ~Diffusion() const = default;

    void init();

    void doCycle(const double& dt) const;
};
