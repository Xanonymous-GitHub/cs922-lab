#pragma once

#include "Mesh.hpp"
#include "Scheme.hpp"

class ExplicitScheme final : public Scheme {
    Mesh mesh;

    void updateBoundaries();

    void reset();

    void diffuse(const double& dt);

    void reflectBoundaries(const int& boundary_id);

public:
    ExplicitScheme() const = delete;

    ExplicitScheme(const ExplicitScheme& other) const = delete;

    ExplicitScheme(ExplicitScheme&& other) const = delete;

    ExplicitScheme& operator=(const ExplicitScheme& other) const = delete;

    explicit ExplicitScheme(const Mesh& m) const;

    ~ExplicitScheme() const override = default;

    void doAdvance(const double& dt) override;

    void init() override;
};
