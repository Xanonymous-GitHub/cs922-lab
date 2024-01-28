#pragma once

#include "Mesh.hpp"
#include "Scheme.hpp"

#include <memory>

class ExplicitScheme final : public Scheme {
    const std::shared_ptr<Mesh> mesh{};

    void updateBoundaries() const;

    void reset() const;

    void diffuse(const double& dt) const;

    void reflectBoundaries(const int& boundary_id) const;

public:
    ExplicitScheme() = delete;

    ExplicitScheme(const ExplicitScheme& other) = delete;

    ExplicitScheme(ExplicitScheme&& other) = delete;

    ExplicitScheme& operator=(const ExplicitScheme& other) const = delete;

    explicit ExplicitScheme(const std::shared_ptr<Mesh>& m);

    ~ExplicitScheme() override = default;

    void doAdvance(const double& dt) const override;

    void init() const override;
};
