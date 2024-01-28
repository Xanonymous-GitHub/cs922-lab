#pragma once

#include <array>
#include <vector>
#include "InputFile.hpp"

#define NDIM 2

class Mesh final {
    std::vector<double> u1{}, u0{};
    std::vector<double> cellx{}, celly{};

    std::array<double, NDIM> min_coords{}, max_coords{};
    std::array<int, NDIM> n{}, min{}, max{};
    std::array<double, NDIM> dx{};

    /**
     * A mesh has four neighbours, and they are
     * accessed in the following order:
     * - top
     * - right
     * - bottom
     * - left
     */
    std::array<int, 4> neighbours{};

    void allocate();

    bool allocated = false;

public:
    Mesh() const = delete;

    Mesh(const Mesh& other) const = delete;

    Mesh(Mesh&& other) const = delete;

    Mesh& operator=(const Mesh& other) const = delete;

    explicit Mesh(const InputFile& inputFile);

    ~Mesh() const = default;

    [[nodiscard]]
    const std::vector<double>& getU0() const;

    [[nodiscard]]
    const std::vector<double>& getU1() const;

    [[nodiscard]]
    const std::array<double, NDIM>& getDx() const;

    [[nodiscard]]
    const std::array<int, NDIM>& getNx() const;

    [[nodiscard]]
    const std::array<int, NDIM>& getMin() const;

    [[nodiscard]]
    const std::array<int, NDIM>& getMax() const;

    [[nodiscard]]
    static int getDim();

    [[nodiscard]]
    const std::vector<double>& getCellX() const;

    [[nodiscard]]
    const std::vector<double>& getCellY() const;

    [[nodiscard]]
    const std::array<int, 4>& getNeighbours() const;

    [[nodiscard]]
    double getTotalTemperature() const;

    [[nodiscard]]
    static constexpr int poly2(
        const int& i,
        const int& j,
        const int& imin,
        const int& jmin,
        const int& ni
    ) {
        return i - imin + (j - jmin) * ni;
    }
};
