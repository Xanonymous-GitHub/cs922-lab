#pragma once

#include "Diffusion.hpp"
#include "InputFile.hpp"
#include "Mesh.hpp"
#include "VtkWriter.hpp"

#include <string>

class Driver final {
    Mesh mesh;
    Diffusion diffusion;
    VtkWriter writer;

    double t_start, t_end, dt_max, dt;

    std::string _problem_name;

    int vis_frequency, summary_frequency;

public:
    Driver() = delete;

    Driver(const Driver& other) = delete;

    Driver(Driver&& other) = delete;

    Driver& operator=(const Driver& other) const = delete;

    Driver(const InputFile& input, const std::string& problem_name);

    ~Driver() = default;

    void run();
};
