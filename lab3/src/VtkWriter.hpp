#pragma once

#include "Mesh.hpp"

#include <string>

class VtkWriter final {
    std::string dump_basename;

    std::string vtk_header;

    Mesh mesh;

    void writeVtk(const int& step, const double& time);

public:
    VtkWriter(std::string basename, const Mesh& mesh);

    void write(const int& step, const double& time);
};
