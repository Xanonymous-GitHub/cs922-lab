#pragma once

#include "Mesh.hpp"

#include <memory>
#include <string>

class VtkWriter final {
    std::string dump_basename;

    std::string vtk_header;

    std::shared_ptr<Mesh> mesh;

    void writeVtk(const int& step, const double& time) const;

public:
    VtkWriter(std::string basename, const std::shared_ptr<Mesh>& mesh);

    void write(const int& step, const double& time);
};
