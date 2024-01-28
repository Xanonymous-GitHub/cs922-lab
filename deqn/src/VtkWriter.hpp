#pragma once

#include "Mesh.hpp"

#include <memory>
#include <string>

class VtkWriter final {
    const std::string dump_basename;

    const std::shared_ptr<Mesh> mesh;

    inline static const std::string vtk_header = "# vtk DataFile Version 3.0\nvtk output\nASCII\n";
    inline static const std::string _path_prefix = "out/";

    void writeVtk(const int& step, const double& time) const;

public:
    VtkWriter(std::string basename, const std::shared_ptr<Mesh>& mesh);

    void write(const int& step, const double& time) const;
};