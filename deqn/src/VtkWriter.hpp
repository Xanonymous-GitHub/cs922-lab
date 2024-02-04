#ifndef VTKWRITER_H_
#define VTKWRITER_H_

#include "Mesh.hpp"
#include <string>

class VtkWriter {
private:
    std::string dump_basename;

    std::string vtk_header;

    Mesh *mesh;

    inline static const std::string _path_prefix = "out/";

    void writeVtk(int step, double time);

public:
    VtkWriter(std::string basename, Mesh *mesh);

    void write(int step, double time);
};
#endif
