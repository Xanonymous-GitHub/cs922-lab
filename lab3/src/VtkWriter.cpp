#include "VtkWriter.hpp"

#include <fstream>
#include <iostream>
#include <sstream>

VtkWriter::VtkWriter(std::string basename, const Mesh& mesh)
    : dump_basename(std::move(basename)),
      vtk_header("# vtk DataFile Version 3.0\nvtk output\nASCII\n"),
      mesh(mesh) {
    std::ofstream file;
    std::stringstream fname;

    fname << dump_basename << ".visit";

    std::string file_name = fname.str();

    file.open(file_name.c_str());

    file << "!NBLOCKS " << 1 << std::endl;

    file.close();
}

void VtkWriter::write(const int& step, const double& time) {
    // Master process writes out the .visit file to coordinate the .vtk files
    std::ofstream file;
    std::stringstream fname;

    fname << dump_basename << ".visit";

    const std::string file_name = fname.str();

    // Open file in append mode
    file.open(file_name.c_str(), std::ofstream::out | std::ofstream::in | std::ofstream::app);

    file << dump_basename
            << "."
            << step
            << "."
            << 1
            << ".vtk" << std::endl;

    writeVtk(step, time);

    file.close();
}

void VtkWriter::writeVtk(const int& step, const double& time) {
    std::ofstream file;
    std::stringstream fname;

    fname << dump_basename
            << "."
            << step
            << "."
            << 1
            << ".vtk";

    const std::string file_name = fname.str();

    file.open(file_name.c_str());

    file.setf(std::ios::fixed, std::ios::floatfield);
    file.precision(8);

    file << vtk_header;

    file << "DATASET RECTILINEAR_GRID" << std::endl;
    file << "FIELD FieldData 2" << std::endl;
    file << "TIME 1 1 double" << std::endl;
    file << time << std::endl;
    file << "CYCLE 1 1 int" << std::endl;
    file << step << std::endl;

    if (Mesh::getDim() == 2) {
        file
                << "DIMENSIONS " << mesh.getNx()[0] + 1
                << " " << mesh.getNx()[1] + 1
                << " 1"
                << std::endl;
    } else if (Mesh::getDim() == 3) {
        file
                << "DIMENSIONS " << mesh.getNx()[0] + 1
                << " " << mesh.getNx()[1] + 1
                << " " << mesh.getNx()[2] + 1
                << std::endl;
    }

    file << "X_COORDINATES " << mesh.getNx()[0] + 1 << " float" << std::endl;
    for (int i = 1; i <= mesh.getNx()[0] + 1; i++) {
        file << mesh.getCellX()[i] << " ";
    }

    file << std::endl;

    file << "Y_COORDINATES " << mesh.getNx()[1] + 1 << " float" << std::endl;
    for (int j = 1; j <= mesh.getNx()[1] + 1; j++) {
        file << mesh.getCellY()[j] << " ";
    }

    file << std::endl;

    file << "Z_COORDINATES 1 float" << std::endl;
    file << "0.0000" << std::endl;
    file << "CELL_DATA " << mesh.getNx()[0] * mesh.getNx()[1] << std::endl;
    file << "FIELD FieldData 1" << std::endl;
    file << "u 1 " << mesh.getNx()[0] * mesh.getNx()[1] << " double" << std::endl;

    for (int j = 1; j <= mesh.getNx()[1]; j++) {
        for (int i = 1; i <= mesh.getNx()[0]; i++) {
            file << mesh.getU0()[i + j * (mesh.getNx()[0] + 2)] << " ";
        }
        file << std::endl;
    }

    file.close();
}
