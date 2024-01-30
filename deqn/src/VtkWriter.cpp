#include "VtkWriter.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>

VtkWriter::VtkWriter(std::string basename, const std::shared_ptr<Mesh>& mesh)
    : dump_basename{std::move(basename)},
      mesh{mesh} {
    std::ofstream file;
    std::stringstream fname;

    if (!std::filesystem::exists(_path_prefix) || !std::filesystem::is_directory(_path_prefix)) {
        std::filesystem::create_directory(_path_prefix);
    }

    fname << _path_prefix << dump_basename << ".visit";

    file.open(fname.str());

    file << "!NBLOCKS " << 1 << '\n';

    file.close();
}

void VtkWriter::write(const int& step, const double& time) const {
    // Master process writes out the .visit file to coordinate the .vtk files
    std::ofstream file;
    std::stringstream fname;

    fname << _path_prefix << dump_basename << ".visit";

    // Open file in append mode
    file.open(
        fname.str(),
        std::ofstream::out | std::ofstream::in | std::ofstream::app
    );

    file
            << dump_basename
            << "."
            << step
            << "."
            << 1
            << ".vtk" << '\n';

    writeVtk(step, time);

    file.close();
}

void VtkWriter::writeVtk(const int& step, const double& time) const {
    std::ofstream file;
    std::stringstream fname;

    fname
            << _path_prefix
            << dump_basename
            << "."
            << step
            << "."
            << 1
            << ".vtk";

    file.open(fname.str());

    file.setf(std::ios::fixed, std::ios::floatfield);
    file.precision(8);

    file << vtk_header;

    file << "DATASET RECTILINEAR_GRID" << '\n';
    file << "FIELD FieldData 2" << '\n';
    file << "TIME 1 1 double" << '\n';
    file << time << '\n';
    file << "CYCLE 1 1 int" << '\n';
    file << step << '\n';

    if (Mesh::getDim() == 2) {
        file
                << "DIMENSIONS " << mesh->getNx()[0] + 1
                << " " << mesh->getNx()[1] + 1
                << " 1"
                << '\n';
    } else if (Mesh::getDim() == 3) {
        file
                << "DIMENSIONS " << mesh->getNx()[0] + 1
                << " " << mesh->getNx()[1] + 1
                << " " << mesh->getNx()[2] + 1
                << '\n';
    }

    file << "X_COORDINATES " << mesh->getNx()[0] + 1 << " float" << '\n';
    for (int i = 1; i <= mesh->getNx()[0] + 1; i++) {
        file << mesh->getCellX()[i] << " ";
    }

    file << '\n';

    file << "Y_COORDINATES " << mesh->getNx()[1] + 1 << " float" << '\n';
    for (int j = 1; j <= mesh->getNx()[1] + 1; j++) {
        file << mesh->getCellY()[j] << " ";
    }

    file << '\n';

    file << "Z_COORDINATES 1 float" << '\n';
    file << "0.0000" << '\n';
    file << "CELL_DATA " << mesh->getNx()[0] * mesh->getNx()[1] << '\n';
    file << "FIELD FieldData 1" << '\n';
    file << "u 1 " << mesh->getNx()[0] * mesh->getNx()[1] << " double" << '\n';

    for (int j = 1; j <= mesh->getNx()[1]; j++) {
        for (int i = 1; i <= mesh->getNx()[0]; i++) {
            file << mesh->getU0()[i + j * (mesh->getNx()[0] + 2)] << " ";
        }
        file << '\n';
    }

    file.close();
}
