#ifndef DRIVER_H_
#define DRIVER_H_

#include <string>

#include "Diffusion.hpp"
#include "InputFile.hpp"
#include "Mesh.hpp"
#include "VtkWriter.hpp"

class Driver {
private:
    InputFile *input;
    Mesh *mesh;
    Diffusion *diffusion;
    VtkWriter *writer;

    double t_start;
    double t_end;
    double dt_max;

    double dt;

    std::string problem_name;

    int vis_frequency;
    int summary_frequency;

public:
    Driver(const InputFile *input, const std::string &problem_name);

    ~Driver();

    void run();
};
#endif
