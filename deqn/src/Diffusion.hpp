#ifndef DIFFUSION_H_
#define DIFFUSION_H_

#include "InputFile.hpp"
#include "Mesh.hpp"
#include "Scheme.hpp"

#include <vector>

class Diffusion {
private:
    Mesh *mesh;

    Scheme *scheme;

    std::vector<double> subregion;

public:
    Diffusion(const InputFile *input, Mesh *m);

    ~Diffusion();

    void init();
    void doCycle(const double dt);
};
#endif
