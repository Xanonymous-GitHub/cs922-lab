#ifndef EXPLICIT_SCHEME_H_
#define EXPLICIT_SCHEME_H_

#include "InputFile.h"
#include "Scheme.h"

class ExplicitScheme : public Scheme {
private:
    Mesh *mesh;

    void updateBoundaries();
    void reset();
    void diffuse(double dt);
    void reflectBoundaries(int boundary_id);

public:
    ExplicitScheme(const InputFile *input, Mesh *m);

    void doAdvance(const double dt);

    void init();
};
#endif
