#include <string>
#include "datadef.hpp"
#include "matrix.hpp"

/* Given the boundary conditions defined by the flag matrix, update
 * the u and v velocities. Also enforce the boundary conditions at the
 * edges of the matrix.
 */
void apply_boundary_conditions(
    const matrix<float>& u,
    const matrix<float>& v,
    const matrix<char>& flag,
    const int& imax,
    const int& jmax,
    const float& ui,
    const float& vi
) {
    int i, j;

    for (j = 0; j <= jmax + 1; j++) {
        /* Fluid freely flows in from the west */
        u->at(0, j) = u->at(1, j);
        v->at(0, j) = v->at(1, j);

        /* Fluid freely flows out to the east */
        u->at(imax, j) = u->at(imax - 1, j);
        v->at(imax + 1, j) = v->at(imax, j);
    }

    for (i = 0; i <= imax + 1; i++) {
        /* The vertical velocity approaches 0 at the north and south
         * boundaries, but fluid flows freely in the horizontal direction */
        v->at(i, jmax) = 0.0;
        u->at(i, jmax + 1) = u->at(i, jmax);

        v->at(i, 0) = 0.0;
        u->at(i, 0) = u->at(i, 1);
    }

    /* Apply no-slip boundary conditions to cells that are adjacent to
     * internal obstacle cells. This forces the u and v velocity to
     * tend towards zero in these cells.
     */
    for (i = 1; i <= imax; i++) {
        for (j = 1; j <= jmax; j++) {
            if (flag->at(i, j) & B_NSEW) {
                switch (flag->at(i, j)) {
                    case B_N:
                        v->at(i, j) = 0.0;
                        u->at(i, j) = -u->at(i, j + 1);
                        u->at(i - 1, j) = -u->at(i - 1, j + 1);
                        break;
                    case B_E:
                        u->at(i, j) = 0.0;
                        v->at(i, j) = -v->at(i + 1, j);
                        v->at(i, j - 1) = -v->at(i + 1, j - 1);
                        break;
                    case B_S:
                        v->at(i, j - 1) = 0.0;
                        u->at(i, j) = -u->at(i, j - 1);
                        u->at(i - 1, j) = -u->at(i - 1, j - 1);
                        break;
                    case B_W:
                        u->at(i - 1, j) = 0.0;
                        v->at(i, j) = -v->at(i - 1, j);
                        v->at(i, j - 1) = -v->at(i - 1, j - 1);
                        break;
                    case B_NE:
                        v->at(i, j) = 0.0;
                        u->at(i, j) = 0.0;
                        v->at(i, j - 1) = -v->at(i + 1, j - 1);
                        u->at(i - 1, j) = -u->at(i - 1, j + 1);
                        break;
                    case B_SE:
                        v->at(i, j - 1) = 0.0;
                        u->at(i, j) = 0.0;
                        v->at(i, j) = -v->at(i + 1, j);
                        u->at(i - 1, j) = -u->at(i - 1, j - 1);
                        break;
                    case B_SW:
                        v->at(i, j - 1) = 0.0;
                        u->at(i - 1, j) = 0.0;
                        v->at(i, j) = -v->at(i - 1, j);
                        u->at(i, j) = -u->at(i, j - 1);
                        break;
                    case B_NW:
                        v->at(i, j) = 0.0;
                        u->at(i - 1, j) = 0.0;
                        v->at(i, j - 1) = -v->at(i - 1, j - 1);
                        u->at(i, j) = -u->at(i, j + 1);
                        break;
                }
            }
        }
    }

    /* Finally, fix the horizontal velocity at the  western edge to have
     * a continual flow of fluid into the simulation.
     */
    v->at(0, 0) = 2 * vi - v->at(1, 0);
    for (j = 1; j <= jmax; j++) {
        u->at(0, j) = ui;
        v->at(0, j) = 2 * vi - v->at(1, j);
    }
}
