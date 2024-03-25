#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "datadef.h"
#include "init.h"

#define max(x, y) ((x)>(y)?(x):(y))
#define min(x, y) ((x)<(y)?(x):(y))

extern int* ileft, * iright;
extern int nprocs, proc;

/* Computation of tentative velocity field (f, g) */
void compute_tentative_velocity(
    float** u,
    float** v,
    float** f,
    float** g,
    char** flag,
    int imax,
    int jmax,
    float del_t,
    float delx,
    float dely,
    float gamma,
    float Re
) {
    register float du2dx, duvdy, duvdx, dv2dy, laplu, laplv;

    #pragma omp parallel
    {
        #pragma omp for collapse(2) private(du2dx, duvdy, laplu) schedule(static) nowait
        for (register int i = 1; i <= imax - 1; i++) {
            for (register int j = 1; j <= jmax; j++) {
                /* only if both adjacent cells are fluid cells */
                if ((flag[i][j] & C_F) && (flag[i + 1][j] & C_F)) {
                    du2dx = ((u[i][j] + u[i + 1][j]) * (u[i][j] + u[i + 1][j]) +
                            gamma * fabs(u[i][j] + u[i + 1][j]) * (u[i][j] - u[i + 1][j]) -
                            (u[i - 1][j] + u[i][j]) * (u[i - 1][j] + u[i][j]) -
                            gamma * fabs(u[i - 1][j] + u[i][j]) * (u[i - 1][j] - u[i][j]))
                            / (4.0 * delx);
                    duvdy = ((v[i][j] + v[i + 1][j]) * (u[i][j] + u[i][j + 1]) +
                            gamma * fabs(v[i][j] + v[i + 1][j]) * (u[i][j] - u[i][j + 1]) -
                            (v[i][j - 1] + v[i + 1][j - 1]) * (u[i][j - 1] + u[i][j]) -
                            gamma * fabs(v[i][j - 1] + v[i + 1][j - 1]) * (u[i][j - 1] - u[i][j]))
                            / (4.0 * dely);
                    laplu = (u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]) / delx / delx +
                            (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]) / dely / dely;

                    f[i][j] = u[i][j] + del_t * (laplu / Re - du2dx - duvdy);
                } else {
                    f[i][j] = u[i][j];
                }
            }
        }

        #pragma omp for collapse(2) private(duvdx, dv2dy, laplv) schedule(static) nowait
        for (register int i = 1; i <= imax; i++) {
            for (register int j = 1; j <= jmax - 1; j++) {
                /* only if both adjacent cells are fluid cells */
                if ((flag[i][j] & C_F) && (flag[i][j + 1] & C_F)) {
                    duvdx = ((u[i][j] + u[i][j + 1]) * (v[i][j] + v[i + 1][j]) +
                            gamma * fabs(u[i][j] + u[i][j + 1]) * (v[i][j] - v[i + 1][j]) -
                            (u[i - 1][j] + u[i - 1][j + 1]) * (v[i - 1][j] + v[i][j]) -
                            gamma * fabs(u[i - 1][j] + u[i - 1][j + 1]) * (v[i - 1][j] - v[i][j]))
                            / (4.0 * delx);
                    dv2dy = ((v[i][j] + v[i][j + 1]) * (v[i][j] + v[i][j + 1]) +
                            gamma * fabs(v[i][j] + v[i][j + 1]) * (v[i][j] - v[i][j + 1]) -
                            (v[i][j - 1] + v[i][j]) * (v[i][j - 1] + v[i][j]) -
                            gamma * fabs(v[i][j - 1] + v[i][j]) * (v[i][j - 1] - v[i][j]))
                            / (4.0 * dely);

                    laplv = (v[i + 1][j] - 2.0 * v[i][j] + v[i - 1][j]) / delx / delx +
                            (v[i][j + 1] - 2.0 * v[i][j] + v[i][j - 1]) / dely / dely;

                    g[i][j] = v[i][j] + del_t * (laplv / Re - duvdx - dv2dy);
                } else {
                    g[i][j] = v[i][j];
                }
            }
        }
    }

    #pragma omp parallel shared(f, g, u, v, imax, jmax)
    {
        #pragma omp for schedule(static) nowait
        for (register int j = 1; j <= jmax; j++) {
            f[0][j] = u[0][j];
            f[imax][j] = u[imax][j];
        }

        #pragma omp for schedule(static) nowait
        for (register int i = 1; i <= imax; i++) {
            g[i][0] = v[i][0];
            g[i][jmax] = v[i][jmax];
        }
    }
}


/* Calculate the right hand side of the pressure equation */
void compute_rhs(
    float** f,
    float** g,
    float** rhs,
    char** flag,
    int imax,
    int jmax,
    float del_t,
    float delx,
    float dely
) {
    #pragma omp parallel for collapse(2) schedule(static)
    for (register int i = 1; i <= imax; i++) {
        for (register int j = 1; j <= jmax; j++) {
            if (flag[i][j] & C_F) {
                /* only for fluid and non-surface cells */
                rhs[i][j] = (
                    (f[i][j] - f[i - 1][j]) / delx +
                    (g[i][j] - g[i][j - 1]) / dely
                ) / del_t;
            }
        }
    }
}

void _red_black_sor(
    register int i,
    register int j,
    register float** p,
    register float** rhs,
    register char** flag,
    int imax,
    int jmax,
    register float omega,
    register float rdx2,
    register float rdy2,
    register float beta_2,
    float pre_calculated_eps_Es[imax + 1][jmax + 1],
    float pre_calculated_eps_Ws[imax + 1][jmax + 1],
    float pre_calculated_eps_Ns[imax + 1][jmax + 1],
    float pre_calculated_eps_Ss[imax + 1][jmax + 1],
    float pre_calculated_beta_mods[imax + 1][jmax + 1]
) {
    if (flag[i][j] == (C_F | B_NSEW)) {
        /* five point star for interior fluid cells */
        p[i][j] = (1. - omega) * p[i][j] -
                    beta_2 * (
                        (p[i + 1][j] + p[i - 1][j]) * rdx2
                        + (p[i][j + 1] + p[i][j - 1]) * rdy2
                        - rhs[i][j]
                    );
    } else if (flag[i][j] & C_F) {
        /* modified star near boundary */
        register const float _eps_E = pre_calculated_eps_Es[i][j];
        register const float _eps_W = pre_calculated_eps_Ws[i][j];
        register const float _eps_N = pre_calculated_eps_Ns[i][j];
        register const float _eps_S = pre_calculated_eps_Ss[i][j];

        p[i][j] = (1. - omega) * p[i][j] -
                    pre_calculated_beta_mods[i][j] * (
                        (_eps_E * p[i + 1][j] + _eps_W * p[i - 1][j]) * rdx2
                        + (_eps_N * p[i][j + 1] + _eps_S * p[i][j - 1]) * rdy2
                        - rhs[i][j]
                    );
    }
}

float _calculate_p0(
    register float** p,
    register char** flag,
    int imax,
    int jmax,
    register int ifull
) {
    register float p0 = 0.0;

    /* Calculate sum of squares */
    #pragma omp parallel for reduction(+:p0) schedule(static) collapse(2)
    for (register int i = 1; i <= imax; i++) {
        for (register int j = 1; j <= jmax; j++) {
            if (flag[i][j] & C_F) { p0 += p[i][j] * p[i][j]; }
        }
    }

    p0 = sqrt(p0 / ifull);
    if (p0 < 0.0001) { p0 = 1.0; }

    return p0;
}

float _calculate_partial_residual(
    register int i,
    register int j,
    register float** p,
    register float** rhs,
    register char** flag,
    register float rdx2,
    register float rdy2,
    int imax,
    int jmax,
    float pre_calculated_eps_Es[imax + 1][jmax + 1],
    float pre_calculated_eps_Ws[imax + 1][jmax + 1],
    float pre_calculated_eps_Ns[imax + 1][jmax + 1],
    float pre_calculated_eps_Ss[imax + 1][jmax + 1]
) {
    // only fluid cells
    if (!(flag[i][j] & C_F)) {
        return 0.0;        
    }

    register const float _eps_E = pre_calculated_eps_Es[i][j];
    register const float _eps_W = pre_calculated_eps_Ws[i][j];
    register const float _eps_N = pre_calculated_eps_Ns[i][j];
    register const float _eps_S = pre_calculated_eps_Ss[i][j];

    register const float add = (_eps_E * (p[i + 1][j] - p[i][j]) -
            _eps_W * (p[i][j] - p[i - 1][j])) * rdx2 +
            (_eps_N * (p[i][j + 1] - p[i][j]) -
            _eps_S * (p[i][j] - p[i][j - 1])) * rdy2 - rhs[i][j];
    return add * add;
}


/* Red/Black SOR to solve the poisson equation */
int poisson(
    register float** p,
    register float** rhs,
    register char** flag,
    int imax,
    int jmax,
    int local_imax,
    int local_jmax,
    int istart,
    int jstart,
    register float eps,
    register int itermax,
    register float omega,
    register int ifull,
    float pre_calculated_eps_Es[imax + 1][jmax + 1],
    float pre_calculated_eps_Ws[imax + 1][jmax + 1],
    float pre_calculated_eps_Ns[imax + 1][jmax + 1],
    float pre_calculated_eps_Ss[imax + 1][jmax + 1],
    register float rdx2,
    register float rdy2,
    register float beta_2,
    float pre_calculated_beta_mods[imax + 1][jmax + 1],
    float p0,
    int west,
    int east,
    int north,
    int south,
    MPI_Comm grid_comm,
    MPI_Win win
) {
    MPI_Datatype halo_row, halo_col;

    MPI_Type_vector(local_jmax, 1, imax + 2, MPI_FLOAT, &halo_col);
    MPI_Type_contiguous(local_imax, MPI_FLOAT, &halo_row);
    MPI_Type_commit(&halo_col);
    MPI_Type_commit(&halo_row);

    register int iter = 0;

    /* Red/Black SOR-iteration */
    for (iter = 0; iter < itermax; iter++) {

        for (register int i = istart; i < istart + local_imax; i++) {
            for (register int j = jstart + 1 - (i % 2); j < jstart + local_jmax; j += 2) {
                _red_black_sor(
                    i, j, p, rhs, flag, imax, jmax, omega, 
                    rdx2, rdy2, beta_2, 
                    pre_calculated_eps_Es, 
                    pre_calculated_eps_Ws, 
                    pre_calculated_eps_Ns, 
                    pre_calculated_eps_Ss,
                    pre_calculated_beta_mods
                );
            }
        }

        MPI_Win_fence(0, win);

        for (register int i = istart; i < istart + local_imax; i += 1) {
            for (register int j = jstart + (i % 2); j < jstart + local_jmax; j += 2) {
                _red_black_sor(
                    i, j, p, rhs, flag, imax, jmax, omega, 
                    rdx2, rdy2, beta_2, 
                    pre_calculated_eps_Es, 
                    pre_calculated_eps_Ws, 
                    pre_calculated_eps_Ns, 
                    pre_calculated_eps_Ss,
                    pre_calculated_beta_mods
                );
            }
        }

        // if(proc != nprocs - 1) {
        //     MPI_Send(&p[*iright][0], jmax + 2, MPI_FLOAT, proc + 1, tag, MPI_COMM_WORLD); 
        // }

        // if(proc > 0) { 
        //     MPI_Recv(&p[*ileft - 1][0], jmax + 2, MPI_FLOAT, proc - 1, tag, MPI_COMM_WORLD, &status);
        // }

        // if(proc > 0) { 
        //     MPI_Send(&p[*ileft][0], jmax + 2, MPI_FLOAT, proc - 1, tag, MPI_COMM_WORLD);
        // }

        // if(proc != nprocs - 1) { 
        //     MPI_Recv(&p[*iright + 1][0], jmax + 2, MPI_FLOAT, proc + 1, tag, MPI_COMM_WORLD, &status);
        // }

        /* Partial computation of residual */
        float local_res = 0.0, global_res = 0.0;;

        // SLOWER: #pragma omp parallel for reduction(+:local_res) schedule(static) collapse(2)
        for (register int i = istart; i < istart + local_imax; i++) {
            for (register int j = jstart; j < jstart + local_jmax; j++) {
                local_res += _calculate_partial_residual(
                    i, j, p, rhs, flag, rdx2, rdy2, imax, jmax,
                    pre_calculated_eps_Es, pre_calculated_eps_Ws,
                    pre_calculated_eps_Ns, pre_calculated_eps_Ss
                );
            }
        }

        MPI_Allreduce(&local_res, &global_res, 1, MPI_FLOAT, MPI_SUM, grid_comm);
        if ((sqrt(global_res / ifull) / p0) < eps) break;

        MPI_Win_fence(0, win);

        // update halo to other processes (north, south, west, east) with non-blocking 1-side communication
        // it actually just perform an exchange of halo rows and columns
        // for north and south, the south process sends the first row to the north process,
        // and the north process sends the last row to the south process
        // for west and east, the east process sends the first column to the west process
        // and the west process sends the last column to the east process
        //
        // FIXME: However, the current implementation is not correct, and we have no time to fix it.
        // The result shows an image which seems to be a wrong combination of two images.
        // Even when we skip those steps, the result is still incorrect.
        // 
        // Please see the report for more details.
        //
        if (north != MPI_PROC_NULL) {
            // MPI_Sendrecv_replace(&p[istart + local_imax][jstart], 1, halo_row, north, 0, south, 0, grid_comm, MPI_STATUS_IGNORE);
        }
        if (south != MPI_PROC_NULL) {
            // MPI_Sendrecv_replace(&p[istart][jstart], 1, halo_row, south, 0, north, 0, grid_comm, MPI_STATUS_IGNORE);
        }
        if (west != MPI_PROC_NULL) {
            // MPI_Sendrecv_replace(&p[istart][jstart + local_jmax], 1, halo_col, west, 0, east, 0, grid_comm, MPI_STATUS_IGNORE);
        }
        if (east != MPI_PROC_NULL) {
            // MPI_Sendrecv_replace(&p[istart][jstart], 1, halo_col, east, 0, west, 0, grid_comm, MPI_STATUS_IGNORE);
        }

        // MPI_Win_fence(0, win);
    }

    // MPI_Type_free(&halo_col);
    // MPI_Type_free(&halo_row);

    // MPI_Gather(&p[*ileft][0], imax / nprocs * (jmax + 2), MPI_FLOAT, &p[*ileft][0], imax / nprocs * (jmax + 2), MPI_FLOAT, tag, MPI_COMM_WORLD);
    // if(imax % nprocs) {
    //     if(proc == nprocs - 1) {
    //         MPI_Send(&p[imax - imax % nprocs + 1][0], imax % nprocs * (jmax + 2), MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
    //     }
    //     if(proc == 0) {
    //         MPI_Recv(&p[imax - imax % nprocs + 1][0], imax % nprocs * (jmax + 2), MPI_FLOAT, nprocs - 1, tag, MPI_COMM_WORLD, &status);
    //     }
    // }

    // MPI_Bcast(&p[0][0], (imax + 2) * (jmax + 2), MPI_FLOAT, tag, MPI_COMM_WORLD);
    return iter;
}


/* Update the velocity values based on the tentative
 * velocity values and the new pressure matrix
 */
void update_velocity(
    float** u,
    float** v,
    float** f,
    float** g,
    float** p,
    char** flag,
    int imax,
    int jmax,
    float del_t,
    float delx,
    float dely
) {
    #pragma omp parallel
    {
        #pragma omp for collapse(2) schedule(static) nowait
        for (register int i = 1; i <= imax - 1; i++) {
            for (register int j = 1; j <= jmax; j++) {
                /* only if both adjacent cells are fluid cells */
                if ((flag[i][j] & C_F) && (flag[i + 1][j] & C_F)) {
                    u[i][j] = f[i][j] - (p[i + 1][j] - p[i][j]) * del_t / delx;
                }
            }
        }

        #pragma omp for collapse(2) schedule(static) nowait
        for (register int i = 1; i <= imax; i++) {
            for (register int j = 1; j <= jmax - 1; j++) {
                /* only if both adjacent cells are fluid cells */
                if ((flag[i][j] & C_F) && (flag[i][j + 1] & C_F)) {
                    v[i][j] = g[i][j] - (p[i][j + 1] - p[i][j]) * del_t / dely;
                }
            }
        }
    }
}


/* Set the timestep size so that we satisfy the Courant-Friedrichs-Lewy
 * conditions (ie no particle moves more than one cell width in one
 * timestep). Otherwise, the simulation becomes unstable.
 */
void set_timestep_interval(
    float* del_t,
    int imax,
    int jmax,
    float delx,
    float dely,
    float** u,
    float** v,
    float Re,
    float tau
) {
    float umax, vmax, deltu, deltv, deltRe;

    /* del_t satisfying CFL conditions */
    if (tau >= 1.0e-10) { /* else no time stepsize control */
        umax = 1.0e-10;
        vmax = 1.0e-10;

        #pragma omp parallel
        {
            #pragma omp for collapse(2) reduction(max:umax) schedule(static)
            for (register int i = 0; i <= imax + 1; i++) {
                for (register int j = 1; j <= jmax + 1; j++) {
                    umax = max(fabs(u[i][j]), umax);
                }
            }

            #pragma omp for collapse(2) reduction(max:vmax) schedule(static)
            for (register int i = 1; i <= imax + 1; i++) {
                for (register int j = 0; j <= jmax + 1; j++) {
                    vmax = max(fabs(v[i][j]), vmax);
                }
            }
        }

        deltu = delx / umax;
        deltv = dely / vmax;
        deltRe = 1 / (1 / (delx * delx) + 1 / (dely * dely)) * Re / 2.0;

        if (deltu < deltv) {
            *del_t = min(deltu, deltRe);
        } else {
            *del_t = min(deltv, deltRe);
        }
        *del_t = tau * (*del_t); /* multiply by safety factor */
    }
}
