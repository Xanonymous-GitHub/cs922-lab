#pragma once

#include <mpi.h>

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
);

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
);

int poisson(
    register float** p,
    register float** rhs,
    register char** flag,
    int imax,
    int jmax,
    int istart,
    int jstart,
    register float eps,
    register int itermax,
    register float omega,
    float* res,
    register int ifull,
    float pre_calculated_eps_Es[imax + 1][jmax + 1],
    float pre_calculated_eps_Ws[imax + 1][jmax + 1],
    float pre_calculated_eps_Ns[imax + 1][jmax + 1],
    float pre_calculated_eps_Ss[imax + 1][jmax + 1],
    register float rdx2,
    register float rdy2,
    register float beta_2,
    float pre_calculated_beta_mods[imax + 1][jmax + 1],
    MPI_Win win,
    MPI_Comm grid_comm
);

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
);

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
);
