#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
    int i, j;
    float du2dx, duvdy, duvdx, dv2dy, laplu, laplv;

    for (i = 1; i <= imax - 1; i++) {
        for (j = 1; j <= jmax; j++) {
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

    for (i = 1; i <= imax; i++) {
        for (j = 1; j <= jmax - 1; j++) {
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

    /* f & g at external boundaries */
    for (j = 1; j <= jmax; j++) {
        f[0][j] = u[0][j];
        f[imax][j] = u[imax][j];
    }
    for (i = 1; i <= imax; i++) {
        g[i][0] = v[i][0];
        g[i][jmax] = v[i][jmax];
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
    int i, j;

    for (i = 1; i <= imax; i++) {
        for (j = 1; j <= jmax; j++) {
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
    int i,
    int j,
    float** p,
    float** rhs,
    char** flag,
    int imax,
    int jmax,
    float omega,
    float rdx2,
    float rdy2,
    float beta_2,
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
        const float _eps_E = pre_calculated_eps_Es[i][j];
        const float _eps_W = pre_calculated_eps_Ws[i][j];
        const float _eps_N = pre_calculated_eps_Ns[i][j];
        const float _eps_S = pre_calculated_eps_Ss[i][j];

        p[i][j] = (1. - omega) * p[i][j] -
                    pre_calculated_beta_mods[i][j] * (
                        (_eps_E * p[i + 1][j] + _eps_W * p[i - 1][j]) * rdx2
                        + (_eps_N * p[i][j + 1] + _eps_S * p[i][j - 1]) * rdy2
                        - rhs[i][j]
                    );
    }
}


/* Red/Black SOR to solve the poisson equation */
int poisson(
    float** p,
    float** rhs,
    char** flag,
    int imax,
    int jmax,
    float delx,
    float dely,
    float eps,
    int itermax,
    float omega,
    float* res,
    int ifull,
    float pre_calculated_eps_Es[imax + 1][jmax + 1],
    float pre_calculated_eps_Ws[imax + 1][jmax + 1],
    float pre_calculated_eps_Ns[imax + 1][jmax + 1],
    float pre_calculated_eps_Ss[imax + 1][jmax + 1],
    float rdx2,
    float rdy2,
    float beta_2,
    float pre_calculated_beta_mods[imax + 1][jmax + 1]
) {
    int i = 1, j = 1, iter = 0;
    float p0 = 0.0;

    /* Calculate sum of squares */
    for (i = 1; i <= imax; i++) {
        for (j = 1; j <= jmax; j++) {
            if (flag[i][j] & C_F) { p0 += p[i][j] * p[i][j]; }
        }
    }

    p0 = sqrt(p0 / ifull);
    if (p0 < 0.0001) { p0 = 1.0; }

    /* Red/Black SOR-iteration */
    for (iter = 0; iter < itermax; iter++) {

        for (i = 1; i <= imax; i++) {
            for (j = 2 - (i % 2); j <= jmax; j += 2) {
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

        for (i = 1; i <= imax; i += 1) {
            for (j = 1 + (i % 2); j <= jmax; j += 2) {
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

        /* Partial computation of residual */
        *res = 0.0;
        for (i = 1; i <= imax; i++) {
            for (j = 1; j <= jmax; j++) {
                if (flag[i][j] & C_F) {
                    /* only fluid cells */
                    const float _eps_E = pre_calculated_eps_Es[i][j];
                    const float _eps_W = pre_calculated_eps_Ws[i][j];
                    const float _eps_N = pre_calculated_eps_Ns[i][j];
                    const float _eps_S = pre_calculated_eps_Ss[i][j];

                    const float add = (_eps_E * (p[i + 1][j] - p[i][j]) -
                           _eps_W * (p[i][j] - p[i - 1][j])) * rdx2 +
                          (_eps_N * (p[i][j + 1] - p[i][j]) -
                           _eps_S * (p[i][j] - p[i][j - 1])) * rdy2 - rhs[i][j];
                    *res += add * add;
                }
            }
        }
        *res = sqrt((*res) / ifull) / p0;

        /* convergence? */
        if (*res < eps) break;
    } /* end of iter */

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
    int i, j;

    for (i = 1; i <= imax - 1; i++) {
        for (j = 1; j <= jmax; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i + 1][j] & C_F)) {
                u[i][j] = f[i][j] - (p[i + 1][j] - p[i][j]) * del_t / delx;
            }
        }
    }
    for (i = 1; i <= imax; i++) {
        for (j = 1; j <= jmax - 1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i][j + 1] & C_F)) {
                v[i][j] = g[i][j] - (p[i][j + 1] - p[i][j]) * del_t / dely;
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
    int i, j;
    float umax, vmax, deltu, deltv, deltRe;

    /* del_t satisfying CFL conditions */
    if (tau >= 1.0e-10) { /* else no time stepsize control */
        umax = 1.0e-10;
        vmax = 1.0e-10;
        for (i = 0; i <= imax + 1; i++) {
            for (j = 1; j <= jmax + 1; j++) {
                umax = max(fabs(u[i][j]), umax);
            }
        }
        for (i = 1; i <= imax + 1; i++) {
            for (j = 0; j <= jmax + 1; j++) {
                vmax = max(fabs(v[i][j]), vmax);
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
