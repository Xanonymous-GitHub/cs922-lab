#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>

#include <mpi.h>
#include "alloc.h"
#include "boundary.h"
#include "datadef.h"
#include "init.h"
#include "simulation.h"

void write_bin(
    float** u,
    float** v,
    float** p,
    char** flag,
    int imax,
    int jmax,
    float xlength,
    float ylength,
    char* file
);

int read_bin(
    float** u,
    float** v,
    float** p,
    char** flag,
    int imax,
    int jmax,
    float xlength,
    float ylength,
    char* file
);

static void print_usage(void);

static void print_version(void);

static void print_help(void);

static char* progname;

int proc = 0;                       /* Rank of the current process */
int nprocs = 0;                /* Number of processes in communicator */

int* ileft, * iright;           /* Array bounds for each processor */

#define PACKAGE "karman"
#define VERSION "1.0"
#define SHOULD_RECORD_TIME 1

/* Command line options */
static struct option long_opts[] = {
    {"del-t",   1, NULL, 'd'},
    {"help",    0, NULL, 'h'},
    {"imax",    1, NULL, 'x'},
    {"infile",  1, NULL, 'i'},
    {"jmax",    1, NULL, 'y'},
    {"outfile", 1, NULL, 'o'},
    {"t-end",   1, NULL, 't'},
    {"verbose", 1, NULL, 'v'},
    {"version", 1, NULL, 'V'},
    {0,         0, 0,    0}
};
#define GETOPTS "d:hi:o:t:v:Vx:y:"

int main(int argc, char* argv[]) {
    int verbose = 1;          /* Verbosity level */
    const float xlength = 22.0;     /* Width of simulated domain */
    const float ylength = 4.1;      /* Height of simulated domain */
    register int imax = 660;           /* Number of cells horizontally */
    register int jmax = 120;           /* Number of cells vertically */

    char* infile;             /* Input raw initial conditions */
    char* outfile;            /* Output raw simulation results */

    float t_end = 2.1;        /* Simulation runtime */
    float del_t = 0.003;      /* Duration of each timestep */
    register const float tau = 0.5;          /* Safety factor for timestep control */

    register int itermax = 100;        /* Maximum number of iterations in SOR */
    register float eps = 0.001;        /* Stopping error threshold for SOR */
    register const float omega = 1.7;        /* Relaxation parameter for SOR */
    register const float gamma = 0.9;        /* Upwind differencing factor in PDE
                                 discretisation */

    register const float Re = 150.0;         /* Reynolds number */
    register const float ui = 1.0;           /* Initial X velocity */
    register const float vi = 0.0;           /* Initial Y velocity */

    register float t, delx, dely;
    register int i, j, itersor = 0, ifluid = 0;
    int ibound = 0;
    float res;
    register float** u, ** v, ** p, ** rhs, ** f, ** g;
    register char** flag;
    int init_case, iters = 0;
    int show_help = 0, show_usage = 0, show_version = 0;

    progname = argv[0];
    infile = strdup("karman.bin");
    outfile = strdup("karman.bin");

    int optc;
    while ((optc = getopt_long(argc, argv, GETOPTS, long_opts, NULL)) != -1) {
        switch (optc) {
            case 'h':
                show_help = 1;
                break;
            case 'V':
                show_version = 1;
                break;
            case 'v':
                verbose = atoi(optarg);
                break;
            case 'x':
                imax = atoi(optarg);
                break;
            case 'y':
                jmax = atoi(optarg);
                break;
            case 'i':
                free(infile);
                infile = strdup(optarg);
                break;
            case 'o':
                free(outfile);
                outfile = strdup(optarg);
                break;
            case 'd':
                del_t = atof(optarg);
                break;
            case 't':
                t_end = atof(optarg);
                break;
            default:
                show_usage = 1;
        }
    }
    if (show_usage || optind < argc) {
        print_usage();
        return 1;
    }

    if (show_version) {
        print_version();
        if (!show_help) {
            return 0;
        }
    }

    if (show_help) {
        print_help();
        return 0;
    }

    MPI_Init(&argc, &argv);

    delx = xlength / imax;
    dely = ylength / jmax;

    register const float rdx2 = 1.0 / (delx * delx);
    register const float rdy2 = 1.0 / (dely * dely);
    register const float beta_2 = -omega / (2.0 * (rdx2 + rdy2));

    /* Allocate arrays */
    u = alloc_floatmatrix(imax + 2, jmax + 2);
    v = alloc_floatmatrix(imax + 2, jmax + 2);
    f = alloc_floatmatrix(imax + 2, jmax + 2);
    g = alloc_floatmatrix(imax + 2, jmax + 2);
    p = alloc_floatmatrix(imax + 2, jmax + 2);
    rhs = alloc_floatmatrix(imax + 2, jmax + 2);
    flag = alloc_charmatrix(imax + 2, jmax + 2);

    if (!u || !v || !f || !g || !p || !rhs || !flag) {
        fprintf(stderr, "Couldn't allocate memory for matrices.\n");
        MPI_Finalize();
        return 1;
    }

    /* Read in initial values from a file if it exists */
    init_case = read_bin(u, v, p, flag, imax, jmax, xlength, ylength, infile);

    if (init_case > 0) {
        /* Error while reading file */
        MPI_Finalize();
        return 1;
    }

    #if SHOULD_RECORD_TIME
        double last_time_record = MPI_Wtime();
        double duration_of_1 = 0.0;
        double duration_of_2 = 0.0;
        double duration_of_3 = 0.0;
        double duration_of_4 = 0.0;
        double duration_of_5 = 0.0;
        double duration_of_6 = 0.0;
        double duration_of_7 = 0.0;
        double duration_of_8 = 0.0;
    #endif

    if (init_case < 0) {
        /* Set initial values if file doesn't exist */
        for (i = 0; i <= imax + 1; i++) {
            for (j = 0; j <= jmax + 1; j++) {
                u[i][j] = ui;
                v[i][j] = vi;
                p[i][j] = 0.0;
            }
        }

        // 1: init_flag
        #if SHOULD_RECORD_TIME
            if (verbose > 1) {
                last_time_record = MPI_Wtime();
            }
        #endif
        init_flag(flag, imax, jmax, delx, dely, &ibound);
        #if SHOULD_RECORD_TIME
            if (verbose > 1) {
                duration_of_1 += MPI_Wtime() - last_time_record;
            }
        #endif
        
        // 2: apply_boundary_conditions
        #if SHOULD_RECORD_TIME
            if (verbose > 1) {
                last_time_record = MPI_Wtime();
            }
        #endif
        apply_boundary_conditions(u, v, flag, imax, jmax, ui, vi);
        #if SHOULD_RECORD_TIME
            if (verbose > 1) {
                duration_of_2 += MPI_Wtime() - last_time_record;
            }
        #endif
    }

    ifluid = (imax * jmax) - ibound;

    float _pre_calculated_eps_Es[imax + 1][jmax + 1];
    float _pre_calculated_eps_Ws[imax + 1][jmax + 1];
    float _pre_calculated_eps_Ns[imax + 1][jmax + 1];
    float _pre_calculated_eps_Ss[imax + 1][jmax + 1];


    if (ifluid > 0) {
        #pragma omp parallel for private(i, j) collapse(2)
        for (i = 1; i <= imax; i++) {
            for (j = 1; j <= jmax; j++) {
                _pre_calculated_eps_Es[i][j] = flag[i + 1][j] & C_F ? 1 : 0;
                _pre_calculated_eps_Ws[i][j] = flag[i - 1][j] & C_F ? 1 : 0;
                _pre_calculated_eps_Ns[i][j] = flag[i][j + 1] & C_F ? 1 : 0;
                _pre_calculated_eps_Ss[i][j] = flag[i][j - 1] & C_F ? 1 : 0;
            }
        }
    }

    float _beta_mods[imax + 1][jmax + 1];

    #pragma omp parallel for private(i, j) collapse(2)
    for (i = 1; i <= imax; i++) {
        for (j = 1; j <= jmax; j++) {
            _beta_mods[i][j] = -omega / (
                (_pre_calculated_eps_Es[i][j] + _pre_calculated_eps_Ws[i][j]) * rdx2 + 
                (_pre_calculated_eps_Ns[i][j] + _pre_calculated_eps_Ss[i][j]) * rdy2
            );
        }
    }

    /* Main loop */
    for (t = 0.0; t < t_end; t += del_t, iters++) {
        // 3: set_timestep_interval
        #if SHOULD_RECORD_TIME
            if (verbose > 1) {
                last_time_record = MPI_Wtime();
            }
        #endif
        set_timestep_interval(&del_t, imax, jmax, delx, dely, u, v, Re, tau);
        #if SHOULD_RECORD_TIME
            if (verbose > 1) {
                duration_of_3 += MPI_Wtime() - last_time_record;
            }
        #endif

        // D: compute_tentative_velocity
        #if SHOULD_RECORD_TIME
            if (verbose > 1) {
                last_time_record = MPI_Wtime();
            }
        #endif
        compute_tentative_velocity(
            u, v, f, g, flag, 
            imax, jmax, del_t, delx, dely, gamma, Re
        );
        #if SHOULD_RECORD_TIME
            if (verbose > 1) {
                duration_of_4 += MPI_Wtime() - last_time_record;
            }
        #endif

        // E: compute_rhs
        #if SHOULD_RECORD_TIME
            if (verbose > 1) {
                last_time_record = MPI_Wtime();
            }
        #endif
        compute_rhs(f, g, rhs, flag, imax, jmax, del_t, delx, dely);
        #if SHOULD_RECORD_TIME
            if (verbose > 1) {
                duration_of_5 += MPI_Wtime() - last_time_record;
            }
        #endif

        if (ifluid > 0) {
            // F: poisson
            #if SHOULD_RECORD_TIME
                if (verbose > 1) {
                    last_time_record = MPI_Wtime();
                }
            #endif
            itersor = poisson(
                p, rhs, flag, imax, 
                jmax, eps,
                itermax, omega, &res, ifluid,
                _pre_calculated_eps_Es,
                _pre_calculated_eps_Ws,
                _pre_calculated_eps_Ns,
                _pre_calculated_eps_Ss,
                rdx2, rdy2, beta_2, _beta_mods
            );
            #if SHOULD_RECORD_TIME
                if (verbose > 1) {
                    duration_of_6 += MPI_Wtime() - last_time_record;
                }
            #endif
        } else {
            itersor = 0;
        }

        if (proc == 0 && verbose > 1) {
            printf("%d t:%g, del_t:%g, SOR iters:%3d, res:%e, bcells:%d\n",
                   iters, t + del_t, del_t, itersor, res, ibound);
        }

        // G: update_velocity
        #if SHOULD_RECORD_TIME
            if (verbose > 1) {
                last_time_record = MPI_Wtime();
            }
        #endif
        update_velocity(u, v, f, g, p, flag, imax, jmax, del_t, delx, dely);
        #if SHOULD_RECORD_TIME
            if (verbose > 1) {
                duration_of_7 += MPI_Wtime() - last_time_record;
            }
        #endif

        // H: apply_boundary_conditions
        #if SHOULD_RECORD_TIME
            if (verbose > 1) {
                last_time_record = MPI_Wtime();
            }
        #endif
        apply_boundary_conditions(u, v, flag, imax, jmax, ui, vi);
        #if SHOULD_RECORD_TIME
            if (verbose > 1) {
                duration_of_8 += MPI_Wtime() - last_time_record;
            }
        #endif
    } /* End of main loop */

    if (outfile != NULL && strcmp(outfile, "") != 0 && proc == 0) {
        write_bin(u, v, p, flag, imax, jmax, xlength, ylength, outfile);
    }

    free_matrix(u);
    free_matrix(v);
    free_matrix(f);
    free_matrix(g);
    free_matrix(p);
    free_matrix(rhs);
    free_matrix(flag);

    #if SHOULD_RECORD_TIME
        if (verbose > 1 && proc == 0) {
            printf("\n>>> init_flag [1]: %f s\n\n", duration_of_1);
            printf("\n>>> apply_boundary_conditions [2]: %f s\n\n", duration_of_2);
            printf("\n>>> set_timestep_interval [3]: %f s\n\n", duration_of_3);
            printf("\n>>> compute_tentative_velocity [4]: %f s\n\n", duration_of_4);
            printf("\n>>> compute_rhs [5]: %f s\n\n", duration_of_5);
            printf("\n>>> poisson [6]: %f s\n\n", duration_of_6);
            printf("\n>>> update_velocity [7]: %f s\n\n", duration_of_7);
            printf("\n>>> apply_boundary_conditions [8]: %f s\n\n", duration_of_8);
        }
    #endif

    MPI_Finalize();
    return 0;
}

/* Save the simulation state to a file */
void write_bin(
    float** u,
    float** v,
    float** p,
    char** flag,
    int imax,
    int jmax,
    float xlength,
    float ylength,
    char* file
) {
    int i;
    FILE* fp;

    fp = fopen(file, "wb");

    if (fp == NULL) {
        fprintf(stderr, "Could not open file '%s': %s\n", file,
                strerror(errno));
        return;
    }

    fwrite(&imax, sizeof(int), 1, fp);
    fwrite(&jmax, sizeof(int), 1, fp);
    fwrite(&xlength, sizeof(float), 1, fp);
    fwrite(&ylength, sizeof(float), 1, fp);

    for (i = 0; i < imax + 2; i++) {
        fwrite(u[i], sizeof(float), jmax + 2, fp);
        fwrite(v[i], sizeof(float), jmax + 2, fp);
        fwrite(p[i], sizeof(float), jmax + 2, fp);
        fwrite(flag[i], sizeof(char), jmax + 2, fp);
    }
    fclose(fp);
}

/* Read the simulation state from a file */
int read_bin(
    float** u,
    float** v,
    float** p,
    char** flag,
    int imax,
    int jmax,
    float xlength,
    float ylength,
    char* file
) {
    int i, j;
    FILE* fp;

    if (file == NULL) return -1;

    if ((fp = fopen(file, "rb")) == NULL) {
        fprintf(stderr, "Could not open file '%s': %s\n", file,
                strerror(errno));
        fprintf(stderr, "Generating default state instead.\n");
        return -1;
    }

    fread(&i, sizeof(int), 1, fp);
    fread(&j, sizeof(int), 1, fp);
    float xl, yl;
    fread(&xl, sizeof(float), 1, fp);
    fread(&yl, sizeof(float), 1, fp);

    if (i != imax || j != jmax) {
        fprintf(stderr, "Warning: imax/jmax have wrong values in %s\n", file);
        fprintf(stderr, "%s's imax = %d, jmax = %d\n", file, i, j);
        fprintf(stderr, "Program's imax = %d, jmax = %d\n", imax, jmax);
        return 1;
    }
    if (xl != xlength || yl != ylength) {
        fprintf(stderr, "Warning: xlength/ylength have wrong values in %s\n", file);
        fprintf(stderr, "%s's xlength = %g,  ylength = %g\n", file, xl, yl);
        fprintf(stderr, "Program's xlength = %g, ylength = %g\n", xlength,
                ylength);
        return 1;
    }

    for (i = 0; i < imax + 2; i++) {
        fread(u[i], sizeof(float), jmax + 2, fp);
        fread(v[i], sizeof(float), jmax + 2, fp);
        fread(p[i], sizeof(float), jmax + 2, fp);
        fread(flag[i], sizeof(char), jmax + 2, fp);
    }
    fclose(fp);
    return 0;
}

static void print_usage(void) {
    fprintf(stderr, "Try '%s --help' for more information.\n", progname);
}

static void print_version(void) {
    fprintf(stderr, "%s %s\n", PACKAGE, VERSION);
}

static void print_help(void) {
    fprintf(stderr, "%s. A simple computational fluid dynamics tutorial.\n\n", PACKAGE);
    fprintf(stderr, "Usage: %s [OPTIONS]...\n\n", progname);
    fprintf(stderr, "  -h, --help            Print a summary of the options\n");
    fprintf(stderr, "  -V, --version         Print the version number\n");
    fprintf(stderr, "  -v, --verbose=LEVEL   Set the verbosity level. 0 is silent\n");
    fprintf(stderr, "  -x, --imax=IMAX       Set the number of interior cells in the X direction\n");
    fprintf(stderr, "  -y, --jmax=JMAX       Set the number of interior cells in the Y direction\n");
    fprintf(stderr, "  -t, --t-end=TEND      Set the simulation end time\n");
    fprintf(stderr, "  -d, --del-t=DELT      Set the simulation timestep size\n");
    fprintf(stderr, "  -i, --infile=FILE     Read the initial simulation state from this file\n");
    fprintf(stderr, "                        (default is 'karman.bin')\n");
    fprintf(stderr, "  -o, --outfile=FILE    Write the final simulation state to this file\n");
    fprintf(stderr, "                        (default is 'karman.bin')\n");
}
