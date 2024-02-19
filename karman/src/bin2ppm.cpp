#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <getopt.h>
#include <string>
#include "alloc.hpp"
#include "datadef.hpp"

#define max(x,y) (((x)>(y))?(x):(y))
#define min(x,y) (((x)<(y))?(x):(y))

void calc_psi_zeta(float** u, float** v, float** psi, float** zeta,
                   char** flag, int imax, int jmax, float delx, float dely);

static void print_usage(void);

static void print_version(void);

static void print_help(void);

static char* progname;

#define PACKAGE "bin2ppm"
#define VERSION "1.0"

/* Command line options */
static struct option long_opts[] = {
    {"help", 0, nullptr, 'h'},
    {"infile", 1, nullptr, 'i'},
    {"outfile", 1, nullptr, 'o'},
    {"plot-psi", 0, nullptr, 'p'},
    {"plot-zeta", 0, nullptr, 'z'},
    {"version", 0, nullptr, 'V'},
    {"verbose", 1, nullptr, 'v'},
    {0, 0, 0, 0}
};
#define GETOPTS "hi:o:pv:Vz"

/* Output modes */
#define ZETA 0
#define PSI  1

int main(int argc, char** argv) {
    int i, j, imax, jmax;
    int outmode = ZETA, verbose = 1;
    float xlength, ylength;
    float zmax = -1e10, zmin = 1e10;
    float pmax = -1e10, pmin = 1e10;
    float umax = -1e10, umin = 1e10;
    float vmax = -1e10, vmin = 1e10;

    int show_help = 0, show_usage = 0, show_version = 0;
    char *infile = nullptr, *outfile = nullptr;
    FILE *fin = stdin, *fout = stdout;

    progname = argv[0];
    int optc;
    while ((optc = getopt_long(argc, argv, GETOPTS, long_opts, nullptr)) != -1) {
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
            case 'i':
                if (infile != nullptr) {
                    free(infile);
                }
                infile = strdup(optarg);
                break;
            case 'o':
                if (outfile != nullptr) {
                    free(outfile);
                }
                outfile = strdup(optarg);
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

    if (infile != nullptr) {
        fin = fopen(infile, "rb");
        if (!fin) {
            fprintf(stderr, "Could not open '%s'\n", infile);
            return 1;
        }
    }

    if (outfile != nullptr) {
        fout = fopen(outfile, "wb");
        if (!fout) {
            fprintf(stderr, "Could not open '%s'\n", outfile);
            return 1;
        }
    }
    fread(&imax, sizeof(int), 1, fin);
    fread(&jmax, sizeof(int), 1, fin);

    float** u = alloc_floatmatrix(imax + 2, jmax + 2);
    float** v = alloc_floatmatrix(imax + 2, jmax + 2);
    float** p = alloc_floatmatrix(imax + 2, jmax + 2);
    float** psi = alloc_floatmatrix(imax + 2, jmax + 2);
    float** zeta = alloc_floatmatrix(imax + 2, jmax + 2);
    char** flag = alloc_charmatrix(imax + 2, jmax + 2);

    if (!u || !v || !p || !psi || !zeta || !flag) {
        fprintf(stderr, "Couldn't allocate memory for matrices.\n");
        return 1;
    }

    fread(&xlength, sizeof(float), 1, fin);
    fread(&ylength, sizeof(float), 1, fin);
    float delx = xlength / imax;
    float dely = ylength / jmax;

    if (verbose > 1) {
        printf("imax: %d\n", imax);
        printf("jmax: %d\n", jmax);
        printf("xlength: %g\n", xlength);
        printf("ylength: %g\n", ylength);
    }
    for (i = 0; i <= imax + 2; i++) {
        fread(u[i], sizeof(float), jmax + 2, fin);
        fread(v[i], sizeof(float), jmax + 2, fin);
        fread(p[i], sizeof(float), jmax + 2, fin);
        fread(flag[i], sizeof(char), jmax + 2, fin);
    }

    calc_psi_zeta(u, v, psi, zeta, flag, imax, jmax, delx, dely);
    fprintf(fout, "P6 %d %d 255\n", imax, jmax);

    for (j = 1; j < jmax + 1; j++) {
        for (i = 1; i < imax + 1; i++) {
            int r, g, b;
            if (!(flag[i][j] & C_F)) {
                r = 0;
                b = 0;
                g = 255;
            } else {
                zmax = max(zmax, zeta[i][j]);
                zmin = min(zmin, zeta[i][j]);
                pmax = max(pmax, psi[i][j]);
                pmin = min(pmin, psi[i][j]);
                umax = max(umax, u[i][j]);
                umin = min(umin, u[i][j]);
                vmax = max(vmax, v[i][j]);
                vmin = min(vmin, v[i][j]);
                if (outmode == ZETA) {
                    float z = (i < imax && j < jmax) ? zeta[i][j] : 0.0;
                    r = g = b = pow(fabs(z / 12.6), .4) * 255;
                } else if (outmode == PSI) {
                    float p = (i < imax && j < jmax) ? psi[i][j] : 0.0;
                    r = g = b = (p + 3.0) / 7.5 * 255;
                }
            }
            fprintf(fout, "%c%c%c", r, g, b);
        }
    }
    if (verbose > 0) {
        printf("u:    % .5e -- % .5e\n", umin, umax);
        printf("v:    % .5e -- % .5e\n", vmin, vmax);
        printf("psi:  % .5e -- % .5e\n", pmin, pmax);
        printf("zeta: % .5e -- % .5e\n", zmin, zmax);
    }
    fclose(fin);
    fclose(fout);

    free_matrix(u);
    free_matrix(v);
    free_matrix(p);
    free_matrix(psi);
    free_matrix(zeta);
    free_matrix(flag);

    return 0;
}

/* Computation of stream function and vorticity */
void calc_psi_zeta(float** u, float** v, float** psi, float** zeta,
                   char** flag, int imax, int jmax, float delx, float dely) {
    int i, j;

    /* Computation of the vorticity zeta at the upper right corner     */
    /* of cell (i,j) (only if the corner is surrounded by fluid cells) */
    for (i = 1; i <= imax - 1; i++) {
        for (j = 1; j <= jmax - 1; j++) {
            if ((flag[i][j] & C_F) && (flag[i + 1][j] & C_F) &&
                (flag[i][j + 1] & C_F) && (flag[i + 1][j + 1] & C_F)) {
                zeta[i][j] = (u[i][j + 1] - u[i][j]) / dely
                             - (v[i + 1][j] - v[i][j]) / delx;
            } else {
                zeta[i][j] = 0.0;
            }
        }
    }

    /* Computation of the stream function at the upper right corner    */
    /* of cell (i,j) (only if bother lower cells are fluid cells)      */
    for (i = 0; i <= imax; i++) {
        psi[i][0] = 0.0;
        for (j = 1; j <= jmax; j++) {
            psi[i][j] = psi[i][j - 1];
            if ((flag[i][j] & C_F) || (flag[i + 1][j] & C_F)) {
                psi[i][j] += u[i][j] * dely;
            }
        }
    }
}

static void print_usage(void) {
    fprintf(stderr, "Try '%s --help' for more information.\n", progname);
}

static void print_version(void) {
    fprintf(stderr, "%s %s\n", PACKAGE, VERSION);
}

static void print_help(void) {
    fprintf(stderr, "%s. Converts karman output into portable pixmaps.\n\n",
            PACKAGE);
    fprintf(stderr, "Usage: %s [OPTIONS]...\n\n", progname);
    fprintf(stderr, "  -h, --help            Print a summary of the options\n");
    fprintf(stderr, "  -V, --version         Print the version number\n");
    fprintf(stderr, "  -v, --verbose=LEVEL   Set the verbosity level. 0 is silent\n");
    fprintf(stderr, "  -i, --infile=FILE     Read the simulation state from this file\n");
    fprintf(stderr, "                        (defaults to standard input)\n");
    fprintf(stderr, "  -o, --outfile=FILE    Write the image to this file\n");
    fprintf(stderr, "                        (defaults to standard output)\n");
    fprintf(stderr, "  -p, --plot-psi        Plot psi values in the image\n");
    fprintf(stderr, "  -z, --plot-zeta       Plot zeta (vorticity) in the image\n");
}
