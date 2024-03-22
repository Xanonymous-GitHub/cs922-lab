#include <cmath>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <string_view>
#include <memory>
#include "alloc.hpp"
#include "datadef.hpp"
#include "matrix.hpp"

#define PACKAGE "bin2ppm"
#define VERSION "1.0"
#define GETOPTS "hi:o:pv:Vz"

// Command line options
static const std::array<option, 8> long_opts = {
    option{"help", no_argument, nullptr, 'h'},
    option{"infile", required_argument, nullptr, 'i'},
    option{"outfile", required_argument, nullptr, 'o'},
    option{"plot-psi", no_argument, nullptr, 'p'},
    option{"plot-zeta", no_argument, nullptr, 'z'},
    option{"verbose", required_argument, nullptr, 'v'},
    option{"version", no_argument, nullptr, 'V'},
    option{nullptr, no_argument, nullptr, 0}
};

// Computation of stream function and vorticity
void calc_psi_zeta(
    const matrix<float>& u,
    const matrix<float>& v,
    const matrix<float>& psi,
    const matrix<float>& zeta,
    const matrix<char>& flag,
    const int& imax,
    const int& jmax,
    const float& delx,
    const float& dely
) {
    // Computation of the vorticity zeta in the upper right corner
    // of cell (i,j) (only if the corner is surrounded by fluid cells)
    for (int i = 1; i <= imax - 1; i++) {
        for (int j = 1; j <= jmax - 1; j++) {
            if (flag->at(i, j) & C_F && flag->at(i + 1, j) & C_F &&
                flag->at(i, j + 1) & C_F && flag->at(i + 1, j + 1) & C_F) {
                zeta->at(i, j) = (u->at(i, j + 1) - u->at(i, j)) / dely - (v->at(i + 1, j) - v->at(i, j)) / delx;
            } else {
                zeta->at(i, j) = 0.0;
            }
        }
    }

    // Computation of the stream function in the upper right corner
    // of cell (i,j) (only if bother lower cells are fluid cells)
    for (int i = 0; i <= imax; i++) {
        psi->at(i, 0) = 0.0;
        for (int j = 1; j <= jmax; j++) {
            psi->at(i, j) = psi->at(i, j - 1);
            if (flag->at(i, j) & C_F || flag->at(i + 1, j) & C_F) {
                psi->at(i, j) += u->at(i, j) * dely;
            }
        }
    }
}

static void print_usage(const std::string_view& progname) {
    std::cout << "Try '" << progname << " --help' for more information.\n";
}

static void print_version() {
    std::cout << PACKAGE << " " << VERSION << "\n";
}

static void print_help(const std::string_view& progname) {
    std::cerr
        << PACKAGE
        << ". Converts karman output into portable pixmaps.\n\n"
           "Usage: "
        << progname
        << " [OPTIONS]...\n\n"
           "  -h, --help            Print a summary of the options\n"
           "  -V, --version         Print the version number\n"
           "  -v, --verbose=LEVEL   Set the verbosity level. 0 is silent\n"
           "  -i, --infile=FILE     Read the simulation state from this file\n"
           "                        (defaults to standard input)\n"
           "  -o, --outfile=FILE    Write the image to this file\n"
           "                        (defaults to standard output)\n"
           "  -p, --plot-psi        Plot psi values in the image\n"
           "  -z, --plot-zeta       Plot zeta (vorticity) in the image\n";
}

int main(const int argc, char* const* argv) {
    float zmax = -1e10, zmin = 1e10;
    float pmax = -1e10, pmin = 1e10;
    float umax = -1e10, umin = 1e10;
    float vmax = -1e10, vmin = 1e10;

    int show_help = 0;
    int show_usage = 0;
    int show_version = 0;
    int verbose = 1;

    std::string input_file_name;
    std::string output_file_name;

    const auto progname = argv[0];

    while (true) {
        const auto option_char = getopt_long(
            argc,
            argv,
            GETOPTS,
            long_opts.cbegin(),
            nullptr
        );
        if (option_char == -1) {
            break;
        }

        switch (option_char) {
            case 'h':
                show_help = 1;
                break;
            case 'V':
                show_version = 1;
                break;
            case 'v':
                verbose = std::stoi(optarg);
                break;
            case 'i':
                input_file_name = std::string(optarg);
                break;
            case 'o':
                output_file_name = std::string(optarg);
                break;
            default:
                show_usage = 1;
        }
    }

    if (show_usage || optind < argc) {
        print_usage(progname);
        return 1;
    }

    if (show_version) {
        print_version();
        if (!show_help) {
            return 0;
        }
    }

    if (show_help) {
        print_help(progname);
        return 0;
    }

    std::ifstream fin;
    if (!input_file_name.empty()) {
        fin.open(input_file_name, std::ios::binary);
        if (!fin.is_open()) {
            std::cerr << "Could not open '" << input_file_name << "'\n";
            fin.close();
            return 1;
        }
    }

    std::ofstream fout;
    if (!output_file_name.empty()) {
        fout.open(output_file_name, std::ios::binary);
        if (!fout.is_open()) {
            std::cerr << "Could not open '" << output_file_name << "'\n";
            fout.close();
            return 1;
        }
    }

    int imax = 0, jmax = 0;
    fin.read(reinterpret_cast<char*>(&imax), sizeof(int));
    fin.read(reinterpret_cast<char*>(&jmax), sizeof(int));

    float x_length, y_length;
    fin.read(reinterpret_cast<char*>(&x_length), sizeof(float));
    fin.read(reinterpret_cast<char*>(&y_length), sizeof(float));

    const float delx = x_length / static_cast<float>(imax);
    const float dely = y_length / static_cast<float>(jmax);

    const auto& u = alloc_float_matrix(imax + 2, jmax + 2);
    const auto& v = alloc_float_matrix(imax + 2, jmax + 2);
    const auto& p = alloc_float_matrix(imax + 2, jmax + 2);
    const auto& psi = alloc_float_matrix(imax + 2, jmax + 2);
    const auto& zeta = alloc_float_matrix(imax + 2, jmax + 2);
    const auto& flag = alloc_char_matrix(imax + 2, jmax + 2);

    if (verbose > 1) {
        std::cout << "imax: " << imax << "\n";
        std::cout << "jmax: " << jmax << "\n";
        std::cout << "x_length: " << x_length << "\n";
        std::cout << "y_length: " << y_length << "\n";
    }

    const auto _float_jmax_size = static_cast<std::streamsize>((jmax + 2) * sizeof(float));
    const auto _char_jmax_size = static_cast<std::streamsize>((jmax + 2) * sizeof(char));
    for (int i = 0; i <= imax + 2; i++) {
        float _float_tmp;
        char _char_tmp;

        fin.read(reinterpret_cast<char*>(&_float_tmp), _float_jmax_size);
        u->at(i, 0) = _float_tmp;

        fin.read(reinterpret_cast<char*>(&_float_tmp), _float_jmax_size);
        v->at(i, 0) = _float_tmp;

        fin.read(reinterpret_cast<char*>(&_float_tmp), _float_jmax_size);
        p->at(i, 0) = _float_tmp;

        fin.read(reinterpret_cast<char*>(&_char_tmp), _char_jmax_size);
        flag->at(i, 0) = _char_tmp;
    }

    calc_psi_zeta(u, v, psi, zeta, flag, imax, jmax, delx, dely);
    fout << "P6 " << imax << " " << jmax << " 255\n";

    for (int j = 1; j < jmax + 1; j++) {
        for (int i = 1; i < imax + 1; i++) {
            int r, g, b;
            if (!(flag->at(i, j) & C_F)) {
                r = 0;
                b = 0;
                g = 255;
            } else {
                zmax = std::max(zmax, zeta->at(i, j));
                zmin = std::min(zmin, zeta->at(i, j));
                pmax = std::max(pmax, psi->at(i, j));
                pmin = std::min(pmin, psi->at(i, j));
                umax = std::max(umax, u->at(i, j));
                umin = std::min(umin, u->at(i, j));
                vmax = std::max(vmax, v->at(i, j));
                vmin = std::min(vmin, v->at(i, j));

                // when output mode == ZETA
                const auto _z = (i < imax && j < jmax) ? zeta->at(i, j) : 0.0;
                r = g = b = static_cast<int>(std::pow(std::fabs(_z / 12.6), .4) * 255);

                // when output mode == PSI
                // const auto _p = i < imax && j < jmax ? psi->at(i, j) : 0.0;
                // r = g = b = static_cast<int>((_p + 3.0) / 7.5 * 255);
            }
            fout << static_cast<char>(r) << static_cast<char>(g) << static_cast<char>(b);
        }
    }
    if (verbose > 0) {
        std::cout << "u:    " << umin << " -- " << umax << "\n";
        std::cout << "v:    " << vmin << " -- " << vmax << "\n";
        std::cout << "psi:  " << pmin << " -- " << pmax << "\n";
        std::cout << "zeta: " << zmin << " -- " << zmax << "\n";
    }

    fin.close();
    fout.close();

    return 0;
}
