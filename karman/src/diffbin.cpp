#include <cerrno>
#include <cmath>
#include <getopt.h>
#include <cstring>
#include <string>
#include <string_view>
#include <array>
#include <iostream>
#include <fstream>
#include <vector>

#define PACKAGE "diffbin"
#define VERSION "1.0"
#define GETOPTS "e:hm:V"
#define MODE_DIFF 0
#define MODE_OUTPUT_U 1
#define MODE_OUTPUT_V 2
#define MODE_OUTPUT_P 3
#define MODE_OUTPUT_FLAGS 4

/* Command line options */
static const std::array<option, 5> long_opts = {
    option{"help", no_argument, nullptr, 'h'},
    option{"version", no_argument, nullptr, 'V'},
    option{"epsilon", required_argument, nullptr, 'e'},
    option{"mode", required_argument, nullptr, 'm'},
    option{nullptr, 0, nullptr, 0}
};

static void print_usage(const std::string_view& progname) {
    std::cerr << "Try '" << progname << " --help' for more information.\n";
}

static void print_version() {
    std::cout << PACKAGE << " " << VERSION << '\n';
}

static void print_help(const std::string_view& progname) {
    std::cerr << PACKAGE << ". A utility to compare karman state files.\n\n";
    std::cerr << "Usage: " << progname << " [OPTIONS] FILE1 FILE2\n\n";
    std::cerr << "  -h, --help            Print a summary of the options\n";
    std::cerr << "  -V, --version         Print the version number\n";
    std::cerr << "  -e, --epsilon=EPSILON Set epsilon: the maximum allowed difference\n";
    std::cerr << "  -m, --mode=MODE       Set the mode, may be one of 'diff', 'plot-u',\n";
    std::cerr << "                        'plot-v', 'plot-p', or 'plot-flags'. The plot\n";
    std::cerr << "                        modes produce output ready to be used by the\n";
    std::cerr << "                        'splot matrix' command in gnuplot.\n";
}

int main(int argc, char** argv) {
    float epsilon = 1e-7;
    int mode = MODE_DIFF;
    int show_help = 0, show_usage = 0, show_version = 0;
    const std::string progname = argv[0];

    while (true) {
        const auto option_char = getopt_long(argc, argv, GETOPTS, long_opts.cbegin(), nullptr);
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
            case 'e':
                epsilon = std::stof(optarg);
                break;
            case 'm':
                if (strncasecmp(optarg, "diff", 4) == 0) {
                    mode = MODE_DIFF;
                } else if (strncasecmp(optarg, "plot-u", 6) == 0) {
                    mode = MODE_OUTPUT_U;
                } else if (strncasecmp(optarg, "plot-v", 6) == 0) {
                    mode = MODE_OUTPUT_V;
                } else if (strncasecmp(optarg, "plot-p", 6) == 0) {
                    mode = MODE_OUTPUT_P;
                } else if (strncasecmp(optarg, "plot-flags", 10) == 0) {
                    mode = MODE_OUTPUT_FLAGS;
                } else {
                    std::cerr << progname << ": Invalid mode '" << optarg << "'\n";
                    show_usage = 1;
                }
                break;
            default:
                show_usage = 1;
        }
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

    if (show_usage || optind != (argc - 2)) {
        print_usage(progname);
        return 1;
    }


    std::ifstream f1(argv[optind], std::ios::binary);
    if (!f1) {
        std::cerr << "Could not open '" << argv[optind] << "': " << std::strerror(errno) << "\n";
        return 1;
    }

    std::ifstream f2(argv[optind + 1], std::ios::binary);
    if (!f2) {
        std::cerr << "Could not open '" << argv[optind + 1] << "': " << std::strerror(errno) << "\n";
        return 1;
    }

    int imax, jmax, ii, jj;

    f1.read(reinterpret_cast<char*>(&imax), sizeof(int));
    f1.read(reinterpret_cast<char*>(&jmax), sizeof(int));
    f2.read(reinterpret_cast<char*>(&ii), sizeof(int));
    f2.read(reinterpret_cast<char*>(&jj), sizeof(int));

    if (ii != imax || jj != jmax) {
        std::cout << "Number of cells differ! (" << imax << "x" << jmax << " vs " << ii << "x" << jj << ")\n";
        return 1;
    }

    float x_length1, y_length1, x_length2, y_length2;

    f1.read(reinterpret_cast<char*>(&x_length1), sizeof(float));
    f1.read(reinterpret_cast<char*>(&y_length1), sizeof(float));
    f2.read(reinterpret_cast<char*>(&x_length2), sizeof(float));
    f2.read(reinterpret_cast<char*>(&y_length2), sizeof(float));

    if (x_length1 != x_length2 || y_length1 != y_length2) {
        std::cout << "Image domain dimensions differ! ("
                  << std::setprecision(6)
                  << x_length1 << "x"
                  << y_length1 << " vs "
                  << x_length2 << "x"
                  << y_length2 << ")\n";
        return 1;
    }

    const auto _jmax_plus_2 = jmax + 2;

    std::vector<float> u1(_jmax_plus_2), v1(_jmax_plus_2), p1(_jmax_plus_2);
    std::vector<char> flags1(_jmax_plus_2);

    std::vector<float> u2(_jmax_plus_2), v2(_jmax_plus_2), p2(_jmax_plus_2);
    std::vector<char> flags2(_jmax_plus_2);

    int diff_found = 0;
    const auto _float_jmax_size = static_cast<std::streamsize>(_jmax_plus_2 * sizeof(float));
    const auto _char_jmax_size = static_cast<std::streamsize>(_jmax_plus_2 * sizeof(char));
    for (int i = 0; i < imax + 2 && !diff_found; i++) {
        f1.read(reinterpret_cast<char*>(u1.data()), _float_jmax_size);
        f1.read(reinterpret_cast<char*>(v1.data()), _float_jmax_size);
        f1.read(reinterpret_cast<char*>(p1.data()), _float_jmax_size);
        f1.read(reinterpret_cast<char*>(flags1.data()), _char_jmax_size);

        f2.read(reinterpret_cast<char*>(u2.data()), _float_jmax_size);
        f2.read(reinterpret_cast<char*>(v2.data()), _float_jmax_size);
        f2.read(reinterpret_cast<char*>(p2.data()), _float_jmax_size);
        f2.read(reinterpret_cast<char*>(flags2.data()), _char_jmax_size);

        for (int j = 0; j < _jmax_plus_2 && !diff_found; j++) {
            const auto du = u1[j] - u2[j];
            const auto dv = v1[j] - v2[j];
            const auto dp = p1[j] - p2[j];
            const auto d_flags = flags1[j] - flags2[j];

            switch (mode) {
                case MODE_DIFF:
                    if (std::isnan(du) || std::isnan(dv) || std::isnan(dp) ||
                        std::isinf(du) || std::isinf(dv) || std::isinf(dp)) {
                        diff_found = 1;
                    }
                    if (std::fabs(du) > epsilon || std::fabs(dv) > epsilon ||
                        std::fabs(dp) > epsilon || std::fabs(d_flags) > epsilon) {
                        diff_found = 1;
                    }
                    break;
                case MODE_OUTPUT_U:
                    std::cout << du << ((j == jmax + 1) ? '\n' : ' ');
                    break;
                case MODE_OUTPUT_V:
                    std::cout << dv << ((j == jmax + 1) ? '\n' : ' ');
                    break;
                case MODE_OUTPUT_P:
                    std::cout << dp << ((j == jmax + 1) ? '\n' : ' ');
                    break;
                case MODE_OUTPUT_FLAGS:
                    std::cout << static_cast<int>(d_flags) << ((j == jmax + 1) ? '\n' : ' ');
                    break;
            }
        }
    }

    if (diff_found) {
        std::cout << "Files differ.\n";
        return 1;
    }

    if (mode == MODE_DIFF) {
        std::cout << "Files identical.\n";
    }

    return 0;
}
