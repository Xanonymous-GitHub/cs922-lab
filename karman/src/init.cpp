#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <limits>
#include <memory>

#include "init.hpp"
#include "datadef.hpp"

void load_flag_from_pgm(
    const matrix<char>& flag,
    const int& imax,
    const int& jmax,
    const std::string_view& filename
) {
    std::ifstream fp(filename, std::ios::binary);
    if (!fp) {
        std::cout << "Couldn't open file '" << filename << "'\n";
        return;
    }

    // Read and validate the header
    std::string line, type;
    if (!std::getline(fp, type)) {
        std::cerr << "Error reading the type from '" << filename << "'.\n";
        return;
    }

    // Check if the file starts with the P5 magic number
    if (type.substr(0, 2) != "P5") {
        std::cerr << "'" << filename << "' is not a PGM file.\n";
        return;
    }

    // Skip comments
    while (fp.peek() == '#') {
        if (!std::getline(fp, line)) {
            std::cerr << "Error reading past comments in '" << filename << "'.\n";
            return;
        }
    }

    int width, height, maxVal;
    fp >> width >> height >> maxVal;

    // Basic validation
    if (!fp || width < 1 || height < 1 || maxVal < 1 || maxVal > 255) {
        std::cerr << "'" << filename << "' has invalid headers.\n";
        return;
    }

    // Important to ignore one or more whitespace characters after maxVal before reading pixel data
    fp.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    std::vector<char> pix(width * height);
    fp.read(pix.data(), static_cast<std::streamsize>(pix.size()));
    if (!fp) {
        std::cerr << "Error reading pixel data from '" << filename << "'.\n";
        return;
    }

    for (int j = 1; j < jmax + 2; j++) {
        if (j <= height) {
            fp.read(pix.data(), width);
        }

        for (int i = 1; i < imax + 2; i++) {
            if (j >= height + 1 || i >= width + 1) {
                flag->at(i, j) = C_F;
            } else {
                if (pix[i - 1] == 0) {
                    flag->at(i, j) = C_B;
                } else {
                    flag->at(i, j) = C_F;
                }
            }
        }
    }

    fp.close();
}

/* Initialize the flag array, marking any obstacle cells and the edge cells
 * as boundaries. The cells adjacent to boundary cells have their relevant
 * flags set too.
 */
void init_flag(
    const matrix<char>& flag,
    const int& imax,
    const int& jmax,
    const float& delx,
    const float& dely,
    const std::shared_ptr<int>& ibound
) {
    /* Mark a circular obstacle as boundary cells, the rest as fluid */
    const auto mx = 20.0 / 41.0 * jmax * dely;
    const auto my = mx;
    const auto rad1 = 5.0 / 41.0 * jmax * dely;

    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            const auto x = (i - 0.5) * delx - mx;
            const auto y = (j - 0.5) * dely - my;
            flag->at(i, j) = (x * x + y * y <= rad1 * rad1) ? C_B : C_F;
        }
    }

    /* Mark the north & south boundary cells */
    for (int i = 0; i <= imax + 1; i++) {
        flag->at(i, 0) = C_B;
        flag->at(i, jmax + 1) = C_B;
    }
    /* Mark the east and west boundary cells */
    for (int j = 1; j <= jmax; j++) {
        flag->at(0, j) = C_B;
        flag->at(imax + 1, j) = C_B;
    }

    /* flags for boundary cells */
    *(ibound) = 0;

    for (int i = 1; i <= imax; i++) {
        for (int j = 1; j <= jmax; j++) {
            if (!(flag->at(i, j) & C_F)) {
                (*ibound)++;
                if (flag->at(i - 1, j) & C_F) flag->at(i, j) |= B_W;
                if (flag->at(i + 1, j) & C_F) flag->at(i, j) |= B_E;
                if (flag->at(i, j - 1) & C_F) flag->at(i, j) |= B_S;
                if (flag->at(i, j + 1) & C_F) flag->at(i, j) |= B_N;
            }
        }
    }
}
