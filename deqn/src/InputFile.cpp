#include "InputFile.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

InputFile::InputFile(const std::string& filename) {
    std::ifstream ifs{filename};

    if (!ifs.good()) {
        std::cerr << "File " << filename << " not found!" << std::endl;
        std::exit(1);
    }

    while (true) {
        std::string line;

        std::getline(ifs, line);

        if (ifs.eof()) {
            break;
        }

        std::istringstream iss{line};
        std::string key;

        iss >> key;
        if (key.empty() || key[0] == '#') {
            continue;
        }

        // TODO: migrate to `contains` when using c++20 or above.
        if (const auto& found_iter = pairs.find(key); found_iter != pairs.cend()) {
            std::cerr << "Duplicate key " << key << " in input file" << std::endl;
            std::exit(1);
        }

        std::string val;
        std::getline(iss, val);
        pairs[key] = val;
    }

    ifs.close();
}

template<typename T>
T InputFile::get(const std::string& name, const T& dfault) const {
    const auto itr = pairs.find(name);

    // TODO: migrate to `contains` when using c++20 or above.
    if (itr == pairs.cend()) {
        return dfault;
    }

    std::istringstream iss{itr->second};

    T val;
    iss >> val;

    return val;
}

int InputFile::getInt(const std::string& name, const int& dfault) const {
    return get(name, dfault);
}

double InputFile::getDouble(const std::string& name, const double& dfault) const {
    return get(name, dfault);
}

std::string InputFile::getString(const std::string& name, const std::string& dfault) const {
    return get(name, dfault);
}

std::vector<double> InputFile::getDoubleList(
    const std::string& name,
    const std::vector<double>& dfault
) const {
    const auto itr = pairs.find(name);

    // TODO: migrate to `contains` when using c++20 or above.
    if (itr == pairs.cend()) {
        return dfault;
    }

    std::istringstream iss{itr->second};

    std::vector<double> vallist;
    double val;

    while (iss >> val) {
        vallist.push_back(val);
    }

    return vallist;
}
