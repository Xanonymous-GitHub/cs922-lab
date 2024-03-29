#include "InputFile.hpp"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

InputFile::InputFile(const std::string& filename) {
    std::ifstream ifs{filename};

    if (!ifs.good()) {
        std::cerr << "File " << filename << " not found!" << '\n';
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

        if (pairs.contains(key)) {
            std::cerr << "Duplicate key " << key << " in input file" << '\n';
        }

        std::string val;
        std::getline(iss, val);
        pairs[key] = val;
    }

    ifs.close();
}

template<typename T>
T InputFile::get(const std::string& name, const T& dfault) const {
    if (!pairs.contains(name)) {
        return dfault;
    }

    const auto itr = pairs.find(name);
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
    if (!pairs.contains(name)) {
        return dfault;
    }

    const auto itr = pairs.find(name);
    std::istringstream iss{itr->second};

    std::vector<double> vallist;
    double val;

    while (iss >> val) {
        vallist.push_back(val);
    }

    return vallist;
}
