#include "Driver.hpp"
#include "InputFile.hpp"

#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include <omp.h>
#include <string>

int main(const int argc, const char* argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: deqn <filename>" << '\n';
        std::exit(1);
    }

    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    const std::string filename = argv[1];

    const InputFile input{filename};

    std::string problem_name(filename);

    if (const auto len = problem_name.length(); problem_name.substr(len - 3, 3) == ".in") {
        problem_name = problem_name.substr(0, len - 3);
    }

    // Strip out leading path
    size_t last_sep = problem_name.find_last_of('/');

    if (last_sep != std::string::npos) {
        last_sep = last_sep + 1;
    } else {
        last_sep = 0;
    }

    problem_name = problem_name.substr(last_sep, problem_name.size());

    MPI_Init(&const_cast<int &>(argc), const_cast<char ***>(&argv));
    const auto& has_initialized = std::make_unique<int>(0);
    MPI_Initialized(has_initialized.get());

    if (static_cast<bool>(*has_initialized)) {
        const Driver driver{input, problem_name};
        const auto& start_time = omp_get_wtime();
        driver.run();
        const auto& end_time = omp_get_wtime();
        MPI_Finalize();
        std::cout << "Elapsed time: " << end_time - start_time << " seconds" << '\n';
    } else {
        std::cerr << "Error: MPI has not been initialized" << '\n';
        std::exit(1);
    }

    return 0;
}
