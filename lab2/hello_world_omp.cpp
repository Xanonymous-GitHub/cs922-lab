#include <omp.h>

#include <iostream>

int main(int argc, char const *argv[]) {
#pragma omp parallel
    {
        const int thread_id = omp_get_thread_num();
        std::cout << "Hello World from thread " << thread_id << std::endl;

#pragma omp barrier
        if (thread_id == 0) {
            std::cout << "There are " << omp_get_num_threads() << " threads" << std::endl;
        }
    }
    return 0;
}
