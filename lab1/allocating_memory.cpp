#include <array>
#include <cstdlib>
#include <iostream>

#define SIZE 5

int main(int argc, char const *argv[]) {
    const auto my_array = std::array<int, 5>{10, 11, 12, 14, 15};
    std::cout << "The first element is: " << my_array.front() << std::endl;

    const int *my_malloc_array;
    my_malloc_array = reinterpret_cast<const int *>(std::malloc(SIZE * sizeof(const int)));

    if (my_malloc_array != nullptr) {
        for (int i = 0; i < SIZE; i++) {
            std::construct_at(my_malloc_array + i, i);
        }

        for (int i = 0; i < SIZE; i++) {
            std::cout << my_malloc_array[i] << std::endl;
        }
    }

    std::free((void *)my_malloc_array);

    return 0;
}
