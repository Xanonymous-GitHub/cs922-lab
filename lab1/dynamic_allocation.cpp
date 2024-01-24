#include <array>
#include <iostream>
#include <string>

#define SIZE 5

int main(int argc, char const* argv[]) {
    std::array<int, SIZE> arr = {};

    std::string input = "";
    int i = 0;

    // Read in the value
    // s until user enters nothing
    while (i < SIZE && std::getline(std::cin, input) && !input.empty()) {
        arr[i++] = std::stoi(input);
    }

    auto it = arr.rbegin();
    it += SIZE - i;
    for (; it != arr.rend(); ++it) {
        std::cout << *it << std::endl;
    }

    return 0;
}
