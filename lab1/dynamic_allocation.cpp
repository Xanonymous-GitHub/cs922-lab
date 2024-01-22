#include <array>
#include <iostream>
#include <string>

#define SIZE 5

int main(int argc, char const *argv[]) {
    std::array<int, SIZE> arr;
    int temp;
    int i = 0;
    
    // Read in the value
    // s until user enters nothing
    while (std::cin >> temp) {
        arr[i++] = temp;
    }

    i--;
    for (; i >= 0; i--) {
        std::cout << arr[i] << '\n';
    }

    return 0;
}
