#include <cstdint>
#include <iostream>

uint_fast64_t factorial(const int& n) noexcept {
    uint_fast64_t result = 1;
    for (int i = 1; i <= n; i++) {
        result *= i;
    }
    return result;
}

int main(int argc, char const* argv[]) {
    const int factorial_input = 55;
    const auto factorial_output = factorial(factorial_input);
    std::cout << factorial_input << "! = " << factorial_output << '\n';
    return 0;
}
