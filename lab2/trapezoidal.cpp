#include <cmath>
#include <format>
#include <iostream>

inline constexpr float f(const float& x) noexcept {
    return 2.0f / (1.0f + std::pow(x, 4));
}

int main(int argc, char const* argv[]) {
    const int a = 0;
    const int b = 1;
    const long n = 1024000000;
    const double h = (b - a) / (static_cast<const double>(n));

    double integral = (f(a) + f(b)) / 2.0f;
    double x = static_cast<const double>(a);

#pragma omp parallel for reduction(+ : integral)
    for (long i = 1; i <= n; ++i) {
        x += h;
        integral += f(x);
    }

    integral *= h;

    std::cout << std::format("With n = {} trapezoids, estimate: {}\n", n, integral);

    return 0;
}
