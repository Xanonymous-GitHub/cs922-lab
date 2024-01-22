#include <array>
#include <iostream>

constexpr float fahrenheit_to_celsius(const float &fahrenheit) noexcept {
    return (fahrenheit - 32) * 5 / 9;
}

int main(int argc, char const *argv[]) {
    // An array of temperatures in Fahrenheit, form 0 to 300.
    std::array<float, 301> fahrenheit_temperatures;
    for (int i = 0; i < fahrenheit_temperatures.size(); ++i) {
        fahrenheit_temperatures[i] = i;
    }

    // Convert the temperatures to Celsius.
    std::array<float, 301> celsius_temperatures;
    for (int i = 0; i < celsius_temperatures.size(); ++i) {
        celsius_temperatures[i] = fahrenheit_to_celsius(fahrenheit_temperatures[i]);
    }

    // Print the results.
    for (int i = 0; i < celsius_temperatures.size(); ++i) {
        std::cout << fahrenheit_temperatures[i] << "F = " << celsius_temperatures[i] << "C\n";
    }

    return 0;
}
