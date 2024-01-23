#include <iostream>
#include <memory>
#include <utility>

void swap(int *a, int *b) {
    std::swap(*a, *b);
}

void show_value_and_address(const int &value) {
    std::cout << "Value: " << value << '\n';
    std::cout << "Address: " << &value << '\n';
}

int main(int argc, const char *argv[]) {
    int value_a = 50;
    int value_b = 5;

    /** Print out the values and memory locations of value_a and value_b here **/
    show_value_and_address(value_a);
    show_value_and_address(value_b);

    /** Use swap here **/
    swap(&value_a, &value_b);

    /** Print out the values and memory locations of value_a and value_b here **/
    show_value_and_address(value_a);
    show_value_and_address(value_b);
}