#include <iostream>
#include <memory>

int main(int argc, char const* argv[]) {
    const int my_int = 6;
    const auto my_int_ptr = std::make_unique<int>(my_int);

    std::cout << "The value of my_int is " << my_int << ", the value of pntr is " << my_int_ptr.get() << '\n';
    std::cout << "pntr is pointing at " << *my_int_ptr << '\n';

    const int a = 50;
    const auto a_ptr = std::make_unique<int>(a);

    std::cout << "pointer is " << a_ptr.get() << ", the value of " << *a_ptr << '\n';

    return 0;
}
