#include "matrix.hpp"

// Allocate memory for a rows*cols tuple of floats.
// The elements within a column are contiguous in memory, and columns
// themselves are also contiguous in memory.
matrix<float> alloc_float_matrix(const int& cols, const int& rows) {
    return make_matrix<float>(cols, rows);
}

// Allocate memory for a rows*cols array of chars.
matrix<char> alloc_char_matrix(const int& cols, const int& rows) {
    return make_matrix<char>(cols, rows);
}
