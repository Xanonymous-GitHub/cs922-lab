#include <mpi.h>
#include <string_view>
#include <iostream>
#include <iomanip>
#include "alloc.hpp"

#define ROWS 9
#define COLS 7

/* A sample program that demonstrates how to time an MPI program, and how to
 * send (parts of) the rows and columns in a two-dimensional array between
 * two processes using MPI.
 * In the matrix representation used here, columns are simple to copy, as
 * column elements are contiguous in memory. Row copying is more complex as
 * consecutive row elements are not contiguous in memory and require the
 * use of derived MPI_Datatypes to perform the copy.
 */


void zero_matrix(const matrix<float>& matrix) {
    for (int i = 0; i < COLS; i++) {
        for (int j = 0; j < ROWS; j++) {
            matrix->at(i, j) = 0;
        }
    }
}

void print_matrix(const std::string_view& title, const matrix<float>& matrix) {
    std::cout << title << "\n";
    for (int i = 0; i < COLS; i++) {
        for (int j = 0; j < ROWS; j++) {
            std::cout
                << std::fixed
                << std::setprecision(2)
                << matrix->at(i, j)
                << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}

int main(int argc, char** argv) {
    int n = 0, p = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &n);
    MPI_Comm_rank(MPI_COMM_WORLD, &p);

    if (n != 2) {
        std::cout << "This program must be run with exactly 2 processes.\n";
        MPI_Finalize();
        return 1;
    }

    auto t = MPI_Wtime();

    /* Allocate an COLS * ROWS array. */
    const auto matrix = alloc_float_matrix(COLS, ROWS);

    /* Fill processor 1's matrix with numbers */
    for (int i = 0; i < COLS; i++) {
        for (int j = 0; j < ROWS; j++) {
            matrix->at(i, j) = static_cast<float>((i * 10) + j);
        }
    }

    /* Define two MPI_Datatypes for rows that we use later */
    MPI_Datatype full_row_type, part_row_type;
    MPI_Type_vector(ROWS, 1, ROWS, MPI_FLOAT, &full_row_type);
    MPI_Type_commit(&full_row_type);
    MPI_Type_vector(3, 1, ROWS, MPI_FLOAT, &part_row_type);
    MPI_Type_commit(&part_row_type);

    if (p == 0) {
        MPI_Status s;

        zero_matrix(matrix);
        MPI_Recv(&matrix->at(4, 0), ROWS, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, &s);
        print_matrix("After receiving column 4:", matrix);

        zero_matrix(matrix);
        MPI_Recv(&matrix->at(6, 2), 4, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, &s);
        print_matrix("After receiving column 6, rows 3-5:", matrix);

        zero_matrix(matrix);
        MPI_Recv(&matrix->at(3, 0), ROWS * 2, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, &s);
        print_matrix("After receiving column 3 and 4:", matrix);

        zero_matrix(matrix);
        MPI_Recv(&matrix->at(0, 0), ROWS * COLS, MPI_FLOAT, 1, 0, MPI_COMM_WORLD, &s);
        print_matrix("After receiving all columns:", matrix);

        zero_matrix(matrix);
        MPI_Recv(&matrix->at(0, 6), 1, full_row_type, 1, 0, MPI_COMM_WORLD, &s);
        print_matrix("After receiving row 6:", matrix);

        zero_matrix(matrix);
        MPI_Recv(&matrix->at(0, 1), 1, part_row_type, 1, 0, MPI_COMM_WORLD, &s);
        print_matrix("After receiving row 1 cols 0-2:", matrix);

        zero_matrix(matrix);
        MPI_Recv(&matrix->at(4, 1), 1, part_row_type, 1, 0, MPI_COMM_WORLD, &s);
        print_matrix("After receiving row 1 cols 4-6:", matrix);
    } else {
        /* Send all of column 4 to processor 0 */
        MPI_Send(&matrix->at(4, 0), ROWS, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

        /* Send column 6 rows 2-5 to processor 0 */
        MPI_Send(&matrix->at(6, 2), 4, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

        /* Send columns 3 and 4 to processor 0 */
        MPI_Send(&matrix->at(3, 0), ROWS * 2, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

        /* Send the entire matrix (ie all columns) to processor 0 */
        MPI_Send(&matrix->at(0, 0), ROWS * COLS, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);

        /* Send row 6 to processor 0 */
        MPI_Send(&matrix->at(0, 6), 1, full_row_type, 0, 0, MPI_COMM_WORLD);

        /* Send row 1 col 0-2 to processor 0 */
        MPI_Send(&matrix->at(0, 1), 1, part_row_type, 0, 0, MPI_COMM_WORLD);

        /* Send row 1 col 4-6 to processor 0 */
        MPI_Send(&matrix->at(4, 1), 1, part_row_type, 0, 0, MPI_COMM_WORLD);
    }
    if (p == 0) {
        t = MPI_Wtime() - t;
        std::cout << "Program took " << t << " secs to run.\n";
    }

    /* Free the derived MPI_Datatypes */
    MPI_Type_free(&full_row_type);
    MPI_Type_free(&part_row_type);

    MPI_Finalize();
    return 0;
}
