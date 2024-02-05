#include <mpi.h>

int main(int argc, char* argv[]) {

    int size, rank, tag=0;
    int myarray[3];

    MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(rank == 0){
        myarray[0] = 1; 
        myarray[1] = 2; 
        myarray[2] = 3;
        MPI_Send(myarray, 3, MPI_INT, 1, tag, MPI_COMM_WORLD);
    } else {
        MPI_Recv(myarray, 3, MPI_INT, 0, tag, MPI_COMM_WORLD, 
                &status);
    }
    MPI_Finalize();

}
