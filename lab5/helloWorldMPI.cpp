#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
    char idstr[32];
    int numprocs;
    int myid;
    int i;
    MPI_Status stat; 

    /* all MPI programs start with MPI_Init */
       MPI_Init(&argc,&argv);     

    /* Comm_size tells us how many processes there are  */
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs); 

    /* Comm_rank finds the rank of the process */
    MPI_Comm_rank(MPI_COMM_WORLD,&myid); 

    /* Print out a message */
    printf("Hello world, from process %d of %d\n", myid, numprocs);  

    /* MPI Programs end with MPI Finalize; this is a weak synchronization point */
    MPI_Finalize(); 

    return 0; 
} 
