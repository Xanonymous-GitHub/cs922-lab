#include <stdio.h>
#include <mpi.h>

int main (int argc, char** argv) {
  int my_rank;     // My rank
  int p;           // Total number of processes
  int a = 0.0;     // Left endpoint
  int b = 1.0;     // Right endpoint
  int n = 1024;    // Number of trapezoids
  float h;         // Base length
  float local_a;   // My left endpoint
  float local_b;   // My right endpoint
  float local_n;    // My number of trapezoids
  float integral;  // My integral
  float total;     // Total integral
  int source;      // Sending process id
  int dest = 0;    // Send partial totals to 0
  int tag = 0;     // Message tag

  // Status of the received message
  MPI_Status status;

  // This function will calculate the local integral
  float trap(float local_a, float local_b, int local_n,
      float h);

  MPI_Init(&argc, &argv);

  // Find out my rank
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  // How many processes are there?
  MPI_Comm_size(MPI_COMM_WORLD, &p);

  h = (b - a)/(float)n; // width of trapezoid

  // p processes, so each process has n/p trapezoids
  local_n = n/p;

  // Work out local bounds, based on rank
  local_a = a + my_rank*local_n*h;
  local_b = local_a + local_n*h;

  // Find partial integral estimate
  integral = trap(local_a, local_b, local_n, h);

  /* BEGIN SUM */
  // The root process (0) will receive the results from all
  // the other processess
  if(my_rank ==  0) {
    total = integral;
    for (source = 1; source < p; source++) {
      MPI_Recv(&integral, 1, MPI_FLOAT, source, tag, MPI_COMM_WORLD,
          &status);

      total = total + integral;
    }
  } else {
    // Other processes just send their result to root
    MPI_Send(&integral, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
  }
  /* END SUM */

  if (my_rank == 0)
    printf("With n = %d trapezoids, estimate is %f\n", n, total);

  MPI_Finalize();
}

// The trap function just implements the integral estimation, like
// the serial code
float trap(float local_a, float local_b, int local_n, float h) {
  float integral;
  float x;
  int i;

  float f(float x);

  integral = (f(local_a) + f(local_b))/2.0;
  x = local_a;

  for(i = 1; i <= local_n; i++) {
    x = x + h;
    integral = integral + f(x);
  }

  integral = integral*h;

  return integral;
}

// The function to integrate, here we use x^2
float f(float x) {
  return x*x;
}
