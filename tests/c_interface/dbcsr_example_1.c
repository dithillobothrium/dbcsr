#include "dbcsr.h"
#include "mpi.h"

#include <stdio.h>

int main(int argc, char** argv)
{
  MPI_Init(NULL, NULL);

  int mpi_size;
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  // Make 2D grid
  int dims[2] = {0};
  MPI_Dims_create(mpi_size, 2, dims);
  const int periods[2] = {1};
  int reorder = 0;
  MPI_Comm group;
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &group);
  
  int coord[2];
  MPI_Cart_coords(group, mpi_rank, 2, coord);

  printf("I'm processor %d over %d proc, (%d, %d) in the 2D grid \n",
	 mpi_rank, mpi_size, coord[0], coord[1]);

  c_dbcsr_init_lib();
  c_dbcsr_finalize_lib(group);
  
  MPI_Comm_free(&group);
  MPI_Finalize();
  
  return 0;
}
