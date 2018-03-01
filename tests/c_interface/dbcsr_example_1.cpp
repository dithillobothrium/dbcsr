#include <mpi.h>

#include <cstdint>
#include <vector>
#include <iostream>

#include <stdio.h>

#include "dbcsr.h"

using namespace std;


vector<int> random_dist(int dist_size, int nbins)
{
    vector<int> dist(dist_size);

    for(int i=0; i < dist_size; i++)
    {
        dist[i] = (nbins-i+1) % nbins;
    }

    return std::move(dist);
}


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

    std::cout << "I'm processor " 
         << mpi_rank 
         << " over "
         << mpi_size 
         << " proc, ("
         << coord[0] << ", " << coord[1] 
         << ") in the 2D grid" << std::endl;

    test((char*)"test pass string");

    c_dbcsr_init_lib();


   // dbcsr::init_lib();

    int nblkrows_total = 4;
    int nblkcols_total = 4;

    vector<int> row_blk_sizes(nblkrows_total, 2), col_blk_sizes(nblkcols_total, 2); 
   
    auto row_dist = random_dist(nblkrows_total, dims[0]);
    auto col_dist = random_dist(nblkcols_total, dims[1]);
  
    for (auto a: row_dist) cout<<a<<"\n";
   // dbcsr::finalize_lib(group);

    void* dist = nullptr;

    c_dbcsr_distribution_new(&dist, MPI_Comm_c2f(group), row_dist.data(), row_dist.size(), 
                             col_dist.data(), col_dist.size());
    
    void* matrix = nullptr;

    c_dbcsr_create_new_d(&matrix, (char*)"fish chips", dist, 'N', row_blk_sizes.data(), row_blk_sizes.size(),
                         col_blk_sizes.data(), col_blk_sizes.size());

    c_dbcsr_release(&matrix);
    c_dbcsr_distribution_release(&dist);
    c_dbcsr_finalize_lib(group);
    
    MPI_Comm_free(&group);
    MPI_Finalize();
  
  return 0;
}
