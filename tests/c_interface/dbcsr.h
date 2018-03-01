#ifndef DBCSR
#define DBCSR

#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#endif
    void c_dbcsr_init_lib();
    
    void c_dbcsr_finalize_lib_aux(MPI_Fint fcomm);
    
    void c_dbcsr_finalize_lib(MPI_Comm comm)
    {
        c_dbcsr_finalize_lib_aux(MPI_Comm_c2f(comm));
    }
  
    void c_dbcsr_distribution_new(void** dist, MPI_Fint fcomm, int* row_dist, int row_dist_size,
                                  int* col_dist, int col_dist_size);

    void c_dbcsr_distribution_release(void** dist);

    void test(char* str);

    void c_dbcsr_create_new_d(void** matrix, char* name, void* dist, char matrix_type, int* row_blk_sizes,
                              int row_blk_sizes_length, int* col_blk_sizes, int col_blk_sizes_length);

    void c_dbcsr_release(void** matrix);

#ifdef __cplusplus
}

namespace dbcsr {
  constexpr auto init_lib = c_dbcsr_init_lib;
  constexpr auto finalize_lib = c_dbcsr_finalize_lib;
}

#endif


#endif // DBCSR
