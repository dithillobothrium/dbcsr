

module carma

   USE dbcsr_dist_methods,              ONLY: &
        dbcsr_distribution_col_dist, dbcsr_distribution_get_num_images_1d, &
        dbcsr_distribution_has_threads, dbcsr_distribution_hold, dbcsr_distribution_make_threads, &
        dbcsr_distribution_mp, dbcsr_distribution_ncols, dbcsr_distribution_no_threads, &
        dbcsr_distribution_nrows,  dbcsr_distribution_row_dist, dbcsr_distribution_new, &
        dbcsr_distribution_release


   USE dbcsr_methods,                   ONLY: &
        dbcsr_col_block_offsets, dbcsr_col_block_sizes, dbcsr_destroy_array, dbcsr_distribution, &
        dbcsr_get_matrix_type, dbcsr_has_symmetry, dbcsr_image_dist_release, dbcsr_nblkcols_total, &
        dbcsr_nfullcols_total, dbcsr_nfullrows_total,  dbcsr_release_locals, &
        dbcsr_row_block_offsets, dbcsr_release, dbcsr_mp_release

     USE dbcsr_types,                     ONLY: &
              dbcsr_type, dbcsr_data_obj, dbcsr_distribution_obj, dbcsr_mp_obj,dbcsr_type_no_symmetry

   USE dbcsr_mp_methods,                ONLY: dbcsr_mp_group,&
                                              dbcsr_mp_npcols,&
                                              dbcsr_mp_nprows,&
                                              dbcsr_mp_numnodes,&
                                              dbcsr_mp_new

  USE message_passing,                 ONLY: mp_cart_create,&
                                             mp_comm_free,&
                                             mp_environ,&
                                             mp_world_finalize,&
                                             mp_world_init,&
                                             mp_sync,&
                                             mp_cart_rank,&
                                             mp_max, &
                                             mp_cart_sub
                                             
    use dbcsr_work_operations

    use dbcsr_io, only: dbcsr_print

    use dbcsr_index_operations
    use dbcsr_multiply_api
contains
    
  !
  ! CARMA-like multiplication 
  !  --------   -----
  !  |      | x |   |
  !  --------   |   |
  !             |   |
  !             -----
  !
  subroutine carma_mm(matrix_a, matrix_b, matrix_c)
    type(dbcsr_type), intent(inout) :: matrix_a, matrix_b, matrix_c  
    logical, dimension(2)           :: free_coords
    integer, dimension(1,1)         :: loc_grid
    integer                         :: group, row_group
    
    integer                         :: dist_nrows_a, dist_ncols_a, &
                                       dist_nrows_b, dist_ncols_b, &
                                       dist_nrows_c, dist_ncols_c

    integer, dimension(:), pointer  :: local_row_dist_a, local_col_dist_a, &
                                       local_row_dist_b, local_col_dist_b, &
                                       local_row_dist_c, local_col_dist_c, &
                                       row_blk_size_c, col_blk_size_c
    
    integer, dimension(:), pointer  :: row_dist, col_dist
    type(dbcsr_distribution_obj)    :: dist, loc_dist_a, loc_dist_b, loc_dist_c
    integer                         :: ierr
    type(dbcsr_type)                :: loc_matr_a, loc_matr_b, loc_matr_c
    type(dbcsr_mp_obj)              :: mp_env

    type(dbcsr_data_obj)            :: data_buf_a, data_buf_b


    dist = matrix_a % dist
    dist_nrows_a = dbcsr_distribution_nrows(dist)
    dist_ncols_a = dbcsr_distribution_ncols(dist)

    print *,"from carma"
    !print *,dist_nrows_a    
    !print *,dist_ncols_a
    
    !print *,size(matrix_a%row_blk_size%low%data) 
    !print *,size(matrix_a%col_blk_size%low%data) 

        
    group = dbcsr_mp_group(dbcsr_distribution_mp(dist))

    ! set subgrids keeping communication only along rows, communicating different columns (2nd dimension)   
    ! we consider 1 rank per row.
    ! TODO: any number of ranks per row
    free_coords(1) = .false.
    free_coords(2) = .true.

    call mp_cart_sub(group, free_coords, row_group)

    ! assume we have 1 rank per communicator
    loc_grid(1,1) = 0

    call dbcsr_mp_new(mp_env, loc_grid, row_group, 0) 
    
    ! create distribution for a matrix
    allocate(local_row_dist_a(dist_nrows_a),local_col_dist_a(dist_ncols_a))

    local_row_dist_a(:) = 0
    local_col_dist_a(:) = 0

    call dbcsr_distribution_new(loc_dist_a, mp_env, &
                                row_dist_block = local_row_dist_a, &
                                col_dist_block = local_col_dist_a, reuse_arrays=.TRUE.)


    ! create distribution for b matrix
    dist = matrix_b % dist
    dist_nrows_b = dbcsr_distribution_nrows(dist)
    dist_ncols_b = dbcsr_distribution_ncols(dist)

    allocate(local_row_dist_b(dist_nrows_b),local_col_dist_b(dist_ncols_b))

    local_row_dist_b(:) = 0
    local_col_dist_b(:) = 0

    call dbcsr_distribution_new(loc_dist_b, mp_env, &
                                row_dist_block = local_row_dist_b, &
                                col_dist_block = local_col_dist_b, reuse_arrays=.TRUE.)

    ! create distribution for c matrix    
    allocate(local_row_dist_c(dist_ncols_a),local_col_dist_c(dist_ncols_b))

    local_row_dist_c(:) = 0
    local_col_dist_c(:) = 0

    call dbcsr_distribution_new(loc_dist_c, mp_env, &
                                row_dist_block = local_row_dist_c, &
                                col_dist_block = local_col_dist_c, reuse_arrays=.TRUE.)
   

    !allocate(row_blk_size_c()) 
    !col_blk_size_c
    !print *,matrix_a%index(1:30) 
    
       
    ! create local matrix a
    call dbcsr_create(matrix = loc_matr_a, template = matrix_a, dist = loc_dist_a, data_buffer = matrix_a % data_area)
    
    allocate(loc_matr_a % index(size(matrix_a % index)))
    loc_matr_a % index(:) = matrix_a % index(:)
  
    call dbcsr_repoint_index(loc_matr_a)

    print *, "matr a matrix%row_blk_size",loc_matr_a % row_blk_size % low % data

    loc_matr_a % valid = .true.

    print *,"--------- local matrix a ---------------"
    call dbcsr_print(loc_matr_a)
     
    ! create local matrix b
    call dbcsr_create(matrix = loc_matr_b, template = matrix_b, dist = loc_dist_b, data_buffer = matrix_b % data_area)

    allocate(loc_matr_b % index(size(matrix_b % index)))
    loc_matr_b % index(:) = matrix_b % index(:)
    
    call dbcsr_repoint_index(loc_matr_b)
    
    print *,"--------- local matrix b ---------------"
    call dbcsr_print(loc_matr_b)
    
    call dbcsr_create(loc_matr_c, "loc_matr_c", loc_dist_c, dbcsr_type_no_symmetry, &
                      row_blk_size_c, col_blk_size_c, matrix_a % col_blk_size, matrix_b % col_blk_size, &
                      data_type = matrix_a % data_type)

    loc_matr_a % valid=.true.
    loc_matr_b % valid=.true.
    loc_matr_c % valid=.true.
    
    ! multiply local matrices
    CALL dbcsr_multiply('T', 'N', 1.0d0, loc_matr_a, loc_matr_b, 0.0d0, loc_matr_c)
   
    print *,"--------- local matrix c ---------------"
    call dbcsr_print(loc_matr_c)
   
    print *,"finita"
     

    call mp_comm_free(row_group)

    print *,"finita zwei"
    
    ! release the matrices
   CALL dbcsr_release(loc_matr_a)
   CALL dbcsr_release(loc_matr_b)
   CALL dbcsr_release(loc_matr_c)

   call dbcsr_mp_release(mp_env)
   CALL dbcsr_distribution_release(loc_dist_a)
   CALL dbcsr_distribution_release(loc_dist_b)
   CALL dbcsr_distribution_release(loc_dist_c)
    
  end subroutine
  
end module
