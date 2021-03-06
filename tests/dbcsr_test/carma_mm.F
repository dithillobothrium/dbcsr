

module carma
 USE machine,                         ONLY: m_walltime,&
                                             default_output_unit
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

     USE dbcsr_types

   USE dbcsr_mp_methods,                ONLY: dbcsr_mp_group,&
                                              dbcsr_mp_npcols,&
                                              dbcsr_mp_nprows,&
                                              dbcsr_mp_numnodes,&
                                              dbcsr_mp_new,  &
                                              dbcsr_mp_myprow

  USE message_passing,                 ONLY: mp_cart_create,&
                                             mp_comm_free,&
                                             mp_environ,&
                                             mp_world_finalize,&
                                             mp_world_init,&
                                             mp_sync,&
                                             mp_cart_rank,&
                                             mp_max, &
                                             mp_cart_sub, mp_waitall,  mp_irecv, mp_isend, &
                                             mp_sum, mp_allgather
    USE dbcsr_data_methods                                     
    use dbcsr_work_operations

    use dbcsr_io, only: dbcsr_print

    use dbcsr_index_operations
    use dbcsr_multiply_api
    
    use dbcsr_operations
    
    use dbcsr_mp_operations,  only:  dbcsr_irecv_any, dbcsr_isend_any
    
    use dbcsr_util
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
    
    type(dbcsr_type)                :: loc_matr_a, loc_matr_b, loc_matr_c, loc_matr_c_store
                                    
    
    type(dbcsr_type), pointer       :: loc_matr_c_recv, loc_matr_c_send, tmp_matr_ptr
    
    type(dbcsr_type), target        :: loc_matr_c_1, loc_matr_c_2
    
    type(dbcsr_mp_obj)              :: mp_env,  total_mp_env

    type(dbcsr_data_obj)            :: send_nze_buffer

    integer                     :: tick, max_sizes(2),  requests(4),  proc_recv,  proc_send, num_rows
    integer(8)                  :: sum_nze_nblks(2)
    real(8)                     :: alpha, beta, t1
    integer, dimension(:), pointer            :: index_sizes, data_sizes

    
    ! for summation
    alpha = 1.0
    beta= 1.0
    
    ! take distribution and communicator
    dist = matrix_a % dist
    dist_nrows_a = dbcsr_distribution_nrows(dist)
    dist_ncols_a = dbcsr_distribution_ncols(dist)

!     print *,"from carma"
    !print *,dist_nrows_a    
    !print *,dist_ncols_a
    
    !print *,size(matrix_a%row_blk_size%low%data) 
    !print *,size(matrix_a%col_blk_size%low%data) 

        
    total_mp_env = dbcsr_distribution_mp(dist)
    group = dbcsr_mp_group(total_mp_env)
    
    ! for communication
    num_rows = dbcsr_mp_nprows(total_mp_env)
    allocate (index_sizes(num_rows), data_sizes(num_rows))
    
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
    ! TODO , make_index=.false.
    call dbcsr_create(matrix = loc_matr_a, template = matrix_a, dist = loc_dist_a, &
                      data_buffer = matrix_a % data_area)
    
    allocate(loc_matr_a % index(size(matrix_a % index)))
    loc_matr_a % index(:) = matrix_a % index(:)
  
    !print *,"sizes   ",loc_matr_a % index(dbcsr_slot_nfullcols_local), loc_matr_a % index(dbcsr_slot_nfullrows_local)
    
    loc_matr_a % index(dbcsr_slot_nfullcols_total) =  loc_matr_a % index(dbcsr_slot_nfullcols_local)
    loc_matr_a % index(dbcsr_slot_nfullrows_total) =  loc_matr_a % index(dbcsr_slot_nfullrows_local)
    
    call dbcsr_repoint_index(loc_matr_a)

    loc_matr_a % nfullrows_total = loc_matr_a % index(dbcsr_slot_nfullrows_local) 
    loc_matr_a % nfullcols_total = loc_matr_a % index(dbcsr_slot_nfullcols_local)
    
    !print *,"asdfsdf", loc_matr_a % nfullrows_local , loc_matr_a % nfullcols_local 
!     print *, "matr a matrix%row_blk_size",loc_matr_a % row_blk_size % low % data

    loc_matr_a % valid = .true.

!     print *,"--------- local matrix a ---------------"
!     call dbcsr_print(loc_matr_a)
     
    ! create local matrix b
    call dbcsr_create(matrix = loc_matr_b, template = matrix_b, dist = loc_dist_b, &
                      data_buffer = matrix_b % data_area)

    allocate(loc_matr_b % index(size(matrix_b % index)))
    loc_matr_b % index(:) = matrix_b % index(:)
    
    loc_matr_b % index(dbcsr_slot_nfullcols_total) =  loc_matr_b % index(dbcsr_slot_nfullcols_local)
    loc_matr_b % index(dbcsr_slot_nfullrows_total) =  loc_matr_b % index(dbcsr_slot_nfullrows_local)
    
    call dbcsr_repoint_index(loc_matr_b)
    
    loc_matr_b % nfullrows_total = loc_matr_b % index(dbcsr_slot_nfullrows_local) 
    loc_matr_b % nfullcols_total = loc_matr_b % index(dbcsr_slot_nfullcols_local)
    
!     print *,"--------- local matrix b ---------------"
!     call dbcsr_print(loc_matr_b)
    
    call dbcsr_create(loc_matr_c, "loc_matr_c", loc_dist_c, dbcsr_type_no_symmetry, &
                      row_blk_size_c, col_blk_size_c, matrix_a % col_blk_size, matrix_b % col_blk_size, &
                      data_type = matrix_a % data_type)

    loc_matr_a % valid=.true.
    loc_matr_b % valid=.true.
    loc_matr_c % valid=.true.
    
    
    !call dbcsr_finalize(loc_matr_c)
    
    !print *,"occ a ,b", dbcsr_get_occupation(loc_matr_a), dbcsr_get_occupation(loc_matr_b)
    
    
    CALL mp_sync(group)
    t1 = m_walltime()
    
    
    ! multiply local matrices
    CALL dbcsr_multiply('T', 'N', 1.0d0, loc_matr_a, loc_matr_b, 0.0d0, loc_matr_c)
   
    !call dbcsr_finalize(loc_matr_c)

   CALL mp_sync(group)
   

   !print *,"+++++ carma dbcsr_multiply :",dbcsr_mp_myprow(total_mp_env), m_walltime()-t1
   
     print *,"--------- local matrix c ---------------"
!     call dbcsr_print(loc_matr_c)
    print *,"index size", loc_matr_c % index(1)
    print *,"nze", loc_matr_c % nze
    
    max_sizes(1) = loc_matr_c % nze
    max_sizes(2) = loc_matr_c % index(1)
    sum_nze_nblks(1) = loc_matr_c % nze
    sum_nze_nblks(2) = loc_matr_c % nblks

    call mp_max(max_sizes, group)
    call mp_sum(sum_nze_nblks, group)
    call mp_allgather(loc_matr_c % nze, data_sizes, group)
    call mp_allgather(loc_matr_c % index(1), index_sizes, group)



    !print *,"index ",loc_matr_c % index

    if (dbcsr_mp_myprow(total_mp_env) == 0) then
        print *,"-- max sizes", max_sizes
        print *,"-- sum b nze, blocks", sum_nze_nblks
            print *,"data_sizes ",data_sizes
    print *,"index_sizes ",index_sizes
    endif
    
    ! TODO , make_index=.false.
    ! create send matrix
    call dbcsr_create(matrix = loc_matr_c_1, template = loc_matr_c, nze = max_sizes(1))
    
    allocate(loc_matr_c_1 % index( max_sizes(2) ))
    
    call dbcsr_repoint_index(loc_matr_c_1)
    loc_matr_c_1 % valid=.true.
    
    call dbcsr_copy(loc_matr_c_1, loc_matr_c)

        
    ! recv matrix
    call dbcsr_create(matrix = loc_matr_c_2, name = "loc_matr_c_2" ,template = loc_matr_c, nze = max_sizes(1) )
    allocate(loc_matr_c_2 % index( max_sizes(2) ))
    
    
    loc_matr_c_send => loc_matr_c_1
    loc_matr_c_recv => loc_matr_c_2
    
    
    proc_recv = mod(dbcsr_mp_myprow(total_mp_env) + 1, dbcsr_mp_nprows(total_mp_env) )
    proc_send = dbcsr_mp_myprow(total_mp_env) - 1
    
    if (proc_send < 0) proc_send = dbcsr_mp_nprows(total_mp_env) - 1
        
    print *,"send proc, recv prc ",     proc_send,  proc_recv
    
    call dbcsr_data_init(send_nze_buffer)
    call dbcsr_data_new(send_nze_buffer, loc_matr_c_send % data_type)
        
    !call dbcsr_print(loc_matr_c)
    ! test batched add
    call dbcsr_batched_add_anytype_begin(loc_matr_c, loc_matr_c_send,  alpha, beta, sum_nze_nblks(1), sum_nze_nblks(2))

    do tick = 1 , dbcsr_mp_nprows(total_mp_env) 
        
        if (tick > 1) then
            call mp_waitall(requests)
            call mp_sync(group)
            call dbcsr_repoint_index(loc_matr_c_send)
            loc_matr_c_send % valid = .true.
            !print *,"received index", loc_matr_c_send % index
        endif
        
        if (tick < dbcsr_mp_nprows(total_mp_env) ) then
            ! recv send
            CALL dbcsr_irecv_any(loc_matr_c_recv % data_area, proc_recv, group, requests(1),0)
            
            CALL mp_irecv(loc_matr_c_recv % index, proc_recv, group, requests(2),1)
            


            CALL dbcsr_data_set_pointer(area = send_nze_buffer, &
                                        rsize = loc_matr_c_send % nze, &
                                        csize = 1, &
                                        pointee = loc_matr_c_send % data_area)
                        
            CALL dbcsr_isend_any(send_nze_buffer, proc_send, group, requests(3),0)
            
            CALL mp_isend(loc_matr_c_send % index(1:loc_matr_c_send % index(1)), proc_send, group, requests(4),1)
        endif
        
        if (tick > 1) then
            
         print *," ^^^^^^^^^^^ tick ", tick 
         print *, " ##### send matrix matrix #########"
         !call dbcsr_print(loc_matr_c_send)
!         print *, "#### receive matrix #####"
!         call dbcsr_print(loc_matr_c_recv)
            
            call dbcsr_batched_add_anytype(loc_matr_c, loc_matr_c_send, beta )

            
            
            !print *,"result matrix"
            !call dbcsr_print(loc_matr_c)
            
        endif
    
        tmp_matr_ptr => loc_matr_c_send
        loc_matr_c_send => loc_matr_c_recv
        loc_matr_c_recv => tmp_matr_ptr
        
    enddo
    
    ! finalize matrix after add
    call dbcsr_batched_add_anytype_end(loc_matr_c)

!     print *, "==================== result matrix ================="
    !call dbcsr_print(loc_matr_c)
    print *,"checksum carma", dbcsr_checksum(loc_matr_c)
    
    call mp_comm_free(row_group)

    
    ! release the matrices
   CALL dbcsr_release(loc_matr_c_1)
   CALL dbcsr_release(loc_matr_c_2)
   
   CALL dbcsr_release(loc_matr_a)
   CALL dbcsr_release(loc_matr_b)
   CALL dbcsr_release(loc_matr_c)

   call dbcsr_mp_release(mp_env)
   CALL dbcsr_distribution_release(loc_dist_a)
   CALL dbcsr_distribution_release(loc_dist_b)
   CALL dbcsr_distribution_release(loc_dist_c)
    
  end subroutine
  
end module