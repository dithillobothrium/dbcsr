
! Matrix-Matrix Multiply, C = AB
! Adapted from http://www.nersc.gov/users/software/programming-libraries/math-libraries/libsci/libsci-example/
!
! input:     ABCp.dat
!              n       number of rows/columns in matrix A
!              nb      matrix distribution block size   

! Run example:  export OMP_NUM_THREADS=2; mpirun -np 9 ./pdgemm_dbcsr

program matrix_matrix_multiply
  USE machine,                         ONLY: m_walltime,&
                                             default_output_unit
  USE dbcsr_api,                       ONLY: dbcsr_type_loc => dbcsr_type,&
                                             dbcsr_distribution_type,&
                                             dbcsr_distribution_new,&
                                             dbcsr_type_no_symmetry,&
                                             dbcsr_type_real_8,&
                                             dbcsr_create,&
                                             dbcsr_distribution_release,&
                                             convert_sizes_to_offsets,&
                                             dbcsr_get_info,&
                                             dbcsr_reserve_all_blocks,&
                                             dbcsr_iterator_type,&
                                             dbcsr_iterator_blocks_left,&
                                             dbcsr_iterator_next_block,&
                                             dbcsr_iterator_start,&
                                             dbcsr_iterator_stop,&
                                             dbcsr_multiply,&
                                             dbcsr_init_lib,&
                                             dbcsr_release,&
                                             dbcsr_finalize_lib,&
                                             dbcsr_get_data_p,&
                                             dbcsr_set_config
  USE message_passing,                 ONLY: mp_cart_create,&
                                             mp_comm_free,&
                                             mp_environ,&
                                             mp_world_finalize,&
                                             mp_world_init,&
                                             mp_sync,&
                                             mp_cart_rank,&
                                             mp_max
  USE kinds,                           ONLY: dp

  use carma, only: carma_mm
  
!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

  implicit none

  integer                             :: k, n, nb    ! problem size and block size
  integer                             :: myrows, mycols   ! size of local subset of global array
  integer                             :: i,j, igrid,jgrid, iproc,jproc, myi,myj, p
  REAL(kind=dp), dimension(:,:), allocatable :: myA, myB, myC_blacs, myC_dbcsr
  REAL(kind=dp)                              :: t1, t2, total_time_dbcsr

  integer                             :: numroc  ! blacs routine
  integer                             :: icontxt, prow, pcol, myrow, mycol  ! blacs data
  integer                             :: info    ! scalapack return value
  integer, dimension(9)               :: ides_a, ides_b, ides_c ! scalapack array desc

  INTEGER, DIMENSION(:, :), POINTER   :: pgrid
  INTEGER, DIMENSION(2)               :: coord, npdims, myploc ! MPI dimensions
  INTEGER                             :: mp_comm, mynode, numnodes, group ! MPI variables
  INTEGER                             :: io_unit, nthreads
  TYPE(dbcsr_type_loc)                    :: matA_dbcsr, matB_dbcsr, matC_dbcsr, matC_carma
  REAL(kind=dp)                       :: diff_val, min_val, max_val, sum_val

  INTEGER                             :: index, comm_thread_load

  LOGICAL, PARAMETER :: verbose=.FALSE.
  character(len=32) :: cdimx, cdimy
  

  CALL dbcsr_set_config(mm_stack_size=2000)

  !
  ! initialize mpi
  CALL mp_world_init(mp_comm)
  !
  ! setup the mp enviroment

  npdims(:) = 0

  if (command_argument_count() >=2) then
    call get_command_argument(1, cdimx)
    call get_command_argument(2, cdimy)
    read(cdimx,*) npdims(1)
    read(cdimy,*) npdims(2)
  endif
  
  
  CALL mp_cart_create(mp_comm, 2, npdims, myploc, group)
  CALL mp_environ(numnodes, mynode, group)
  

  prow = npdims(1)
  pcol = npdims(2)
 

  print *,"prow pcol", prow, pcol
  
  IF (mynode .EQ. 0) THEN
     nthreads = 1
!$OMP PARALLEL SHARED(nthreads)
!$OMP MASTER
!$   nthreads = omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL
     print *, "numnodes =", numnodes, "nthreads =", nthreads
  ENDIF
  if (verbose) then
     print *, 'mynode ', mynode, ' numnodes', numnodes, 'prow', prow, 'pcol', pcol
  endif
  ! Create processor grid
  ALLOCATE(pgrid(0:prow-1, 0:pcol-1))
  DO myrow = 0, prow-1
     DO mycol = 0, pcol-1
        coord = (/myrow, mycol/)
        CALL mp_cart_rank(group, coord, pgrid(myrow, mycol))
     ENDDO
  ENDDO
  !
  ! Read problem description
  open(unit=1,file="ABCp.dat",status="old",form="formatted")
  read(1,*) n
  read(1,*) k
  read(1,*) nb
  close(1)
  !
  if (((n/nb) < prow) .or. ((n/nb) < pcol)) then
     print *, "Problem size too small for processor set!"
     stop 100
  endif
  !
  ! Initialize blacs processor grid
  icontxt = group
  call blacs_gridinit(icontxt, 'R', prow, pcol)
  call blacs_gridinfo(icontxt, prow, pcol, myrow, mycol)
  if (verbose) then
     print *, "--------"
     print *, "Proc ",mynode,": myrow, mycol in p-array is ", &
              myrow, mycol
  endif
  !
  ! Construct local arrays
  ! Global structure:  matrix A of n rows and n columns
  !                    matrix B of n rows and n column
  !                    matrix C of n rows and n column
  myrows = numroc(n, nb, myrow, 0, prow)
  mycols = numroc(k, nb, mycol, 0, pcol)
  if (verbose) then
     print *, "Size of global array is ",n," x ",n
     print *, "Size of block is        ",nb," x ",nb
     print *, "Size of local array is  ",myrows," x ",mycols
  endif
  !
  ! Initialize local arrays    
  allocate(myA(myrows,mycols)) 
  allocate(myB(mycols,myrows)) 
  allocate(myC_blacs(myrows,myrows)) 
  allocate(myC_dbcsr(myrows,myrows)) 
  myC_blacs = 0
  myC_dbcsr = 0
  do i=1,n
     call g2l(i,n,prow,nb,iproc,myi)
     if (myrow==iproc) then
        do j=1,k
           call g2l(j,k,pcol,nb,jproc,myj)
           if (mycol==jproc) then
              myA(myi,myj) = real(i+j)/real(n)
              myB(myj,myi) = real(i-j)/real(n)
           endif
        enddo
     endif
  enddo
  !
  ! Prepare array descriptors for ScaLAPACK 
  ides_a(1) = 1         ! descriptor type
  ides_a(2) = icontxt   ! blacs context
  ides_a(3) = n         ! global number of rows
  ides_a(4) = k         ! global number of columns
  ides_a(5) = nb        ! row block size
  ides_a(6) = nb        ! column block size
  ides_a(7) = 0         ! initial process row
  ides_a(8) = 0         ! initial process column
  ides_a(9) = myrows   ! leading dimension of local array
  !
  do i=1,9
     ides_b(i) = ides_a(i)
     ides_c(i) = ides_a(i)
  enddo
  ides_b(3) = k
  ides_b(4) = n
  ides_b(9) = mycols
  ides_c(3) = n
  ides_c(4) = n
  ides_c(9) = myrows
  !
  ! Call ScaLAPACK library routine
  !
  ! Perform multiply
  CALL mp_sync(group)
  t1 = m_walltime()
  !call pdgemm('N','N', n, n, k, 1.0d0, myA, 1, 1, ides_a, &
  !            myB, 1, 1, ides_b, 0.0d0, &
  !            myC_blacs, 1, 1, ides_c)
  CALL mp_sync(group)
  t2 = m_walltime()
  !
  ! Print results
  if (mynode .EQ. 0) then
     print *, "PDGEMM time (s) =", t2-t1
     IF (verbose) THEN
!     call print_matrix(myA)
!     call print_matrix(myB)
        call print_matrix(myC_blacs)
     ENDIF
  endif
  !
  !***************************************************************************************
  !
  ! initialize libdbcsr
  CALL dbcsr_init_lib()
  total_time_dbcsr = 0
  !
  ! Convert BLACS to DBCSR matrices
  !
  CALL mp_sync(group)
  t1 = m_walltime()
  print *,"test321 myA"

  CALL copy_blacs_to_dbcsr_bc(myA, nb, nb, n, k, &
                              group, pgrid, matA_dbcsr)
  print *,"test321 myB"
  CALL copy_blacs_to_dbcsr_bc(myB, nb, nb, k, n, &
                              group, pgrid, matB_dbcsr)
  print *,"test321 myC"
  CALL copy_blacs_to_dbcsr_bc(myC_dbcsr, nb, nb, n, n, &
                              group, pgrid, matC_dbcsr)
  CALL mp_sync(group)
  t2 = m_walltime()
  total_time_dbcsr = total_time_dbcsr + t2 - t1
  !
  ! Print results
  if (mynode .EQ. 0) then
     print *, "Convert BLACS->DBCSR time (s) =", t2-t1
  endif
  !
  ! multiply the matrices
  CALL mp_sync(group)
  t1 = m_walltime()
   print *,"test321 CALL dbcsr_multiply('N', 'N', 1.0d0, matA_dbcsr, matB_dbcsr, 0.0d0, matC_dbcsr)"
  CALL dbcsr_multiply('N', 'N', 1.0d0, matA_dbcsr, matB_dbcsr, 0.0d0, matC_dbcsr)
  CALL mp_sync(group)
  t2 = m_walltime()
  total_time_dbcsr = total_time_dbcsr + t2 - t1
  !
  ! Print results
  if (mynode .EQ. 0) then
     print *, "Multiply DBCSR time (s) =", t2-t1
  endif






   call carma_mm(matA_dbcsr%prv, matB_dbcsr%prv, matC_dbcsr%prv)






  !
  ! Convert DBCSR to BLACS matrix
  !
  CALL mp_sync(group)
  t1 = m_walltime()
  print *,"test321 copy_dbcsr_to_blacs_bc(matC_dbcsr, myC_dbcsr)"
  CALL copy_dbcsr_to_blacs_bc(matC_dbcsr, myC_dbcsr)
  ! release the matrices
  CALL dbcsr_release(matA_dbcsr)
  CALL dbcsr_release(matB_dbcsr)
  CALL dbcsr_release(matC_dbcsr)
  CALL mp_sync(group)
  t2 = m_walltime()
  total_time_dbcsr = total_time_dbcsr + t2 - t1
  !
  ! Check results. Take maximum relative difference
  diff_val = 0
  do i=1,n
     call g2l(i,n,prow,nb,iproc,myi)
     if (myrow==iproc) then
        do j=1,n
           call g2l(j,n,pcol,nb,jproc,myj)
           if (mycol==jproc) then
              sum_val = abs(myC_blacs(myi,myj)+myC_dbcsr(myi,myj))
              if (sum_val .GT. 0.d0) then
                 diff_val = MAX(diff_val, abs(myC_blacs(myi,myj)-myC_dbcsr(myi,myj))/sum_val)
              endif
           endif
        enddo
     endif
  enddo
  call mp_max(diff_val, group)
  if (mynode .EQ. 0) then
     print *, "Convert DBCSR->BLACS time (s) =", t2-t1
     print *, "DBCSR time (s) =", total_time_dbcsr
     IF (verbose) THEN
        call print_matrix(myC_dbcsr)
     ENDIF
     print *, "Max Diff = ", diff_val
  endif
  !call print_matrix(myC_blacs)
  !call print_matrix(myC_dbcsr)
  !
  ! finalize libdbcsr
  io_unit = 0
  IF (mynode .EQ. 0) io_unit = default_output_unit
  CALL dbcsr_finalize_lib(mp_comm, io_unit)
 
  ! Deallocate the local arrays
  deallocate(myA, myB, myC_blacs, myC_dbcsr)
  deallocate(pgrid)

  ! End blacs for processors that are used
  call blacs_gridexit(icontxt)

  ! free comm
  CALL mp_comm_free(group)

  !
  ! finalize mpi
  CALL mp_world_finalize()

contains
  
  
  
  !
  ! convert global index to local index in block-cyclic distribution
  !
  subroutine g2l(i, n, np, nb, p, il)
    integer :: i    ! global array index, input
    integer :: n    ! global array dimension, input
    integer :: np   ! processor array dimension, input
    integer :: nb   ! block size, input
    integer :: p    ! processor array index, output
    integer :: il   ! local array index, output
    integer :: im1   

    im1 = i-1
    p   = mod((im1/nb),np)
    il  = (im1/(np*nb))*nb + mod(im1,nb) + 1
    return
  end subroutine g2l
  !
  ! convert local index to global index in block-cyclic distribution
  !
  subroutine l2g(il, p, n, np, nb, i)
    integer :: il   ! local array index, input
    integer :: p    ! processor array index, input
    integer :: n    ! global array dimension, input
    integer :: np   ! processor array dimension, input
    integer :: nb   ! block size, input
    integer :: i    ! global array index, output
    integer :: ilm1   

    ilm1 = il-1
    i    = (((ilm1/nb) * np) + p)*nb + mod(ilm1,nb) + 1
    return
  end subroutine l2g

  subroutine print_matrix(matrix)
    REAL(kind=dp), DIMENSION(:,:), INTENT(IN) :: matrix
    
    INTEGER :: i, n

    n = SIZE(matrix, 1)

    ! Column-major
    DO i= 1, n
       print *, matrix(i, 1:n)
     ENDDO
  end subroutine print_matrix

  subroutine print_matrix_dbcsr(matrix_dbcsr, lda)
    TYPE(dbcsr_type_loc), INTENT(IN) :: matrix_dbcsr
    INTEGER, INTENT(IN)          :: lda
    
    REAL(kind=dp), DIMENSION(:), POINTER       :: data_dbcsr
    REAL(kind=dp)                              :: select_data_type
    integer                             :: i

    data_dbcsr => dbcsr_get_data_p(matrix_dbcsr, select_data_type)
    ! Row-major
    DO i = 1, lda
       print *, data_dbcsr((i-1)*lda+1:i*lda)
     ENDDO
  end subroutine print_matrix_dbcsr

  !
  ! convert a BLACS matrix in a DBCSR matrix with block-cyclic distribution
  !
  SUBROUTINE copy_blacs_to_dbcsr_bc(fm, nrow_block, ncol_block, nrow_global, ncol_global, &
                                    group, pgrid, bc_mat)
    REAL(kind=dp), dimension(:,:), allocatable :: fm ! full matrix
    INTEGER, INTENT(IN)                 :: nrow_block, ncol_block, nrow_global, ncol_global, group
    INTEGER, DIMENSION(:, :), POINTER   :: pgrid    
    TYPE(dbcsr_type_loc), INTENT(INOUT)     :: bc_mat
 
    INTEGER                             :: col, row
    INTEGER, ALLOCATABLE, DIMENSION(:)  :: first_col, first_row, last_col, last_row
    INTEGER, DIMENSION(:), POINTER      :: col_blk_size, row_blk_size
    REAL(kind=dp), DIMENSION(:, :), POINTER    :: dbcsr_block
    TYPE(dbcsr_distribution_type)       :: bc_dist
    TYPE(dbcsr_iterator_type)           :: iter

    ! Create a block-cyclic distribution compatible with the FM matrix.
    NULLIFY (col_blk_size, row_blk_size)
    CALL dbcsr_create_dist_block_cyclic(bc_dist, &
         nrows=nrow_global, ncolumns=ncol_global, & ! Actual full matrix size
         nrow_block=nrow_block, ncol_block=ncol_block, & ! BLACS parameters
         group=group, pgrid=pgrid, &
         row_blk_sizes=row_blk_size, col_blk_sizes=col_blk_size) ! block-cyclic row/col sizes
    ! Create the block-cyclic DBCSR matrix
    CALL dbcsr_create(bc_mat, "Block-cyclic ", bc_dist, &
         dbcsr_type_no_symmetry, row_blk_size, col_blk_size, nze=0, &
         reuse_arrays=.TRUE., data_type=dbcsr_type_real_8)
    CALL dbcsr_distribution_release(bc_dist)
    ! allocate all blocks
    CALL dbcsr_reserve_all_blocks(bc_mat)
    
    CALL calculate_fm_block_ranges(bc_mat, first_row, last_row, first_col, last_col)
    !
    ! Copy the FM data to the block-cyclic DBCSR matrix.  This step
    ! could be skipped with appropriate DBCSR index manipulation.
    !$OMP PARALLEL DEFAULT(NONE) PRIVATE(iter, row, col, dbcsr_block) &
    !$OMP SHARED(bc_mat, last_row, first_row, last_col, first_col, fm)
    CALL dbcsr_iterator_start(iter, bc_mat)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, dbcsr_block)
       dbcsr_block(:, :) = fm(first_row(row):last_row(row), first_col(col):last_col(col))
    ENDDO
    CALL dbcsr_iterator_stop(iter)
    !$OMP END PARALLEL

  END SUBROUTINE copy_blacs_to_dbcsr_bc
  !
  ! Copy a DBCSR_BLACS matrix to a BLACS matrix
  !
  SUBROUTINE copy_dbcsr_to_blacs_bc(bc_mat, fm)
    TYPE(dbcsr_type_loc), INTENT(IN)                    :: bc_mat
    REAL(kind=dp), DIMENSION(:, :), ALLOCATABLE            :: fm

    INTEGER :: col, handle, row
    INTEGER, ALLOCATABLE, DIMENSION(:)              :: first_col, first_row, last_col, last_row
    REAL(kind=dp), DIMENSION(:, :), POINTER                :: dbcsr_block
    TYPE(dbcsr_iterator_type)                       :: iter

    CALL calculate_fm_block_ranges(bc_mat, first_row, last_row, first_col, last_col)

    ! Now copy data to the FM matrix
    fm = 0.0d0
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(iter, row, col, dbcsr_block) &
!$OMP SHARED(bc_mat, last_row, first_row, last_col, first_col, fm)
    CALL dbcsr_iterator_start(iter, bc_mat)
    DO WHILE (dbcsr_iterator_blocks_left(iter))
       CALL dbcsr_iterator_next_block(iter, row, col, dbcsr_block)
       fm(first_row(row):last_row(row), first_col(col):last_col(col)) = dbcsr_block(:, :)
    ENDDO
    CALL dbcsr_iterator_stop(iter)
!$OMP END PARALLEL

  END SUBROUTINE copy_dbcsr_to_blacs_bc

  !
  ! Creates a block-cyclic compatible distribution
  !
  SUBROUTINE dbcsr_create_dist_block_cyclic (dist, nrows, ncolumns,&
       nrow_block, ncol_block, group, pgrid, row_blk_sizes, col_blk_sizes)
      TYPE(dbcsr_distribution_type), INTENT(OUT)          :: dist
      INTEGER, INTENT(IN)                                :: nrows, ncolumns, nrow_block, ncol_block
      INTEGER, INTENT(IN)                                :: group
      INTEGER, DIMENSION(:, :), POINTER                  :: pgrid
      INTEGER, DIMENSION(:), INTENT(OUT), POINTER        :: row_blk_sizes, col_blk_sizes

      INTEGER                                            :: nblkcols, nblkrows, npcols, nprows, &
                                                            pdim, sz
      INTEGER, DIMENSION(:), POINTER                     :: cd_data, rd_data

    ! Row sizes
    IF(nrow_block.EQ.0) THEN
       nblkrows = 0
       sz = 0
    ELSE
       nblkrows = nrows / nrow_block
       sz = MOD(nrows, nrow_block)
    ENDIF
    IF (sz .GT. 0) nblkrows = nblkrows + 1
    ALLOCATE (row_blk_sizes(nblkrows), rd_data(nblkrows))
    row_blk_sizes = nrow_block
    IF (sz .NE. 0) row_blk_sizes(nblkrows) = sz

    ! Column sizes
    IF(ncol_block.EQ.0) THEN
       nblkcols = 0
       sz = 0
    ELSE
       nblkcols = ncolumns / ncol_block
       sz = MOD(ncolumns, ncol_block)
    ENDIF
    IF (sz .GT. 0) nblkcols = nblkcols + 1
    ALLOCATE (col_blk_sizes(nblkcols), cd_data(nblkcols))
    col_blk_sizes = ncol_block
    IF (sz .NE. 0) col_blk_sizes(nblkcols) = sz
    !
    ! Calculate process row distribution
    nprows = SIZE(pgrid, 1)
    DO pdim = 0 , MIN(nprows-1,nblkrows-1)
       rd_data(1+pdim:nblkrows:nprows) = pdim
    END DO
    ! Calculate process column distribution
    npcols = SIZE(pgrid, 2)
    DO pdim = 0 , MIN(npcols-1, nblkcols-1)
       cd_data(1+pdim:nblkcols:npcols) = pdim
    END DO
    !
    CALL dbcsr_distribution_new (dist,&
                                 group=group, pgrid=pgrid,&
                                 row_dist=rd_data,&
                                 col_dist=cd_data,&
                                 reuse_arrays=.TRUE.)

  END SUBROUTINE dbcsr_create_dist_block_cyclic
  !
  ! Helper routine used to copy blocks from DBCSR into FM matrices and vice versa
  !
  SUBROUTINE calculate_fm_block_ranges(bc_mat, first_row, last_row, first_col, last_col)
    TYPE(dbcsr_type_loc), INTENT(IN)                    :: bc_mat
    INTEGER :: col, nblkcols_local, nblkcols_total, nblkrows_local, nblkrows_total, row
    INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_col, first_row, last_col, &
         last_row, local_col_sizes, &
         local_row_sizes
    INTEGER, DIMENSION(:), POINTER                     :: col_blk_size, local_cols, local_rows, &
         row_blk_size

    CALL dbcsr_get_info(bc_mat, &
         nblkrows_total=nblkrows_total, &
         nblkcols_total=nblkcols_total, &
         nblkrows_local=nblkrows_local, &
         nblkcols_local=nblkcols_local, &
         local_rows=local_rows, &
         local_cols=local_cols, &
         row_blk_size=row_blk_size, &
         col_blk_size=col_blk_size )

    ! calculate first_row and last_row
    ALLOCATE (local_row_sizes(nblkrows_total))
    local_row_sizes(:) = 0
    IF (nblkrows_local .GE. 1) THEN
       DO row = 1, nblkrows_local
          local_row_sizes(local_rows(row)) = row_blk_size(local_rows(row))
       END DO
    ENDIF
    ALLOCATE (first_row(nblkrows_total), last_row(nblkrows_total))
    CALL convert_sizes_to_offsets(local_row_sizes, first_row, last_row)

    ! calculate first_col and last_col
    ALLOCATE (local_col_sizes(nblkcols_total))
    local_col_sizes(:) = 0
    IF (nblkcols_local .GE. 1) THEN
       DO col = 1, nblkcols_local
          local_col_sizes(local_cols(col)) = col_blk_size(local_cols(col))
       END DO
    ENDIF
    ALLOCATE (first_col(nblkcols_total), last_col(nblkcols_total))
    CALL convert_sizes_to_offsets(local_col_sizes, first_col, last_col)

  END SUBROUTINE calculate_fm_block_ranges

end program matrix_matrix_multiply
