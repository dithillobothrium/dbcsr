!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Common variables and routines for the dbcsr matrix-matrix multiplication algorithms.
!>
!> \author  Alfio Lazzaro
!>
!> <b>Modification history:</b>
!>  - 2016-08    Code organization (Alfio Lazzaro).
! **************************************************************************************************
MODULE dbcsr_mm_common
   USE acc_event,                       ONLY: acc_event_record,&
                                              acc_event_synchronize,&
                                              acc_stream_wait_event
   USE acc_stream,                      ONLY: acc_stream_type
   USE array_types,                     ONLY: array_data,&
                                              array_hold
   USE dbcsr_acc_operations,            ONLY: dbcsr_acc_transpose
   USE dbcsr_data_methods,              ONLY: dbcsr_data_ensure_size,&
                                              dbcsr_data_get_size,&
                                              dbcsr_data_host2dev,&
                                              dbcsr_data_set_size_referenced,&
                                              dbcsr_get_data_p_c,&
                                              dbcsr_get_data_p_d,&
                                              dbcsr_get_data_p_s,&
                                              dbcsr_get_data_p_z
   USE dbcsr_methods,                   ONLY: dbcsr_get_data_type,&
                                              dbcsr_get_index_memory_type,&
                                              dbcsr_nfullcols_local,&
                                              dbcsr_nfullrows_local,&
                                              dbcsr_valid_index
   USE dbcsr_mm_multrec,                ONLY: dbcsr_mm_multrec_type
   USE dbcsr_ptr_util,                  ONLY: ensure_array_size
   USE dbcsr_types,                     ONLY: &
        dbcsr_data_obj, dbcsr_memtype_type, dbcsr_mpi_size_limits, dbcsr_mpi_statistics_type, &
        dbcsr_type, dbcsr_type_complex_4, dbcsr_type_complex_8, dbcsr_type_int_4, &
        dbcsr_type_real_4, dbcsr_type_real_8
   USE dbcsr_work_operations,           ONLY: dbcsr_create
   USE kinds,                           ONLY: int_4,&
                                              int_8,&
                                              real_4,&
                                              real_8,&
                                              sp
#include "../../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_common'

   TYPE dbcsr_memtype_type_p
      TYPE(dbcsr_memtype_type), POINTER :: p => Null()
      ! ensure that array-elements are on different cache lines
      INTEGER(kind=int_4), DIMENSION(64)    :: padding
   END TYPE dbcsr_memtype_type_p

   TYPE(dbcsr_memtype_type_p), DIMENSION(:), POINTER, SAVE :: memtype_product_wm => Null()

   TYPE(dbcsr_mpi_statistics_type), SAVE :: dbcsr_mpi_statistics
   INTEGER, SAVE :: num_multiplications = 0
   REAL, SAVE :: max_memory = 0
   REAL, PARAMETER :: huge_norm = HUGE(1.0)**(1.0/3.0)

   TYPE(dbcsr_memtype_type), SAVE  :: memtype_abpanel_1, memtype_abpanel_2, &
                                      memtype_trsbuffer_1, memtype_trsbuffer_2, &
                                      memtype_mpi_buffer, memtype_mpi_product
   TYPE(acc_stream_type), SAVE     :: stream_1, stream_2
   ! ab-panels and streams are shared between all threads

   TYPE dbcsr_mm_multrec_type_p
      TYPE(dbcsr_mm_multrec_type), POINTER :: p => Null()
      ! ensure that array-elements are on different cache lines
      INTEGER(kind=int_4), DIMENSION(64)       :: padding
   END TYPE dbcsr_mm_multrec_type_p

   PUBLIC :: memtype_product_wm
   PUBLIC :: dbcsr_mpi_statistics, num_multiplications
   PUBLIC :: max_memory

   PUBLIC :: memtype_abpanel_1, memtype_abpanel_2, &
             memtype_trsbuffer_1, memtype_trsbuffer_2, &
             memtype_mpi_buffer, memtype_mpi_product
   PUBLIC :: stream_1, stream_2

   PUBLIC :: dbcsr_mm_multrec_type_p
   PUBLIC :: count_mpi_statistics
   PUBLIC :: setup_buffer_matrix
   PUBLIC :: rec_sort_index
   PUBLIC :: enumerate_blk_sizes

   PUBLIC :: acc_transpose_blocks

   PUBLIC :: product_matrix_size_guess

   PUBLIC :: calculate_norms
   PUBLIC :: huge_norm
   PUBLIC :: local_filter

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mpi_statistics ...
!> \param data_size ...
!> \param element_size_bytes ...
!> \param size_breakdown ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   SUBROUTINE count_mpi_statistics(mpi_statistics, data_size, &
                                   element_size_bytes, size_breakdown)
      REAL, DIMENSION(:), INTENT(INOUT)                  :: mpi_statistics
      INTEGER, INTENT(IN)                                :: data_size
      INTEGER, INTENT(IN)                                :: element_size_bytes
      INTEGER(KIND=int_8), DIMENSION(:, :), &
         INTENT(INOUT), OPTIONAL                         :: size_breakdown

      INTEGER                                            :: ilimit, nlimits
      INTEGER(KIND=int_8)                                :: data_size_bytes, llimit

      ! change in bytes
      data_size_bytes = INT(data_size, KIND=int_8)*INT(element_size_bytes, KIND=int_8)
      !
      mpi_statistics(1) = mpi_statistics(1)+REAL(data_size_bytes)
      mpi_statistics(2) = MIN(mpi_statistics(2), REAL(data_size_bytes))
      mpi_statistics(3) = MAX(mpi_statistics(3), REAL(data_size_bytes))
      IF (PRESENT(size_breakdown)) THEN
         nlimits = SIZE(dbcsr_mpi_size_limits)
         ! check for oversize messages
         IF (data_size_bytes .GT. dbcsr_mpi_size_limits(nlimits)) THEN
            size_breakdown(nlimits+1, 1) = size_breakdown(nlimits+1, 1)+1
            size_breakdown(nlimits+1, 2) = size_breakdown(nlimits+1, 2)+data_size_bytes
            RETURN
         ENDIF
         llimit = 0
         DO ilimit = 1, nlimits
            IF (data_size_bytes .GE. llimit .AND. data_size_bytes .LE. dbcsr_mpi_size_limits(ilimit)) THEN
               size_breakdown(ilimit, 1) = size_breakdown(ilimit, 1)+1
               size_breakdown(ilimit, 2) = size_breakdown(ilimit, 2)+data_size_bytes
               RETURN
            ENDIF
            llimit = dbcsr_mpi_size_limits(ilimit)
         ENDDO
      ENDIF
   END SUBROUTINE count_mpi_statistics

! **************************************************************************************************
!> \brief ...
!> \param matrix ...
!> \param source_matrix ...
!> \param index_size ...
!> \param data_size ...
!> \param data_buffer ...
!> \param data_memory_type ...
! **************************************************************************************************
   SUBROUTINE setup_buffer_matrix(matrix, source_matrix, &
                                  index_size, data_size, data_buffer, data_memory_type)
      TYPE(dbcsr_type), INTENT(INOUT)                    :: matrix
      TYPE(dbcsr_type), INTENT(IN)                       :: source_matrix
      INTEGER, INTENT(IN), OPTIONAL                      :: index_size, data_size
      TYPE(dbcsr_data_obj), INTENT(IN), OPTIONAL         :: data_buffer
      TYPE(dbcsr_memtype_type), INTENT(IN), OPTIONAL     :: data_memory_type

      matrix = dbcsr_type()
      CALL dbcsr_create(matrix, &
                        template=source_matrix, &
                        name=TRIM("Buffer of "//TRIM(source_matrix%name)), &
                        nze=data_size, &
                        data_buffer=data_buffer, &
                        data_memory_type=data_memory_type, &
                        index_memory_type=memtype_mpi_buffer)
      IF (PRESENT(data_size)) THEN
         CALL dbcsr_data_ensure_size( &
            matrix%data_area, &
            data_size, nocopy=.TRUE.)
      ENDIF
      IF (PRESENT(index_size)) THEN
         CALL ensure_array_size( &
            matrix%index, &
            ub=index_size, nocopy=.TRUE., &
            memory_type=dbcsr_get_index_memory_type(matrix))
      ENDIF
      matrix%negate_real = source_matrix%negate_real
      matrix%negate_imaginary = source_matrix%negate_imaginary
      matrix%local_indexing = source_matrix%local_indexing
      matrix%list_indexing = source_matrix%list_indexing
      !
      IF (source_matrix%has_local_rows) THEN
         matrix%local_rows = source_matrix%local_rows
         CALL array_hold(matrix%local_rows)
         matrix%has_local_rows = .TRUE.
      ENDIF
      IF (source_matrix%has_global_rows) THEN
         matrix%global_rows = source_matrix%global_rows
         CALL array_hold(matrix%global_rows)
         matrix%has_global_rows = .TRUE.
      ENDIF
      IF (source_matrix%has_local_cols) THEN
         matrix%local_cols = source_matrix%local_cols
         CALL array_hold(matrix%local_cols)
         matrix%has_local_cols = .TRUE.
      ENDIF
      IF (source_matrix%has_global_cols) THEN
         matrix%global_cols = source_matrix%global_cols
         CALL array_hold(matrix%global_cols)
         matrix%has_global_cols = .TRUE.
      ENDIF

   END SUBROUTINE setup_buffer_matrix

! **************************************************************************************************
!> \brief Sorts index for recursing.
!> \param mi ...
!> \param mf ...
!> \param ni ...
!> \param nf ...
!> \param nele ...
!> \param a ...
!> \param d ...
!> \par History
!> - 2011-02-17 [UB] modified for use in DBCSR; reduced memory usage.
!> \author JV
!> \note Always cut longest first. On a tie cut N
! **************************************************************************************************
   RECURSIVE SUBROUTINE rec_sort_index(mi, mf, ni, nf, nele, a, d)
      INTEGER, INTENT(IN)                                :: mi, mf, ni, nf, nele
      INTEGER, DIMENSION(3, 1:nele), INTENT(inout)       :: a
      INTEGER, INTENT(IN)                                :: d

      LOGICAL, PARAMETER                                 :: dbg = .FALSE.

      INTEGER                                            :: half, M, N, nlow
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: tmp

!   ---------------------------------------------------------------------------

      IF (dbg) THEN
         WRITE (*, *) " rs", mi, mf, "/", ni, nf, "=>", nele, d
         WRITE (*, '(3(1X,I7))') a(:, 1:nele)
      ENDIF
      IF (dbg) THEN
         IF (d .GT. 20) THEN
            WRITE (*, *) a(1, -d*1000)
         ENDIF
      ENDIF
      ALLOCATE (tmp(3, nele))
      M = mf-mi+1
      N = nf-ni+1
      IF (M > N) THEN
         half = M/2
         CALL rec_split(nele, a, tmp, 1, nlow, mi, half)
         a = tmp
         DEALLOCATE (tmp)
         IF (nlow .GT. 1) THEN
            CALL rec_sort_index(mi, mi+half-1, ni, nf, nlow, a(:, 1:nlow), d+1)
         ENDIF
         IF (nele-nlow .GT. 1) THEN
            CALL rec_sort_index(mi+half, mf, ni, nf, nele-nlow, a(:, nlow+1:nele), d+1)
         ENDIF
      ELSE
         half = N/2
         CALL rec_split(nele, a, tmp, 2, nlow, ni, half)
         a = tmp
         DEALLOCATE (tmp)
         IF (nlow .GT. 1) THEN
            CALL rec_sort_index(mi, mf, ni, ni+half-1, nlow, a(:, 1:nlow), d+1)
         ENDIF
         IF (nele-nlow .GT. 1) THEN
            CALL rec_sort_index(mi, mf, ni+half, nf, nele-nlow, a(:, nlow+1:nele), d+1)
         ENDIF
      ENDIF
   END SUBROUTINE rec_sort_index

! **************************************************************************************************
!> \brief ...
!> \param nele ...
!> \param a ...
!> \param split ...
!> \param row_or_col ...
!> \param nlow ...
!> \param mi ...
!> \param half ...
! **************************************************************************************************
   SUBROUTINE rec_split(nele, a, split, row_or_col, nlow, mi, half)
      INTEGER, INTENT(IN)                                :: nele
      INTEGER, DIMENSION(3, nele), INTENT(IN)            :: a
      INTEGER, DIMENSION(3, nele), INTENT(OUT)           :: split
      INTEGER, INTENT(IN)                                :: row_or_col
      INTEGER, INTENT(OUT)                               :: nlow
      INTEGER, INTENT(IN)                                :: mi, half

      CHARACTER(LEN=*), PARAMETER :: routineN = 'rec_split', routineP = moduleN//':'//routineN

      INTEGER                                            :: el, half_m, p_high, p_low

      half_m = mi+half-1
      p_low = 1
      p_high = nele
      DO el = 1, nele
         IF (a(row_or_col, el) <= half_m) THEN
            split(1:3, p_low) = a(1:3, el)
            p_low = p_low+1
         ELSE
            split(1:3, p_high) = a(1:3, el)
            p_high = p_high-1
         ENDIF
      ENDDO
      nlow = p_low-1
      CPASSERT(p_high .EQ. nlow)

   END SUBROUTINE rec_split

! **************************************************************************************************
!> \brief Enumerate all occuring blocksizes
!> \param blk_sizes ...
!> \param max_size ...
!> \param enum ...
!> \param rev_enum ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE enumerate_blk_sizes(blk_sizes, max_size, enum, rev_enum)
      INTEGER, DIMENSION(:), POINTER                     :: blk_sizes
      INTEGER, INTENT(IN)                                :: max_size
      INTEGER, DIMENSION(:), POINTER                     :: enum, rev_enum

      CHARACTER(len=*), PARAMETER :: routineN = 'enumerate_blk_sizes', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, n

      CALL timeset(routineN, handle)

      ALLOCATE (enum(0:max_size))
      enum(:) = 0

      DO i = 1, SIZE(blk_sizes)
         enum(blk_sizes(i)) = 1
      ENDDO

      n = SUM(enum)
      ALLOCATE (rev_enum(n))

      n = 0
      DO i = 0, SIZE(enum)-1
         IF (enum(i) > 0) THEN
            n = n+1
            enum(i) = n
            rev_enum(n) = i
         END IF
      ENDDO

      CALL timestop(handle)

   END SUBROUTINE enumerate_blk_sizes

! **************************************************************************************************
!> \brief write out a stack for transposing the blocks
!> \param matrix ...
!> \param trs_stackbuf ...
!> \param row_blk_sizes ...
!> \param col_blk_sizes ...
!> \param row_blk_sizes2enum ...
!> \param enum2row_blk_sizes ...
!> \param col_blk_sizes2enum ...
!> \param enum2col_blk_sizes ...
!> \param noresize ...
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE acc_transpose_blocks(matrix, trs_stackbuf, &
                                   row_blk_sizes, col_blk_sizes, &
                                   row_blk_sizes2enum, enum2row_blk_sizes, &
                                   col_blk_sizes2enum, enum2col_blk_sizes, &
                                   noresize)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      TYPE(dbcsr_data_obj), INTENT(INOUT)                :: trs_stackbuf
      INTEGER, DIMENSION(:), INTENT(IN), POINTER :: row_blk_sizes, col_blk_sizes, &
         row_blk_sizes2enum, enum2row_blk_sizes, col_blk_sizes2enum, enum2col_blk_sizes
      LOGICAL, INTENT(IN), OPTIONAL                      :: noresize

      CHARACTER(len=*), PARAMETER :: routineN = 'acc_transpose_blocks', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk_p, handle, handle1, i, m, mi, &
                                                            mi_max, n, nblks, ni, ni_max, offset, x, &
                                                            row, col
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: counters, filled, offsets, tmp_stack
      INTEGER, DIMENSION(:), POINTER                     :: blk_index, trs_stack
      LOGICAL                                            :: my_noresize

      CALL timeset(routineN, handle)

      NULLIFY (trs_stack)

      IF (.NOT. matrix%list_indexing) &
         CPABORT("build_trs_stack: only list_indexing supported.")
      IF (.NOT. matrix%local_indexing) &
         CPABORT("build_trs_stack: only local_indexing supported.")
      IF (trs_stackbuf%d%data_type /= dbcsr_type_int_4) &
         CPABORT("build_trs_stack: stac_buf has wrong datatype")
      blk_index => matrix%coo_l
      nblks = matrix%nblks

      ! make sure buffer from previous cannon-tick was uploaded
      CALL timeset(routineN//"_sync", handle1)
      CALL acc_event_synchronize(trs_stackbuf%d%acc_ready)
      CALL timestop(handle1)

      CALL timeset(routineN//"_ensure", handle1)
      my_noresize = .FALSE.
      IF (PRESENT(noresize)) my_noresize = noresize
      IF (my_noresize) THEN
         IF (dbcsr_data_get_size(trs_stackbuf) .LT. nblks) &
            CPABORT("build_trs_stack: trs_stackbuf undersized")
      ELSE
         CALL dbcsr_data_ensure_size(trs_stackbuf, data_size=nblks, nocopy=.TRUE.)
      ENDIF
      CALL dbcsr_data_set_size_referenced(trs_stackbuf, nblks)
      trs_stack => trs_stackbuf%d%i4
      CALL timestop(handle1)

      mi_max = SIZE(enum2row_blk_sizes); ni_max = SIZE(enum2col_blk_sizes)
      ALLOCATE (counters(mi_max, ni_max), offsets(mi_max, ni_max))
      counters(:, :) = 0; offsets(:, :) = 0

      CALL timeset(routineN//"_comp", handle1)

      ! Simplified algorithm for single size blocks
      IF (mi_max .EQ. 1 .AND. ni_max .EQ. 1) THEN
         DO i = 1, nblks
            blk_p = blk_index(3*(i-1)+3)
            IF (blk_p == 0) CYCLE
            counters(1, 1) = counters(1, 1)+1
            trs_stack(counters(1, 1)) = blk_p-1
         ENDDO
      ELSE
         ALLOCATE (tmp_stack(3, nblks))
         ! collect block addresses and dimensions in a temporary stack
         ! while doing so, also count number of blocks per block-dimensions
         DO i = 1, nblks
            row = blk_index(3*(i-1)+1)
            col = blk_index(3*(i-1)+2)
            blk_p = blk_index(3*(i-1)+3)
            IF (blk_p == 0) CYCLE
            m = row_blk_sizes(row)
            n = col_blk_sizes(col)
            mi = row_blk_sizes2enum(m)
            ni = col_blk_sizes2enum(n)
            tmp_stack(1, i) = mi
            tmp_stack(2, i) = ni
            tmp_stack(3, i) = blk_p-1
            counters(mi, ni) = counters(mi, ni)+1
         ENDDO
         ! calculate offsets for first element of each sub-stack
         offset = 0
         DO mi = 1, mi_max
            DO ni = 1, ni_max
               offsets(mi, ni) = offset
               offset = offset+counters(mi, ni)
            ENDDO
         ENDDO
         ! write all sub-stacks into the host-pinned buffer
         ALLOCATE (filled(mi_max, ni_max))
         filled(:, :) = 0
         DO i = 1, nblks
            mi = tmp_stack(1, i)
            ni = tmp_stack(2, i)
            blk_p = tmp_stack(3, i)
            x = offsets(mi, ni)+filled(mi, ni)+1
            trs_stack(x) = blk_p
            filled(mi, ni) = filled(mi, ni)+1
         ENDDO
         !sanity check
         DO ni = 1, ni_max
            DO mi = 1, mi_max
               IF (filled(mi, ni) /= counters(mi, ni)) &
                  CPABORT("acc_transpose_blocks: bug")
            END DO
         END DO
      ENDIF
      CALL timestop(handle1)

      CALL timeset(routineN//"_sync", handle1)
      !transfer all stacks
      CALL dbcsr_data_host2dev(trs_stackbuf)
      ! make sure block-buffer is uploaded befor running the kernels
      CALL acc_stream_wait_event(trs_stackbuf%d%memory_type%acc_stream, matrix%data_area%d%acc_ready)
      CALL timestop(handle1)

      CALL timeset(routineN//"_kernels", handle1)
      ! launch kernels
      DO ni = 1, ni_max
         DO mi = 1, mi_max
            IF (counters(mi, ni) > 0) THEN
               m = enum2row_blk_sizes(mi)
               n = enum2col_blk_sizes(ni)
               CALL dbcsr_acc_transpose( &
                  trs_stack=trs_stackbuf%d%acc_devmem, &
                  offset=offsets(mi, ni), &
                  nblks=counters(mi, ni), &
                  datatype=matrix%data_type, &
                  buffer=matrix%data_area%d%acc_devmem, &
                  m=m, n=n, &
                  stream=trs_stackbuf%d%memory_type%acc_stream)
            END IF
         ENDDO
      ENDDO
      CALL timestop(handle1)

      CALL timeset(routineN//"_sync", handle1)
      ! make sure block-buffer are not used until transpose kernels finnished
      CALL acc_event_record(trs_stackbuf%d%acc_ready, trs_stackbuf%d%memory_type%acc_stream)
      CALL acc_stream_wait_event(matrix%data_area%d%memory_type%acc_stream, trs_stackbuf%d%acc_ready)
      CALL acc_event_record(matrix%data_area%d%acc_ready, matrix%data_area%d%memory_type%acc_stream)
      CALL timestop(handle1)

      CALL timestop(handle)
   END SUBROUTINE acc_transpose_blocks

! **************************************************************************************************
!> \brief Guess the size of the product matrix from the A and B sparsities
!> \param matrix_left ...
!> \param matrix_right ...
!> \param product_matrix ...
!> \param left_data_size ...
!> \param right_data_size ...
!> \param left_col_nimages ...
!> \param right_row_nimages ...
!> \param nthreads ...
!> \retval size_guess ...
!> \author Alfio Lazzaro
! **************************************************************************************************
   FUNCTION product_matrix_size_guess(matrix_left, matrix_right, product_matrix, &
                                      left_data_size, right_data_size, &
                                      left_col_nimages, right_row_nimages, &
                                      nthreads) RESULT(size_guess)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix_left, matrix_right, product_matrix
      INTEGER, INTENT(IN)                                :: left_data_size, right_data_size, &
                                                            left_col_nimages, right_row_nimages, &
                                                            nthreads
      INTEGER                                            :: size_guess

      INTEGER(KIND=int_8)                                :: size8
      REAL(kind=real_8)                                  :: factor, fill_guess, left_fill, right_fill

      ! First we calculate the sparsities
      size8 = INT(dbcsr_nfullrows_local(matrix_left), KIND=int_8)* &
              INT(dbcsr_nfullcols_local(matrix_left), KIND=int_8)
      size8 = MAX(1_int_8, size8)
      left_fill = (REAL(left_data_size, KIND=real_8)*REAL(left_col_nimages, KIND=real_8))/REAL(size8, KIND=real_8)
      size8 = INT(dbcsr_nfullrows_local(matrix_right), KIND=int_8)* &
              INT(dbcsr_nfullcols_local(matrix_right), KIND=int_8)
      size8 = MAX(1_int_8, size8)
      right_fill = (REAL(right_data_size, KIND=real_8)*REAL(right_row_nimages, KIND=real_8))/REAL(size8, KIND=real_8)
      size8 = INT(dbcsr_nfullrows_local(product_matrix), KIND=int_8)* &
              INT(dbcsr_nfullcols_local(product_matrix), KIND=int_8)
      size8 = MAX(1_int_8, size8)
!     factor = 7.0 ! Old guess
      factor = 2.4 ! New guess
      fill_guess = factor*MAX(left_fill, right_fill)
      fill_guess = MIN(1.0_real_8, MAX(0.0_real_8, fill_guess))
      IF (nthreads .GT. 1) THEN
         fill_guess = fill_guess*3.0_real_8/REAL(2*nthreads, KIND=real_8)
      ENDIF
      size_guess = INT(REAL(size8, KIND=real_8)*fill_guess, KIND=int_4)
   END FUNCTION product_matrix_size_guess

! **************************************************************************************************
!> \brief Calculates per-block norms.
!>
!> Rewritten to be very low-level.
!> \param[in,out] matrix     DBCSR matrix for which to calculate norms
!> \param[out] norms         Block norms
! **************************************************************************************************
   SUBROUTINE calculate_norms(matrix, norms, row_blk_sizes, col_blk_sizes)
      TYPE(dbcsr_type), INTENT(IN)                       :: matrix
      REAL(kind=sp), DIMENSION(:), INTENT(OUT)           :: norms
      INTEGER, DIMENSION(:), POINTER, INTENT(IN)         :: row_blk_sizes, col_blk_sizes

      CHARACTER(len=*), PARAMETER :: routineN = 'calculate_norms', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: data_type, handle, nblks

!   ---------------------------------------------------------------------------

      CALL timeset(routineN, handle)
      ! Checks for validity
      IF (.NOT. dbcsr_valid_index(matrix)) &
         CPABORT("The matrix must be valid.")
      data_type = dbcsr_get_data_type(matrix)
      nblks = matrix%nblks
      SELECT CASE (data_type)
      CASE (dbcsr_type_real_4)
         CALL calc_norms_s(norms, nblks, matrix%coo_l, &
                           row_blk_sizes, col_blk_sizes, &
                           dbcsr_get_data_p_s(matrix%data_area))
      CASE (dbcsr_type_real_8)
         CALL calc_norms_d(norms, nblks, matrix%coo_l, &
                           row_blk_sizes, col_blk_sizes, &
                           dbcsr_get_data_p_d(matrix%data_area))
      CASE (dbcsr_type_complex_4)
         CALL calc_norms_c(norms, nblks, matrix%coo_l, &
                           row_blk_sizes, col_blk_sizes, &
                           dbcsr_get_data_p_c(matrix%data_area))
      CASE (dbcsr_type_complex_8)
         CALL calc_norms_z(norms, nblks, matrix%coo_l, &
                           row_blk_sizes, col_blk_sizes, &
                           dbcsr_get_data_p_z(matrix%data_area))
      CASE DEFAULT
         CPABORT("Invalid data type.")
      END SELECT
      !
      CALL timestop(handle)
   END SUBROUTINE calculate_norms

! **************************************************************************************************
!> \brief Gathers the local elements from all data (full_data)
!>
!> \param[in] full_data       All elements
!> \param[in] nle             Number of local elements
!> \param[in] local_elements  List of local elements
!> \param[out] local_data     Local elements obtained from all elements
! **************************************************************************************************
   PURE SUBROUTINE local_filter(full_data, nle, local_elements, local_data)
      INTEGER, DIMENSION(:), INTENT(IN)                  :: full_data
      INTEGER, INTENT(IN)                                :: nle
      INTEGER, DIMENSION(1:nle), INTENT(IN)              :: local_elements
      INTEGER, DIMENSION(1:nle), INTENT(OUT)             :: local_data

      INTEGER                                            :: l

      DO l = 1, nle
         local_data(l) = full_data(local_elements(l))
      END DO
   END SUBROUTINE local_filter

#:include "dbcsr_mm_common.f90"

END MODULE dbcsr_mm_common
