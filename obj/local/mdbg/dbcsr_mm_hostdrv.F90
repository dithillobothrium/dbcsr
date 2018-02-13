# 1 "/data/isivkov/libdbcsr_svn18247/src/dbcsr/mm/dbcsr_mm_hostdrv.F"
!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Stacks of small matrix multiplications
!> \author  Urban Borstnik
!> \date    2011-09-26
!> \version 0.9
!>
!> <b>Modification history:</b>
!  - 2011-09-26 Split dbcsr_internal_operations
! **************************************************************************************************
MODULE dbcsr_mm_hostdrv
   USE dbcsr_config,                    ONLY: dbcsr_cfg,&
                                              has_acc,&
                                              mm_driver_blas,&
                                              mm_driver_matmul,&
                                              mm_driver_smm,&
                                              mm_driver_xsmm
   USE dbcsr_data_methods,              ONLY: dbcsr_data_get_size
   USE dbcsr_mm_types,                  ONLY: dbcsr_ps_width,&
                                              p_a_first,&
                                              p_b_first,&
                                              p_c_first,&
                                              p_k,&
                                              p_m,&
                                              p_n,&
                                              stack_descriptor_type
   USE dbcsr_types,                     ONLY: dbcsr_data_obj,&
                                              dbcsr_type,&
                                              dbcsr_type_complex_4,&
                                              dbcsr_type_complex_8,&
                                              dbcsr_type_real_4,&
                                              dbcsr_type_real_8,&
                                              dbcsr_work_type
   USE kinds,                           ONLY: dp,&
                                              int_8,&
                                              real_4,&
                                              real_8,&
                                              sp
#include "../../base/base_uses.f90"

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_mm_hostdrv'

   CHARACTER(len=*), PARAMETER, PRIVATE :: int_print = "(10(1X,I7))"

   PUBLIC :: dbcsr_mm_hostdrv_lib_init, dbcsr_mm_hostdrv_lib_finalize
   PUBLIC :: dbcsr_mm_hostdrv_process
   PUBLIC :: dbcsr_mm_hostdrv_type
   PUBLIC :: dbcsr_mm_hostdrv_init

   LOGICAL, PARAMETER :: debug_mod = .FALSE.
   LOGICAL, PARAMETER :: careful_mod = .FALSE.

   TYPE dbcsr_mm_hostdrv_type
      TYPE(dbcsr_data_obj)          :: data_area
   END TYPE dbcsr_mm_hostdrv_type

CONTAINS

! **************************************************************************************************
!> \brief Initialize the library
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbcsr_mm_hostdrv_lib_init()
#if defined(__LIBXSMM)
    USE libxsmm,                           ONLY:  libxsmm_init
!$OMP     MASTER
      CALL libxsmm_init()
!$OMP     END MASTER
!$OMP     BARRIER
#endif
   END SUBROUTINE dbcsr_mm_hostdrv_lib_init

! **************************************************************************************************
!> \brief Finalize the library
!> \author Ole Schuett
! **************************************************************************************************
   SUBROUTINE dbcsr_mm_hostdrv_lib_finalize()
#if defined(__LIBXSMM)
    USE libxsmm,                           ONLY:  libxsmm_finalize
!$OMP     BARRIER
!$OMP     MASTER
      CALL libxsmm_finalize()
!$OMP     END MASTER
#endif
   END SUBROUTINE dbcsr_mm_hostdrv_lib_finalize

! **************************************************************************************************
!> \brief Initialize the library
!> \param this ...
!> \param product_wm ...
! **************************************************************************************************
   SUBROUTINE dbcsr_mm_hostdrv_init(this, product_wm)
      TYPE(dbcsr_mm_hostdrv_type), INTENT(INOUT)         :: this
      TYPE(dbcsr_work_type), POINTER                     :: product_wm

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_hostdrv_init', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)

      this%data_area = product_wm%data_area

      CALL timestop(handle)

   END SUBROUTINE dbcsr_mm_hostdrv_init

! **************************************************************************************************
!> \brief Calls the various drivers that process the stack.
!>
!> \param this ...
!> \param[in] left Left-matrix data
!> \param[in] right Right-matrix data
!> \param[in] params           Stack of GEMM parameters
!> \param stack_size ...
!> \param stack_descr ...
!> \param success ...
!> \param used_smm ...
! **************************************************************************************************
   SUBROUTINE dbcsr_mm_hostdrv_process(this, left, right, params, stack_size, &
                                       stack_descr, success, used_smm)
      TYPE(dbcsr_mm_hostdrv_type), INTENT(INOUT)         :: this
      TYPE(dbcsr_type), INTENT(IN)                       :: left, right
      INTEGER, INTENT(IN)                                :: stack_size
      INTEGER, DIMENSION(1:dbcsr_ps_width, stack_size), &
         INTENT(INOUT)                                   :: params
      TYPE(stack_descriptor_type), INTENT(IN)            :: stack_descr
      LOGICAL, INTENT(OUT)                               :: success, used_smm

      CHARACTER(len=*), PARAMETER :: routineN = 'dbcsr_mm_hostdrv_process', &
         routineP = moduleN//':'//routineN
      LOGICAL, PARAMETER                                 :: careful = careful_mod, dbg = .FALSE.

      INTEGER                                            :: error_handle, sp
      REAL(KIND=dp)                                      :: rnd

      IF (has_acc) & !for cpu-only runs this is called too often
         CALL timeset(routineN, error_handle)

      success = .TRUE. !host driver never failes...hopefully
      used_smm = .FALSE.

      IF (dbg) THEN
         CALL RANDOM_NUMBER(rnd)
         IF (rnd < 0.01_dp) THEN
            WRITE (*, *) routineN//" Stack size", stack_size, dbcsr_ps_width
            CALL print_gemm_parameters(params(:, 1:stack_size))
         ENDIF
      ENDIF

      ! Verify stack consistency.  Only the upper bound is verified.
      IF (careful) THEN
         DO sp = 1, stack_size
            IF (params(p_a_first, sp)+params(p_m, sp)*params(p_k, sp)-1 > dbcsr_data_get_size(left%data_area)) &
               CPABORT("A data out of bounds.")
            IF (params(p_b_first, sp)+params(p_k, sp)*params(p_n, sp)-1 > dbcsr_data_get_size(right%data_area)) &
               CPABORT("B data out of bounds.")
            IF (params(p_c_first, sp)+params(p_m, sp)*params(p_n, sp)-1 > dbcsr_data_get_size(this%data_area)) &
               CPABORT("C data out of bounds.")
         ENDDO
      ENDIF

      SELECT CASE (dbcsr_cfg%mm_driver)
      CASE (mm_driver_matmul)
         SELECT CASE (this%data_area%d%data_type)
         CASE (dbcsr_type_real_4)
            CALL internal_process_mm_stack_s(params, &
                                             stack_size, &
                                             left%data_area%d%r_sp, right%data_area%d%r_sp, this%data_area%d%r_sp)
         CASE (dbcsr_type_real_8)
            CALL internal_process_mm_stack_d(params, &
                                             stack_size, &
                                             left%data_area%d%r_dp, right%data_area%d%r_dp, this%data_area%d%r_dp)
         CASE (dbcsr_type_complex_4)
            CALL internal_process_mm_stack_c(params, &
                                             stack_size, &
                                             left%data_area%d%c_sp, right%data_area%d%c_sp, this%data_area%d%c_sp)
         CASE (dbcsr_type_complex_8)
            CALL internal_process_mm_stack_z(params, &
                                             stack_size, &
                                             left%data_area%d%c_dp, right%data_area%d%c_dp, this%data_area%d%c_dp)
         CASE default
            CPABORT("Invalid data type")
         END SELECT
      CASE (mm_driver_smm)
         SELECT CASE (this%data_area%d%data_type)
         CASE (dbcsr_type_real_4)
            CALL smm_process_mm_stack_s(stack_descr, params, &
                                        stack_size, &
                                        left%data_area%d%r_sp, right%data_area%d%r_sp, this%data_area%d%r_sp, used_smm)
         CASE (dbcsr_type_real_8)
            CALL smm_process_mm_stack_d(stack_descr, params, &
                                        stack_size, &
                                        left%data_area%d%r_dp, right%data_area%d%r_dp, this%data_area%d%r_dp, used_smm)
         CASE (dbcsr_type_complex_4)
            CALL smm_process_mm_stack_c(stack_descr, params, &
                                        stack_size, &
                                        left%data_area%d%c_sp, right%data_area%d%c_sp, this%data_area%d%c_sp, used_smm)
         CASE (dbcsr_type_complex_8)
            CALL smm_process_mm_stack_z(stack_descr, params, &
                                        stack_size, &
                                        left%data_area%d%c_dp, right%data_area%d%c_dp, this%data_area%d%c_dp, used_smm)
         CASE default
            CPABORT("Invalid data type")
         END SELECT

      CASE (mm_driver_xsmm)
         SELECT CASE (this%data_area%d%data_type)
         CASE (dbcsr_type_real_4)
            CALL xsmm_process_mm_stack_s(stack_descr, params, stack_size, &
                                         left%data_area%d%r_sp, right%data_area%d%r_sp, this%data_area%d%r_sp, used_smm)
         CASE (dbcsr_type_real_8)
            CALL xsmm_process_mm_stack_d(stack_descr, params, stack_size, &
                                         left%data_area%d%r_dp, right%data_area%d%r_dp, this%data_area%d%r_dp, used_smm)
         CASE (dbcsr_type_complex_4)
            CALL xsmm_process_mm_stack_c(stack_descr, params, stack_size, &
                                         left%data_area%d%c_sp, right%data_area%d%c_sp, this%data_area%d%c_sp, used_smm)
         CASE (dbcsr_type_complex_8)
            CALL xsmm_process_mm_stack_z(stack_descr, params, stack_size, &
                                         left%data_area%d%c_dp, right%data_area%d%c_dp, this%data_area%d%c_dp, used_smm)
         CASE default
            CPABORT("Invalid data type")
         END SELECT

      CASE (mm_driver_blas)
         SELECT CASE (this%data_area%d%data_type)
         CASE (dbcsr_type_real_4)
            CALL blas_process_mm_stack_s(params, &
                                         stack_size, &
                                         left%data_area%d%r_sp, right%data_area%d%r_sp, this%data_area%d%r_sp)
         CASE (dbcsr_type_real_8)
            CALL blas_process_mm_stack_d(params, &
                                         stack_size, &
                                         left%data_area%d%r_dp, right%data_area%d%r_dp, this%data_area%d%r_dp)
         CASE (dbcsr_type_complex_4)
            CALL blas_process_mm_stack_c(params, &
                                         stack_size, &
                                         left%data_area%d%c_sp, right%data_area%d%c_sp, this%data_area%d%c_sp)
         CASE (dbcsr_type_complex_8)
            CALL blas_process_mm_stack_z(params, &
                                         stack_size, &
                                         left%data_area%d%c_dp, right%data_area%d%c_dp, this%data_area%d%c_dp)
         CASE default
            CPABORT("Invalid data type")
         END SELECT
      CASE default
         CPABORT("Invalid multiplication driver")
      END SELECT

      IF (has_acc) & !for cpu-only runs this is called too often
         CALL timestop(error_handle)

   END SUBROUTINE dbcsr_mm_hostdrv_process

! **************************************************************************************************
!> \brief Helper-routine used by dbcsr_mm_hostdrv_process to print debug info.
!> \param params ...
! **************************************************************************************************
   SUBROUTINE print_gemm_parameters(params)
      INTEGER, DIMENSION(:, :), INTENT(in)               :: params

      INTEGER                                            :: sp

      DO sp = 1, SIZE(params, 2)
         WRITE (*, '(1X,A,1X,I7,":",3(1X,I4,","),".",3(1X,I12,","))') &
            "GEMM PARAMETERS", &
            sp, &
            params(p_m, sp), &
            params(p_k, sp), &
            params(p_n, sp), &
            params(p_a_first, sp), &
            params(p_b_first, sp), &
            params(p_c_first, sp)
      ENDDO
   END SUBROUTINE print_gemm_parameters

# 1 "/data/isivkov/libdbcsr_svn18247/src/dbcsr/mm/dbcsr_mm_hostdrv.f90" 1
!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

# 1 "/data/isivkov/libdbcsr_svn18247/src/dbcsr/mm/../data/dbcsr.fypp" 1
# 8 "/data/isivkov/libdbcsr_svn18247/src/dbcsr/mm/../data/dbcsr.fypp"
# 39 "/data/isivkov/libdbcsr_svn18247/src/dbcsr/mm/../data/dbcsr.fypp"
# 7 "/data/isivkov/libdbcsr_svn18247/src/dbcsr/mm/dbcsr_mm_hostdrv.f90" 2
# 8 "/data/isivkov/libdbcsr_svn18247/src/dbcsr/mm/dbcsr_mm_hostdrv.f90"
! **************************************************************************************************
!> \brief Processes MM stack and issues BLAS xGEMM calls
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! **************************************************************************************************
  SUBROUTINE blas_process_mm_stack_d(params,&
       stack_size,&
       a_data, b_data, c_data)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    REAL(kind=real_8), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    REAL(kind=real_8), DIMENSION(*), INTENT(INOUT)      :: c_data

    CHARACTER(len=*), PARAMETER :: routineN = 'blas_process_mm_stack_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL DGEMM('N',&
            'N',&
            params(p_m,sp), params(p_n,sp),& !m, n
            params(p_k,sp),& ! k
            1.0_real_8,& ! alpha
            a_data(params(p_a_first,sp)),& ! A
            params(p_m,sp),& !lda
            b_data(params(p_b_first,sp)),& ! B
            params(p_k,sp),& !ldb
            1.0_real_8,& ! beta
            c_data(params(p_c_first,sp)), params(p_m,sp))
    ENDDO
  END SUBROUTINE blas_process_mm_stack_d

! **************************************************************************************************
!> \brief Processes MM stack and issues internal MM calls.
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! **************************************************************************************************
  SUBROUTINE internal_process_mm_stack_d(params, stack_size,&
       a_data, b_data, c_data)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    REAL(kind=real_8), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    REAL(kind=real_8), DIMENSION(*), INTENT(INOUT)      :: c_data

    CHARACTER(len=*), PARAMETER :: routineN = 'internal_process_mm_stack_d', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL internal_mm_d_nn(&
            params(p_m,sp),&
            params(p_n,sp),&
            params(p_k,sp),&
            a_data(params(p_a_first,sp)),&
            b_data(params(p_b_first,sp)),&
            c_data(params(p_c_first,sp)))
    ENDDO
  END SUBROUTINE internal_process_mm_stack_d


! **************************************************************************************************
!> \brief Processes MM stack and issues SMM library calls
!>
!> \param stack_descr ...
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! **************************************************************************************************
  SUBROUTINE smm_process_mm_stack_d(stack_descr, params,&
       stack_size,&
       a_data, b_data, c_data, used_smm)
    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    REAL(kind=real_8), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    REAL(kind=real_8), DIMENSION(*), INTENT(INOUT)      :: c_data
    LOGICAL, INTENT(OUT)                      :: used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'smm_process_mm_stack_d', &
      routineP = moduleN//':'//routineN

#if defined(__HAS_smm_dnn)

   INTEGER                                   :: sp

   ! TODO we have no way of knowing which calls to libsmm actually resolve to BLAS
   ! Fixing this requires an interface change to libsmm.
   used_smm=.TRUE.

#if defined(__HAS_smm_vec)
    IF(stack_descr%defined_mnk) THEN
       CALL smm_vec_dnn (stack_descr%m, stack_descr%n, stack_descr%k, &
         a_data, b_data, c_data, stack_size, &
         dbcsr_ps_width, params, p_a_first, p_b_first, p_c_first)
       RETURN
    ENDIF
#endif

    DO sp = 1, stack_size
       CALL smm_dnn(&
            params(p_m,sp),&
            params(p_n,sp),&
            params(p_k,sp),&
            a_data(params(p_a_first,sp)),&
            b_data(params(p_b_first,sp)),&
            c_data(params(p_c_first,sp)))
    ENDDO

#else
    ! We do not want to abort here, fall back to BLAS.
    used_smm=.FALSE.
    CALL blas_process_mm_stack_d(params, stack_size,a_data, b_data, c_data)
#endif

    MARK_USED(stack_descr)
  END SUBROUTINE smm_process_mm_stack_d


! **************************************************************************************************
!> \brief Processes MM stack and issues libxsmm calls
!>
!> \param stack_descr ...
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
!> \param[out] used_smm        Flag to signal if an efficient kernel was used
!> \author Ole Schuett
! **************************************************************************************************
  SUBROUTINE xsmm_process_mm_stack_d(stack_descr, params,&
       stack_size, a_data, b_data, c_data, used_smm)

#if defined(__LIBXSMM) && 1
    ! Caution: This dependency is ignored by makedep.py, because libxsmm.F is kinda empty.
    USE libxsmm,                           ONLY: libxsmm_function  => libxsmm_dmmfunction,&
                                                 libxsmm_dispatch  => libxsmm_dmmdispatch,&
                                                 libxsmm_available => libxsmm_dmmavailable,&
                                                 libxsmm_call      => libxsmm_dmmcall,&
                                                 libxsmm_gemm      => libxsmm_dgemm,&
                                                 LIBXSMM_PREFETCH_NONE,&
                                                 LIBXSMM_PREFETCH,&
                                                 LIBXSMM_ROW_MAJOR,&
                                                 LIBXSMM_COL_MAJOR,&
                                                 LIBXSMM_MAX_MNK,&
                                                 LIBXSMM_FLAGS

    INTEGER, PARAMETER :: LIBXSMM_DEFAULT_PREFETCH = LIBXSMM_PREFETCH
    INTEGER, PARAMETER :: LIBXSMM_DEFAULT_FLAGS = LIBXSMM_FLAGS
#endif

    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    REAL(kind=real_8), DIMENSION(*), TARGET, INTENT(IN) :: a_data, b_data
    REAL(kind=real_8), DIMENSION(*), TARGET, &
      INTENT(INOUT)                           :: c_data
    LOGICAL, INTENT(OUT)                      :: used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'libxsmm_process_mm_stack_d', &
      routineP = moduleN//':'//routineN

#if defined(__LIBXSMM) && 1
    REAL(real_8), PARAMETER                  :: one = 1.0_real_8
    LOGICAL                                   :: processed
    INTEGER(int_8)                            :: threshold
    INTEGER                                   :: fa, fb, fc, m, n, k, pa, pb, pc, sp
    REAL(real_8), DIMENSION(:,:), POINTER    :: a_ptr, b_ptr, c_ptr
    TYPE(libxsmm_function)                    :: func

    processed = .FALSE.
    used_smm = .FALSE.

    CPASSERT(LIBXSMM_COL_MAJOR/=0 .AND. LIBXSMM_ROW_MAJOR==0)

    ! check whether the matrix stack is homogeneous or not
    IF (stack_descr%defined_mnk) THEN
       threshold = INT(stack_descr%m, int_8) * &
                   INT(stack_descr%n, int_8) * &
                   INT(stack_descr%k, int_8)

       ! check if matrices are too large for LIBXSMM (BLAS is likely more efficient)
       IF(threshold > LIBXSMM_MAX_MNK) THEN
          CALL blas_process_mm_stack_d(params, stack_size, a_data, b_data, c_data)
          processed = .TRUE.

       ELSE
          ! try to get a function pointer from libxsmm
          CALL libxsmm_dispatch(func, &
               m=stack_descr%m, n=stack_descr%n, k=stack_descr%k, alpha=one, beta=one, &
               flags=LIBXSMM_DEFAULT_FLAGS, prefetch=LIBXSMM_DEFAULT_PREFETCH)

          IF (libxsmm_available(func)) THEN
             ! load first stack entry
             CPASSERT(stack_size > 0)
             pa = params(p_a_first, 1)
             pb = params(p_b_first, 1)
             pc = params(p_c_first, 1)

             DO sp = 1, stack_size-1
                fa = pa; fb = pb; fc = pc
                ! prefetch next blocks
                pa = params(p_a_first, sp + 1)
                pb = params(p_b_first, sp + 1)
                pc = params(p_c_first, sp + 1)

                ! condition evaluates at compile-time (PARAMETERS)
                IF (LIBXSMM_DEFAULT_PREFETCH /= LIBXSMM_PREFETCH_NONE) THEN
                   CALL libxsmm_call(func, &
                        a=a_data(fa), b=b_data(fb), c=c_data(fc), &
                        ! provide locations of the next operand set
                        pa=a_data(pa), pb=b_data(pb), pc=c_data(pc))
                ELSE
                   CALL libxsmm_call(func, &
                        a=a_data(fa), b=b_data(fb), c=c_data(fc))
                ENDIF
             ENDDO

             ! handle last stack entry without out-of-bounds access
             fa = pa; fb = pb; fc = pc

             ! condition evaluates at compile-time (PARAMETERS)
             IF (LIBXSMM_DEFAULT_PREFETCH /= LIBXSMM_PREFETCH_NONE) THEN
                CALL libxsmm_call(func, &
                     a=a_data(fa), b=b_data(fb), c=c_data(fc), &
                     ! prefetch same blocks
                     pa=a_data(pa), pb=b_data(pb), pc=c_data(pc))
             ELSE
                CALL libxsmm_call(func, &
                     a=a_data(fa), b=b_data(fb), c=c_data(fc))
             ENDIF

             processed = .TRUE.
             used_smm  = .TRUE.
          ENDIF
       ENDIF
    ENDIF


    IF(.NOT.processed) THEN
       ! Dispatch interface was not used, call regular interface.
       ! Should only happen for inhomgenous stacks, then prefetching makes no sense.
       DO sp = 1, stack_size
          m = params(p_m,sp)
          n = params(p_n,sp)
          k = params(p_k,sp)
          fa = params(p_a_first,sp)
          fb = params(p_b_first,sp)
          fc = params(p_c_first,sp)
          ! somewhat expensive pointer remapping required by libxsmm interface
          a_ptr(1:m,1:k) => a_data(fa:fa+(m*k))
          b_ptr(1:k,1:n) => b_data(fb:fb+(k*n))
          c_ptr(1:m,1:n) => c_data(fc:fc+(m*n))
          CALL libxsmm_gemm(m=m, n=n, k=k, a=a_ptr, b=b_ptr, c=c_ptr, &
                            alpha=one, beta=one)
       ENDDO
    ENDIF

#else
    MARK_USED(stack_descr)
    ! We do not want to abort here, fall back to BLAS.
    CALL blas_process_mm_stack_d(params, stack_size,a_data, b_data, c_data)
    used_smm = .FALSE.
#endif

  END SUBROUTINE xsmm_process_mm_stack_d

! **************************************************************************************************
!> \brief ...
!> \param M ...
!> \param N ...
!> \param K ...
!> \param A ...
!> \param B ...
!> \param C ...
! **************************************************************************************************
  PURE SUBROUTINE internal_mm_d_nn(&
       M,N,K,A,B,C)
    INTEGER, INTENT(IN)                      :: M, N, K
    REAL(kind=real_8), INTENT(INOUT)                   :: C(M,N)
    REAL(kind=real_8), INTENT(IN)                      :: B(K,N)
    REAL(kind=real_8), INTENT(IN)                      :: A(M,K)
    C(:,:) = C(:,:) + MATMUL (A, B)
  END SUBROUTINE internal_mm_d_nn
# 8 "/data/isivkov/libdbcsr_svn18247/src/dbcsr/mm/dbcsr_mm_hostdrv.f90"
! **************************************************************************************************
!> \brief Processes MM stack and issues BLAS xGEMM calls
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! **************************************************************************************************
  SUBROUTINE blas_process_mm_stack_s(params,&
       stack_size,&
       a_data, b_data, c_data)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    REAL(kind=real_4), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    REAL(kind=real_4), DIMENSION(*), INTENT(INOUT)      :: c_data

    CHARACTER(len=*), PARAMETER :: routineN = 'blas_process_mm_stack_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL SGEMM('N',&
            'N',&
            params(p_m,sp), params(p_n,sp),& !m, n
            params(p_k,sp),& ! k
            1.0_real_4,& ! alpha
            a_data(params(p_a_first,sp)),& ! A
            params(p_m,sp),& !lda
            b_data(params(p_b_first,sp)),& ! B
            params(p_k,sp),& !ldb
            1.0_real_4,& ! beta
            c_data(params(p_c_first,sp)), params(p_m,sp))
    ENDDO
  END SUBROUTINE blas_process_mm_stack_s

! **************************************************************************************************
!> \brief Processes MM stack and issues internal MM calls.
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! **************************************************************************************************
  SUBROUTINE internal_process_mm_stack_s(params, stack_size,&
       a_data, b_data, c_data)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    REAL(kind=real_4), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    REAL(kind=real_4), DIMENSION(*), INTENT(INOUT)      :: c_data

    CHARACTER(len=*), PARAMETER :: routineN = 'internal_process_mm_stack_s', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL internal_mm_s_nn(&
            params(p_m,sp),&
            params(p_n,sp),&
            params(p_k,sp),&
            a_data(params(p_a_first,sp)),&
            b_data(params(p_b_first,sp)),&
            c_data(params(p_c_first,sp)))
    ENDDO
  END SUBROUTINE internal_process_mm_stack_s


! **************************************************************************************************
!> \brief Processes MM stack and issues SMM library calls
!>
!> \param stack_descr ...
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! **************************************************************************************************
  SUBROUTINE smm_process_mm_stack_s(stack_descr, params,&
       stack_size,&
       a_data, b_data, c_data, used_smm)
    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    REAL(kind=real_4), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    REAL(kind=real_4), DIMENSION(*), INTENT(INOUT)      :: c_data
    LOGICAL, INTENT(OUT)                      :: used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'smm_process_mm_stack_s', &
      routineP = moduleN//':'//routineN

#if defined(__HAS_smm_snn)

   INTEGER                                   :: sp

   ! TODO we have no way of knowing which calls to libsmm actually resolve to BLAS
   ! Fixing this requires an interface change to libsmm.
   used_smm=.TRUE.

#if defined(__HAS_smm_vec)
    IF(stack_descr%defined_mnk) THEN
       CALL smm_vec_snn (stack_descr%m, stack_descr%n, stack_descr%k, &
         a_data, b_data, c_data, stack_size, &
         dbcsr_ps_width, params, p_a_first, p_b_first, p_c_first)
       RETURN
    ENDIF
#endif

    DO sp = 1, stack_size
       CALL smm_snn(&
            params(p_m,sp),&
            params(p_n,sp),&
            params(p_k,sp),&
            a_data(params(p_a_first,sp)),&
            b_data(params(p_b_first,sp)),&
            c_data(params(p_c_first,sp)))
    ENDDO

#else
    ! We do not want to abort here, fall back to BLAS.
    used_smm=.FALSE.
    CALL blas_process_mm_stack_s(params, stack_size,a_data, b_data, c_data)
#endif

    MARK_USED(stack_descr)
  END SUBROUTINE smm_process_mm_stack_s


! **************************************************************************************************
!> \brief Processes MM stack and issues libxsmm calls
!>
!> \param stack_descr ...
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
!> \param[out] used_smm        Flag to signal if an efficient kernel was used
!> \author Ole Schuett
! **************************************************************************************************
  SUBROUTINE xsmm_process_mm_stack_s(stack_descr, params,&
       stack_size, a_data, b_data, c_data, used_smm)

#if defined(__LIBXSMM) && 1
    ! Caution: This dependency is ignored by makedep.py, because libxsmm.F is kinda empty.
    USE libxsmm,                           ONLY: libxsmm_function  => libxsmm_smmfunction,&
                                                 libxsmm_dispatch  => libxsmm_smmdispatch,&
                                                 libxsmm_available => libxsmm_smmavailable,&
                                                 libxsmm_call      => libxsmm_smmcall,&
                                                 libxsmm_gemm      => libxsmm_sgemm,&
                                                 LIBXSMM_PREFETCH_NONE,&
                                                 LIBXSMM_PREFETCH,&
                                                 LIBXSMM_ROW_MAJOR,&
                                                 LIBXSMM_COL_MAJOR,&
                                                 LIBXSMM_MAX_MNK,&
                                                 LIBXSMM_FLAGS

    INTEGER, PARAMETER :: LIBXSMM_DEFAULT_PREFETCH = LIBXSMM_PREFETCH
    INTEGER, PARAMETER :: LIBXSMM_DEFAULT_FLAGS = LIBXSMM_FLAGS
#endif

    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    REAL(kind=real_4), DIMENSION(*), TARGET, INTENT(IN) :: a_data, b_data
    REAL(kind=real_4), DIMENSION(*), TARGET, &
      INTENT(INOUT)                           :: c_data
    LOGICAL, INTENT(OUT)                      :: used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'libxsmm_process_mm_stack_s', &
      routineP = moduleN//':'//routineN

#if defined(__LIBXSMM) && 1
    REAL(real_4), PARAMETER                  :: one = 1.0_real_4
    LOGICAL                                   :: processed
    INTEGER(int_8)                            :: threshold
    INTEGER                                   :: fa, fb, fc, m, n, k, pa, pb, pc, sp
    REAL(real_4), DIMENSION(:,:), POINTER    :: a_ptr, b_ptr, c_ptr
    TYPE(libxsmm_function)                    :: func

    processed = .FALSE.
    used_smm = .FALSE.

    CPASSERT(LIBXSMM_COL_MAJOR/=0 .AND. LIBXSMM_ROW_MAJOR==0)

    ! check whether the matrix stack is homogeneous or not
    IF (stack_descr%defined_mnk) THEN
       threshold = INT(stack_descr%m, int_8) * &
                   INT(stack_descr%n, int_8) * &
                   INT(stack_descr%k, int_8)

       ! check if matrices are too large for LIBXSMM (BLAS is likely more efficient)
       IF(threshold > LIBXSMM_MAX_MNK) THEN
          CALL blas_process_mm_stack_s(params, stack_size, a_data, b_data, c_data)
          processed = .TRUE.

       ELSE
          ! try to get a function pointer from libxsmm
          CALL libxsmm_dispatch(func, &
               m=stack_descr%m, n=stack_descr%n, k=stack_descr%k, alpha=one, beta=one, &
               flags=LIBXSMM_DEFAULT_FLAGS, prefetch=LIBXSMM_DEFAULT_PREFETCH)

          IF (libxsmm_available(func)) THEN
             ! load first stack entry
             CPASSERT(stack_size > 0)
             pa = params(p_a_first, 1)
             pb = params(p_b_first, 1)
             pc = params(p_c_first, 1)

             DO sp = 1, stack_size-1
                fa = pa; fb = pb; fc = pc
                ! prefetch next blocks
                pa = params(p_a_first, sp + 1)
                pb = params(p_b_first, sp + 1)
                pc = params(p_c_first, sp + 1)

                ! condition evaluates at compile-time (PARAMETERS)
                IF (LIBXSMM_DEFAULT_PREFETCH /= LIBXSMM_PREFETCH_NONE) THEN
                   CALL libxsmm_call(func, &
                        a=a_data(fa), b=b_data(fb), c=c_data(fc), &
                        ! provide locations of the next operand set
                        pa=a_data(pa), pb=b_data(pb), pc=c_data(pc))
                ELSE
                   CALL libxsmm_call(func, &
                        a=a_data(fa), b=b_data(fb), c=c_data(fc))
                ENDIF
             ENDDO

             ! handle last stack entry without out-of-bounds access
             fa = pa; fb = pb; fc = pc

             ! condition evaluates at compile-time (PARAMETERS)
             IF (LIBXSMM_DEFAULT_PREFETCH /= LIBXSMM_PREFETCH_NONE) THEN
                CALL libxsmm_call(func, &
                     a=a_data(fa), b=b_data(fb), c=c_data(fc), &
                     ! prefetch same blocks
                     pa=a_data(pa), pb=b_data(pb), pc=c_data(pc))
             ELSE
                CALL libxsmm_call(func, &
                     a=a_data(fa), b=b_data(fb), c=c_data(fc))
             ENDIF

             processed = .TRUE.
             used_smm  = .TRUE.
          ENDIF
       ENDIF
    ENDIF


    IF(.NOT.processed) THEN
       ! Dispatch interface was not used, call regular interface.
       ! Should only happen for inhomgenous stacks, then prefetching makes no sense.
       DO sp = 1, stack_size
          m = params(p_m,sp)
          n = params(p_n,sp)
          k = params(p_k,sp)
          fa = params(p_a_first,sp)
          fb = params(p_b_first,sp)
          fc = params(p_c_first,sp)
          ! somewhat expensive pointer remapping required by libxsmm interface
          a_ptr(1:m,1:k) => a_data(fa:fa+(m*k))
          b_ptr(1:k,1:n) => b_data(fb:fb+(k*n))
          c_ptr(1:m,1:n) => c_data(fc:fc+(m*n))
          CALL libxsmm_gemm(m=m, n=n, k=k, a=a_ptr, b=b_ptr, c=c_ptr, &
                            alpha=one, beta=one)
       ENDDO
    ENDIF

#else
    MARK_USED(stack_descr)
    ! We do not want to abort here, fall back to BLAS.
    CALL blas_process_mm_stack_s(params, stack_size,a_data, b_data, c_data)
    used_smm = .FALSE.
#endif

  END SUBROUTINE xsmm_process_mm_stack_s

! **************************************************************************************************
!> \brief ...
!> \param M ...
!> \param N ...
!> \param K ...
!> \param A ...
!> \param B ...
!> \param C ...
! **************************************************************************************************
  PURE SUBROUTINE internal_mm_s_nn(&
       M,N,K,A,B,C)
    INTEGER, INTENT(IN)                      :: M, N, K
    REAL(kind=real_4), INTENT(INOUT)                   :: C(M,N)
    REAL(kind=real_4), INTENT(IN)                      :: B(K,N)
    REAL(kind=real_4), INTENT(IN)                      :: A(M,K)
    C(:,:) = C(:,:) + MATMUL (A, B)
  END SUBROUTINE internal_mm_s_nn
# 8 "/data/isivkov/libdbcsr_svn18247/src/dbcsr/mm/dbcsr_mm_hostdrv.f90"
! **************************************************************************************************
!> \brief Processes MM stack and issues BLAS xGEMM calls
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! **************************************************************************************************
  SUBROUTINE blas_process_mm_stack_z(params,&
       stack_size,&
       a_data, b_data, c_data)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(INOUT)      :: c_data

    CHARACTER(len=*), PARAMETER :: routineN = 'blas_process_mm_stack_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL ZGEMM('N',&
            'N',&
            params(p_m,sp), params(p_n,sp),& !m, n
            params(p_k,sp),& ! k
            CMPLX(1.0, 0.0, real_8),& ! alpha
            a_data(params(p_a_first,sp)),& ! A
            params(p_m,sp),& !lda
            b_data(params(p_b_first,sp)),& ! B
            params(p_k,sp),& !ldb
            CMPLX(1.0, 0.0, real_8),& ! beta
            c_data(params(p_c_first,sp)), params(p_m,sp))
    ENDDO
  END SUBROUTINE blas_process_mm_stack_z

! **************************************************************************************************
!> \brief Processes MM stack and issues internal MM calls.
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! **************************************************************************************************
  SUBROUTINE internal_process_mm_stack_z(params, stack_size,&
       a_data, b_data, c_data)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(INOUT)      :: c_data

    CHARACTER(len=*), PARAMETER :: routineN = 'internal_process_mm_stack_z', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL internal_mm_z_nn(&
            params(p_m,sp),&
            params(p_n,sp),&
            params(p_k,sp),&
            a_data(params(p_a_first,sp)),&
            b_data(params(p_b_first,sp)),&
            c_data(params(p_c_first,sp)))
    ENDDO
  END SUBROUTINE internal_process_mm_stack_z


! **************************************************************************************************
!> \brief Processes MM stack and issues SMM library calls
!>
!> \param stack_descr ...
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! **************************************************************************************************
  SUBROUTINE smm_process_mm_stack_z(stack_descr, params,&
       stack_size,&
       a_data, b_data, c_data, used_smm)
    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    COMPLEX(kind=real_8), DIMENSION(*), INTENT(INOUT)      :: c_data
    LOGICAL, INTENT(OUT)                      :: used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'smm_process_mm_stack_z', &
      routineP = moduleN//':'//routineN

#if defined(__HAS_smm_znn)

   INTEGER                                   :: sp

   ! TODO we have no way of knowing which calls to libsmm actually resolve to BLAS
   ! Fixing this requires an interface change to libsmm.
   used_smm=.TRUE.

#if defined(__HAS_smm_vec)
    IF(stack_descr%defined_mnk) THEN
       CALL smm_vec_znn (stack_descr%m, stack_descr%n, stack_descr%k, &
         a_data, b_data, c_data, stack_size, &
         dbcsr_ps_width, params, p_a_first, p_b_first, p_c_first)
       RETURN
    ENDIF
#endif

    DO sp = 1, stack_size
       CALL smm_znn(&
            params(p_m,sp),&
            params(p_n,sp),&
            params(p_k,sp),&
            a_data(params(p_a_first,sp)),&
            b_data(params(p_b_first,sp)),&
            c_data(params(p_c_first,sp)))
    ENDDO

#else
    ! We do not want to abort here, fall back to BLAS.
    used_smm=.FALSE.
    CALL blas_process_mm_stack_z(params, stack_size,a_data, b_data, c_data)
#endif

    MARK_USED(stack_descr)
  END SUBROUTINE smm_process_mm_stack_z


! **************************************************************************************************
!> \brief Processes MM stack and issues libxsmm calls
!>
!> \param stack_descr ...
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
!> \param[out] used_smm        Flag to signal if an efficient kernel was used
!> \author Ole Schuett
! **************************************************************************************************
  SUBROUTINE xsmm_process_mm_stack_z(stack_descr, params,&
       stack_size, a_data, b_data, c_data, used_smm)

#if defined(__LIBXSMM) && 0
    ! Caution: This dependency is ignored by makedep.py, because libxsmm.F is kinda empty.
    USE libxsmm,                           ONLY: libxsmm_function  => libxsmm_zmmfunction,&
                                                 libxsmm_dispatch  => libxsmm_zmmdispatch,&
                                                 libxsmm_available => libxsmm_zmmavailable,&
                                                 libxsmm_call      => libxsmm_zmmcall,&
                                                 libxsmm_gemm      => libxsmm_zgemm,&
                                                 LIBXSMM_PREFETCH_NONE,&
                                                 LIBXSMM_PREFETCH,&
                                                 LIBXSMM_ROW_MAJOR,&
                                                 LIBXSMM_COL_MAJOR,&
                                                 LIBXSMM_MAX_MNK,&
                                                 LIBXSMM_FLAGS

    INTEGER, PARAMETER :: LIBXSMM_DEFAULT_PREFETCH = LIBXSMM_PREFETCH
    INTEGER, PARAMETER :: LIBXSMM_DEFAULT_FLAGS = LIBXSMM_FLAGS
#endif

    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_8), DIMENSION(*), TARGET, INTENT(IN) :: a_data, b_data
    COMPLEX(kind=real_8), DIMENSION(*), TARGET, &
      INTENT(INOUT)                           :: c_data
    LOGICAL, INTENT(OUT)                      :: used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'libxsmm_process_mm_stack_z', &
      routineP = moduleN//':'//routineN

#if defined(__LIBXSMM) && 0
    REAL(real_8), PARAMETER                  :: one = 1.0_real_8
    LOGICAL                                   :: processed
    INTEGER(int_8)                            :: threshold
    INTEGER                                   :: fa, fb, fc, m, n, k, pa, pb, pc, sp
    REAL(real_8), DIMENSION(:,:), POINTER    :: a_ptr, b_ptr, c_ptr
    TYPE(libxsmm_function)                    :: func

    processed = .FALSE.
    used_smm = .FALSE.

    CPASSERT(LIBXSMM_COL_MAJOR/=0 .AND. LIBXSMM_ROW_MAJOR==0)

    ! check whether the matrix stack is homogeneous or not
    IF (stack_descr%defined_mnk) THEN
       threshold = INT(stack_descr%m, int_8) * &
                   INT(stack_descr%n, int_8) * &
                   INT(stack_descr%k, int_8)

       ! check if matrices are too large for LIBXSMM (BLAS is likely more efficient)
       IF(threshold > LIBXSMM_MAX_MNK) THEN
          CALL blas_process_mm_stack_z(params, stack_size, a_data, b_data, c_data)
          processed = .TRUE.

       ELSE
          ! try to get a function pointer from libxsmm
          CALL libxsmm_dispatch(func, &
               m=stack_descr%m, n=stack_descr%n, k=stack_descr%k, alpha=one, beta=one, &
               flags=LIBXSMM_DEFAULT_FLAGS, prefetch=LIBXSMM_DEFAULT_PREFETCH)

          IF (libxsmm_available(func)) THEN
             ! load first stack entry
             CPASSERT(stack_size > 0)
             pa = params(p_a_first, 1)
             pb = params(p_b_first, 1)
             pc = params(p_c_first, 1)

             DO sp = 1, stack_size-1
                fa = pa; fb = pb; fc = pc
                ! prefetch next blocks
                pa = params(p_a_first, sp + 1)
                pb = params(p_b_first, sp + 1)
                pc = params(p_c_first, sp + 1)

                ! condition evaluates at compile-time (PARAMETERS)
                IF (LIBXSMM_DEFAULT_PREFETCH /= LIBXSMM_PREFETCH_NONE) THEN
                   CALL libxsmm_call(func, &
                        a=a_data(fa), b=b_data(fb), c=c_data(fc), &
                        ! provide locations of the next operand set
                        pa=a_data(pa), pb=b_data(pb), pc=c_data(pc))
                ELSE
                   CALL libxsmm_call(func, &
                        a=a_data(fa), b=b_data(fb), c=c_data(fc))
                ENDIF
             ENDDO

             ! handle last stack entry without out-of-bounds access
             fa = pa; fb = pb; fc = pc

             ! condition evaluates at compile-time (PARAMETERS)
             IF (LIBXSMM_DEFAULT_PREFETCH /= LIBXSMM_PREFETCH_NONE) THEN
                CALL libxsmm_call(func, &
                     a=a_data(fa), b=b_data(fb), c=c_data(fc), &
                     ! prefetch same blocks
                     pa=a_data(pa), pb=b_data(pb), pc=c_data(pc))
             ELSE
                CALL libxsmm_call(func, &
                     a=a_data(fa), b=b_data(fb), c=c_data(fc))
             ENDIF

             processed = .TRUE.
             used_smm  = .TRUE.
          ENDIF
       ENDIF
    ENDIF


    IF(.NOT.processed) THEN
       ! Dispatch interface was not used, call regular interface.
       ! Should only happen for inhomgenous stacks, then prefetching makes no sense.
       DO sp = 1, stack_size
          m = params(p_m,sp)
          n = params(p_n,sp)
          k = params(p_k,sp)
          fa = params(p_a_first,sp)
          fb = params(p_b_first,sp)
          fc = params(p_c_first,sp)
          ! somewhat expensive pointer remapping required by libxsmm interface
          a_ptr(1:m,1:k) => a_data(fa:fa+(m*k))
          b_ptr(1:k,1:n) => b_data(fb:fb+(k*n))
          c_ptr(1:m,1:n) => c_data(fc:fc+(m*n))
          CALL libxsmm_gemm(m=m, n=n, k=k, a=a_ptr, b=b_ptr, c=c_ptr, &
                            alpha=one, beta=one)
       ENDDO
    ENDIF

#else
    MARK_USED(stack_descr)
    ! We do not want to abort here, fall back to BLAS.
    CALL blas_process_mm_stack_z(params, stack_size,a_data, b_data, c_data)
    used_smm = .FALSE.
#endif

  END SUBROUTINE xsmm_process_mm_stack_z

! **************************************************************************************************
!> \brief ...
!> \param M ...
!> \param N ...
!> \param K ...
!> \param A ...
!> \param B ...
!> \param C ...
! **************************************************************************************************
  PURE SUBROUTINE internal_mm_z_nn(&
       M,N,K,A,B,C)
    INTEGER, INTENT(IN)                      :: M, N, K
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: C(M,N)
    COMPLEX(kind=real_8), INTENT(IN)                      :: B(K,N)
    COMPLEX(kind=real_8), INTENT(IN)                      :: A(M,K)
    C(:,:) = C(:,:) + MATMUL (A, B)
  END SUBROUTINE internal_mm_z_nn
# 8 "/data/isivkov/libdbcsr_svn18247/src/dbcsr/mm/dbcsr_mm_hostdrv.f90"
! **************************************************************************************************
!> \brief Processes MM stack and issues BLAS xGEMM calls
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! **************************************************************************************************
  SUBROUTINE blas_process_mm_stack_c(params,&
       stack_size,&
       a_data, b_data, c_data)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_4), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    COMPLEX(kind=real_4), DIMENSION(*), INTENT(INOUT)      :: c_data

    CHARACTER(len=*), PARAMETER :: routineN = 'blas_process_mm_stack_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL CGEMM('N',&
            'N',&
            params(p_m,sp), params(p_n,sp),& !m, n
            params(p_k,sp),& ! k
            CMPLX(1.0, 0.0, real_4),& ! alpha
            a_data(params(p_a_first,sp)),& ! A
            params(p_m,sp),& !lda
            b_data(params(p_b_first,sp)),& ! B
            params(p_k,sp),& !ldb
            CMPLX(1.0, 0.0, real_4),& ! beta
            c_data(params(p_c_first,sp)), params(p_m,sp))
    ENDDO
  END SUBROUTINE blas_process_mm_stack_c

! **************************************************************************************************
!> \brief Processes MM stack and issues internal MM calls.
!>
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! **************************************************************************************************
  SUBROUTINE internal_process_mm_stack_c(params, stack_size,&
       a_data, b_data, c_data)
    INTEGER, INTENT(IN)                       :: stack_size
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_4), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    COMPLEX(kind=real_4), DIMENSION(*), INTENT(INOUT)      :: c_data

    CHARACTER(len=*), PARAMETER :: routineN = 'internal_process_mm_stack_c', &
      routineP = moduleN//':'//routineN

    INTEGER                                   :: sp

!   ---------------------------------------------------------------------------

    DO sp = 1, stack_size
       CALL internal_mm_c_nn(&
            params(p_m,sp),&
            params(p_n,sp),&
            params(p_k,sp),&
            a_data(params(p_a_first,sp)),&
            b_data(params(p_b_first,sp)),&
            c_data(params(p_c_first,sp)))
    ENDDO
  END SUBROUTINE internal_process_mm_stack_c


! **************************************************************************************************
!> \brief Processes MM stack and issues SMM library calls
!>
!> \param stack_descr ...
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
! **************************************************************************************************
  SUBROUTINE smm_process_mm_stack_c(stack_descr, params,&
       stack_size,&
       a_data, b_data, c_data, used_smm)
    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_4), DIMENSION(*), INTENT(IN)         :: a_data, &
                                                 b_data
    COMPLEX(kind=real_4), DIMENSION(*), INTENT(INOUT)      :: c_data
    LOGICAL, INTENT(OUT)                      :: used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'smm_process_mm_stack_c', &
      routineP = moduleN//':'//routineN

#if defined(__HAS_smm_cnn)

   INTEGER                                   :: sp

   ! TODO we have no way of knowing which calls to libsmm actually resolve to BLAS
   ! Fixing this requires an interface change to libsmm.
   used_smm=.TRUE.

#if defined(__HAS_smm_vec)
    IF(stack_descr%defined_mnk) THEN
       CALL smm_vec_cnn (stack_descr%m, stack_descr%n, stack_descr%k, &
         a_data, b_data, c_data, stack_size, &
         dbcsr_ps_width, params, p_a_first, p_b_first, p_c_first)
       RETURN
    ENDIF
#endif

    DO sp = 1, stack_size
       CALL smm_cnn(&
            params(p_m,sp),&
            params(p_n,sp),&
            params(p_k,sp),&
            a_data(params(p_a_first,sp)),&
            b_data(params(p_b_first,sp)),&
            c_data(params(p_c_first,sp)))
    ENDDO

#else
    ! We do not want to abort here, fall back to BLAS.
    used_smm=.FALSE.
    CALL blas_process_mm_stack_c(params, stack_size,a_data, b_data, c_data)
#endif

    MARK_USED(stack_descr)
  END SUBROUTINE smm_process_mm_stack_c


! **************************************************************************************************
!> \brief Processes MM stack and issues libxsmm calls
!>
!> \param stack_descr ...
!> \param[in] params           Stack of MM parameters
!> \param[in] stack_size       Number of parameters
!> \param[in] a_data           Left-matrix data
!> \param[in] b_data           Right-matrix data
!> \param[in,out] c_data       Product data
!> \param[out] used_smm        Flag to signal if an efficient kernel was used
!> \author Ole Schuett
! **************************************************************************************************
  SUBROUTINE xsmm_process_mm_stack_c(stack_descr, params,&
       stack_size, a_data, b_data, c_data, used_smm)

#if defined(__LIBXSMM) && 0
    ! Caution: This dependency is ignored by makedep.py, because libxsmm.F is kinda empty.
    USE libxsmm,                           ONLY: libxsmm_function  => libxsmm_cmmfunction,&
                                                 libxsmm_dispatch  => libxsmm_cmmdispatch,&
                                                 libxsmm_available => libxsmm_cmmavailable,&
                                                 libxsmm_call      => libxsmm_cmmcall,&
                                                 libxsmm_gemm      => libxsmm_cgemm,&
                                                 LIBXSMM_PREFETCH_NONE,&
                                                 LIBXSMM_PREFETCH,&
                                                 LIBXSMM_ROW_MAJOR,&
                                                 LIBXSMM_COL_MAJOR,&
                                                 LIBXSMM_MAX_MNK,&
                                                 LIBXSMM_FLAGS

    INTEGER, PARAMETER :: LIBXSMM_DEFAULT_PREFETCH = LIBXSMM_PREFETCH
    INTEGER, PARAMETER :: LIBXSMM_DEFAULT_FLAGS = LIBXSMM_FLAGS
#endif

    INTEGER, INTENT(IN)                       :: stack_size
    TYPE(stack_descriptor_type), INTENT(IN)   :: stack_descr
    INTEGER, DIMENSION(dbcsr_ps_width,1:stack_size), &
      INTENT(IN)                              :: params
    COMPLEX(kind=real_4), DIMENSION(*), TARGET, INTENT(IN) :: a_data, b_data
    COMPLEX(kind=real_4), DIMENSION(*), TARGET, &
      INTENT(INOUT)                           :: c_data
    LOGICAL, INTENT(OUT)                      :: used_smm

    CHARACTER(len=*), PARAMETER :: routineN = 'libxsmm_process_mm_stack_c', &
      routineP = moduleN//':'//routineN

#if defined(__LIBXSMM) && 0
    REAL(real_4), PARAMETER                  :: one = 1.0_real_4
    LOGICAL                                   :: processed
    INTEGER(int_8)                            :: threshold
    INTEGER                                   :: fa, fb, fc, m, n, k, pa, pb, pc, sp
    REAL(real_4), DIMENSION(:,:), POINTER    :: a_ptr, b_ptr, c_ptr
    TYPE(libxsmm_function)                    :: func

    processed = .FALSE.
    used_smm = .FALSE.

    CPASSERT(LIBXSMM_COL_MAJOR/=0 .AND. LIBXSMM_ROW_MAJOR==0)

    ! check whether the matrix stack is homogeneous or not
    IF (stack_descr%defined_mnk) THEN
       threshold = INT(stack_descr%m, int_8) * &
                   INT(stack_descr%n, int_8) * &
                   INT(stack_descr%k, int_8)

       ! check if matrices are too large for LIBXSMM (BLAS is likely more efficient)
       IF(threshold > LIBXSMM_MAX_MNK) THEN
          CALL blas_process_mm_stack_c(params, stack_size, a_data, b_data, c_data)
          processed = .TRUE.

       ELSE
          ! try to get a function pointer from libxsmm
          CALL libxsmm_dispatch(func, &
               m=stack_descr%m, n=stack_descr%n, k=stack_descr%k, alpha=one, beta=one, &
               flags=LIBXSMM_DEFAULT_FLAGS, prefetch=LIBXSMM_DEFAULT_PREFETCH)

          IF (libxsmm_available(func)) THEN
             ! load first stack entry
             CPASSERT(stack_size > 0)
             pa = params(p_a_first, 1)
             pb = params(p_b_first, 1)
             pc = params(p_c_first, 1)

             DO sp = 1, stack_size-1
                fa = pa; fb = pb; fc = pc
                ! prefetch next blocks
                pa = params(p_a_first, sp + 1)
                pb = params(p_b_first, sp + 1)
                pc = params(p_c_first, sp + 1)

                ! condition evaluates at compile-time (PARAMETERS)
                IF (LIBXSMM_DEFAULT_PREFETCH /= LIBXSMM_PREFETCH_NONE) THEN
                   CALL libxsmm_call(func, &
                        a=a_data(fa), b=b_data(fb), c=c_data(fc), &
                        ! provide locations of the next operand set
                        pa=a_data(pa), pb=b_data(pb), pc=c_data(pc))
                ELSE
                   CALL libxsmm_call(func, &
                        a=a_data(fa), b=b_data(fb), c=c_data(fc))
                ENDIF
             ENDDO

             ! handle last stack entry without out-of-bounds access
             fa = pa; fb = pb; fc = pc

             ! condition evaluates at compile-time (PARAMETERS)
             IF (LIBXSMM_DEFAULT_PREFETCH /= LIBXSMM_PREFETCH_NONE) THEN
                CALL libxsmm_call(func, &
                     a=a_data(fa), b=b_data(fb), c=c_data(fc), &
                     ! prefetch same blocks
                     pa=a_data(pa), pb=b_data(pb), pc=c_data(pc))
             ELSE
                CALL libxsmm_call(func, &
                     a=a_data(fa), b=b_data(fb), c=c_data(fc))
             ENDIF

             processed = .TRUE.
             used_smm  = .TRUE.
          ENDIF
       ENDIF
    ENDIF


    IF(.NOT.processed) THEN
       ! Dispatch interface was not used, call regular interface.
       ! Should only happen for inhomgenous stacks, then prefetching makes no sense.
       DO sp = 1, stack_size
          m = params(p_m,sp)
          n = params(p_n,sp)
          k = params(p_k,sp)
          fa = params(p_a_first,sp)
          fb = params(p_b_first,sp)
          fc = params(p_c_first,sp)
          ! somewhat expensive pointer remapping required by libxsmm interface
          a_ptr(1:m,1:k) => a_data(fa:fa+(m*k))
          b_ptr(1:k,1:n) => b_data(fb:fb+(k*n))
          c_ptr(1:m,1:n) => c_data(fc:fc+(m*n))
          CALL libxsmm_gemm(m=m, n=n, k=k, a=a_ptr, b=b_ptr, c=c_ptr, &
                            alpha=one, beta=one)
       ENDDO
    ENDIF

#else
    MARK_USED(stack_descr)
    ! We do not want to abort here, fall back to BLAS.
    CALL blas_process_mm_stack_c(params, stack_size,a_data, b_data, c_data)
    used_smm = .FALSE.
#endif

  END SUBROUTINE xsmm_process_mm_stack_c

! **************************************************************************************************
!> \brief ...
!> \param M ...
!> \param N ...
!> \param K ...
!> \param A ...
!> \param B ...
!> \param C ...
! **************************************************************************************************
  PURE SUBROUTINE internal_mm_c_nn(&
       M,N,K,A,B,C)
    INTEGER, INTENT(IN)                      :: M, N, K
    COMPLEX(kind=real_4), INTENT(INOUT)                   :: C(M,N)
    COMPLEX(kind=real_4), INTENT(IN)                      :: B(K,N)
    COMPLEX(kind=real_4), INTENT(IN)                      :: A(M,K)
    C(:,:) = C(:,:) + MATMUL (A, B)
  END SUBROUTINE internal_mm_c_nn
# 316 "/data/isivkov/libdbcsr_svn18247/src/dbcsr/mm/dbcsr_mm_hostdrv.f90"
# 289 "/data/isivkov/libdbcsr_svn18247/src/dbcsr/mm/dbcsr_mm_hostdrv.F" 2

END MODULE dbcsr_mm_hostdrv
