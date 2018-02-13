# 1 "/data/isivkov/libdbcsr_svn18247/src/acc/acc_hostmem.F"
!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Accelerator support
!> \author  Ole Schuett
!> \date    2013-04
! **************************************************************************************************
MODULE acc_hostmem
#if defined (__ACC)
    USE ISO_C_BINDING, ONLY: C_INT, C_SIZE_T, C_PTR, C_LOC, C_F_POINTER
#endif
  USE kinds,                           ONLY: int_4,&
                                             int_4_size,&
                                             int_8,&
                                             int_8_size,&
                                             real_4,&
                                             real_4_size,&
                                             real_8,&
                                             real_8_size
  USE acc_stream,                      ONLY: acc_stream_associated,&
                                             acc_stream_cptr,&
                                             acc_stream_type
#include "../base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'acc_hostmem'

   LOGICAL, PARAMETER :: careful_mod = .TRUE.

   PUBLIC :: acc_hostmem_allocate, acc_hostmem_deallocate

   INTERFACE acc_hostmem_allocate
      MODULE PROCEDURE acc_hostmem_alloc_i4, acc_hostmem_alloc_i8
      MODULE PROCEDURE acc_hostmem_alloc_r4, acc_hostmem_alloc_r8
      MODULE PROCEDURE acc_hostmem_alloc_c4, acc_hostmem_alloc_c8
      MODULE PROCEDURE acc_hostmem_alloc_i4_2D, acc_hostmem_alloc_i8_2D
      MODULE PROCEDURE acc_hostmem_alloc_r4_2D, acc_hostmem_alloc_r8_2D
      MODULE PROCEDURE acc_hostmem_alloc_c4_2D, acc_hostmem_alloc_c8_2D
   END INTERFACE

   INTERFACE acc_hostmem_deallocate
      MODULE PROCEDURE acc_hostmem_dealloc_i4, acc_hostmem_dealloc_i8
      MODULE PROCEDURE acc_hostmem_dealloc_r4, acc_hostmem_dealloc_r8
      MODULE PROCEDURE acc_hostmem_dealloc_c4, acc_hostmem_dealloc_c8
      MODULE PROCEDURE acc_hostmem_dealloc_i4_2D, acc_hostmem_dealloc_i8_2D
      MODULE PROCEDURE acc_hostmem_dealloc_r4_2D, acc_hostmem_dealloc_r8_2D
      MODULE PROCEDURE acc_hostmem_dealloc_c4_2D, acc_hostmem_dealloc_c8_2D
   END INTERFACE

#if defined (__ACC)

   INTERFACE
      FUNCTION cuda_host_mem_alloc_cu(mem, n, stream_ptr) RESULT(istat) BIND(C, name="acc_host_mem_allocate")
         IMPORT
         TYPE(C_PTR)                              :: mem
         INTEGER(KIND=C_SIZE_T), INTENT(IN), &
            VALUE                                  :: n
         TYPE(C_PTR), VALUE                       :: stream_ptr
         INTEGER(KIND=C_INT)                      :: istat

      END FUNCTION cuda_host_mem_alloc_cu
   END INTERFACE

   INTERFACE
      FUNCTION cuda_host_mem_dealloc_cu(mem, stream_ptr) RESULT(istat) bind(C, name="acc_host_mem_deallocate")
         IMPORT
         TYPE(C_PTR), VALUE                       :: mem, stream_ptr
         INTEGER(KIND=C_INT)                      :: istat

      END FUNCTION cuda_host_mem_dealloc_cu
   END INTERFACE

#endif

CONTAINS


! **************************************************************************************************
!> \brief Helper-routine performing allocation of host-pinned cuda memory.
!> \param host_mem_c_ptr pointer to allocated memory
!> \param n_bytes number of bytes to allocate
!> \param stream ...
! **************************************************************************************************
#if defined (__ACC)
   SUBROUTINE acc_hostmem_alloc_raw(host_mem_c_ptr, n_bytes, stream)
      TYPE(C_PTR), INTENT(OUT)                           :: host_mem_c_ptr
      INTEGER, INTENT(IN)                                :: n_bytes
      TYPE(acc_stream_type), INTENT(IN)                  :: stream

      CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_alloc_raw', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: istat
      TYPE(C_PTR)                                        :: stream_cptr

      IF (.NOT. acc_stream_associated(stream)) &
         CPABORT("acc_hostmem_alloc_raw: stream not associated")

      stream_cptr = acc_stream_cptr(stream)

      istat = cuda_host_mem_alloc_cu(host_mem_c_ptr, INT(n_bytes, KIND=C_SIZE_T), stream_cptr)
      IF (istat /= 0) &
         CPABORT("acc_hostmem_alloc_raw: Could not allocate host pinned memory")
   END SUBROUTINE acc_hostmem_alloc_raw
#endif

#if defined (__ACC)
! **************************************************************************************************
!> \brief ...
!> \param host_mem_c_ptr ...
!> \param stream ...
! **************************************************************************************************
   SUBROUTINE acc_hostmem_dealloc_raw(host_mem_c_ptr, stream)
      TYPE(C_PTR), INTENT(IN)                            :: host_mem_c_ptr
      TYPE(acc_stream_type), INTENT(IN)                  :: stream

      CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_raw', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: istat
      TYPE(C_PTR)                                        :: stream_cptr

      IF (.NOT. acc_stream_associated(stream)) &
         CPABORT("acc_hostmem_dealloc_raw: stream not associated")

      stream_cptr = acc_stream_cptr(stream)

      istat = cuda_host_mem_dealloc_cu(host_mem_c_ptr, stream_cptr)
      IF (istat /= 0) &
         CPABORT("acc_hostmem_dealloc_raw: Could not deallocate host pinned memory")
   END SUBROUTINE acc_hostmem_dealloc_raw
#endif


# 147 "/data/isivkov/libdbcsr_svn18247/src/acc/acc_hostmem.F"

# 149 "/data/isivkov/libdbcsr_svn18247/src/acc/acc_hostmem.F"

! **************************************************************************************************
!> \brief Allocates 1D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n size given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_alloc_i4 (host_mem, n, stream)
    INTEGER(kind=int_4), DIMENSION(:), POINTER          :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr

    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, MAX(1,n)*int_4_size, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n) /))
#else
    MARK_USED(host_mem)
    MARK_USED(n)
    MARK_USED(stream)
    CPABORT("acc_hostmem_alloc_i4: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_alloc_i4



! **************************************************************************************************
!> \brief Allocates 2D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n1 sizes given in terms of item-count (not bytes!)
!> \param n2 sizes given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_alloc_i4_2D (host_mem, n1, n2, stream)
    INTEGER(kind=int_4), DIMENSION(:,:), POINTER        :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes

    n_bytes = MAX(1,n1)*MAX(1,n2)*int_4_size
    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, n_bytes, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n1),MAX(1,n2) /))
#else
    MARK_USED(host_mem)
    MARK_USED(n1)
    MARK_USED(n2)
    MARK_USED(stream)
    CPABORT("acc_hostmem_alloc_i4_2D: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_alloc_i4_2D


! **************************************************************************************************
!> \brief Deallocates a 1D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_dealloc_i4 (host_mem, stream)
    INTEGER(kind=int_4), DIMENSION(:), POINTER          :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_i4', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1)), stream)
#else
    MARK_USED(host_mem)
    MARK_USED(stream)
    CPABORT("acc_hostmem_dealloc_i4: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_dealloc_i4


! **************************************************************************************************
!> \brief Deallocates a 2D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_dealloc_i4_2D (host_mem, stream)
    INTEGER(kind=int_4), DIMENSION(:,:), POINTER        :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_i4_2D', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1,1)), stream)
#else
    MARK_USED(host_mem)
    MARK_USED(stream)
    CPABORT("acc_hostmem_dealloc_i4: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_dealloc_i4_2D

# 149 "/data/isivkov/libdbcsr_svn18247/src/acc/acc_hostmem.F"

! **************************************************************************************************
!> \brief Allocates 1D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n size given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_alloc_i8 (host_mem, n, stream)
    INTEGER(kind=int_8), DIMENSION(:), POINTER          :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr

    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, MAX(1,n)*int_8_size, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n) /))
#else
    MARK_USED(host_mem)
    MARK_USED(n)
    MARK_USED(stream)
    CPABORT("acc_hostmem_alloc_i8: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_alloc_i8



! **************************************************************************************************
!> \brief Allocates 2D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n1 sizes given in terms of item-count (not bytes!)
!> \param n2 sizes given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_alloc_i8_2D (host_mem, n1, n2, stream)
    INTEGER(kind=int_8), DIMENSION(:,:), POINTER        :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes

    n_bytes = MAX(1,n1)*MAX(1,n2)*int_8_size
    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, n_bytes, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n1),MAX(1,n2) /))
#else
    MARK_USED(host_mem)
    MARK_USED(n1)
    MARK_USED(n2)
    MARK_USED(stream)
    CPABORT("acc_hostmem_alloc_i8_2D: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_alloc_i8_2D


! **************************************************************************************************
!> \brief Deallocates a 1D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_dealloc_i8 (host_mem, stream)
    INTEGER(kind=int_8), DIMENSION(:), POINTER          :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_i8', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1)), stream)
#else
    MARK_USED(host_mem)
    MARK_USED(stream)
    CPABORT("acc_hostmem_dealloc_i8: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_dealloc_i8


! **************************************************************************************************
!> \brief Deallocates a 2D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_dealloc_i8_2D (host_mem, stream)
    INTEGER(kind=int_8), DIMENSION(:,:), POINTER        :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_i8_2D', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1,1)), stream)
#else
    MARK_USED(host_mem)
    MARK_USED(stream)
    CPABORT("acc_hostmem_dealloc_i8: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_dealloc_i8_2D

# 149 "/data/isivkov/libdbcsr_svn18247/src/acc/acc_hostmem.F"

! **************************************************************************************************
!> \brief Allocates 1D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n size given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_alloc_r4 (host_mem, n, stream)
    REAL(kind=real_4), DIMENSION(:), POINTER          :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr

    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, MAX(1,n)*real_4_size, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n) /))
#else
    MARK_USED(host_mem)
    MARK_USED(n)
    MARK_USED(stream)
    CPABORT("acc_hostmem_alloc_r4: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_alloc_r4



! **************************************************************************************************
!> \brief Allocates 2D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n1 sizes given in terms of item-count (not bytes!)
!> \param n2 sizes given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_alloc_r4_2D (host_mem, n1, n2, stream)
    REAL(kind=real_4), DIMENSION(:,:), POINTER        :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes

    n_bytes = MAX(1,n1)*MAX(1,n2)*real_4_size
    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, n_bytes, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n1),MAX(1,n2) /))
#else
    MARK_USED(host_mem)
    MARK_USED(n1)
    MARK_USED(n2)
    MARK_USED(stream)
    CPABORT("acc_hostmem_alloc_r4_2D: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_alloc_r4_2D


! **************************************************************************************************
!> \brief Deallocates a 1D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_dealloc_r4 (host_mem, stream)
    REAL(kind=real_4), DIMENSION(:), POINTER          :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_r4', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1)), stream)
#else
    MARK_USED(host_mem)
    MARK_USED(stream)
    CPABORT("acc_hostmem_dealloc_r4: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_dealloc_r4


! **************************************************************************************************
!> \brief Deallocates a 2D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_dealloc_r4_2D (host_mem, stream)
    REAL(kind=real_4), DIMENSION(:,:), POINTER        :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_r4_2D', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1,1)), stream)
#else
    MARK_USED(host_mem)
    MARK_USED(stream)
    CPABORT("acc_hostmem_dealloc_r4: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_dealloc_r4_2D

# 149 "/data/isivkov/libdbcsr_svn18247/src/acc/acc_hostmem.F"

! **************************************************************************************************
!> \brief Allocates 1D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n size given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_alloc_r8 (host_mem, n, stream)
    REAL(kind=real_8), DIMENSION(:), POINTER          :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr

    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, MAX(1,n)*real_8_size, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n) /))
#else
    MARK_USED(host_mem)
    MARK_USED(n)
    MARK_USED(stream)
    CPABORT("acc_hostmem_alloc_r8: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_alloc_r8



! **************************************************************************************************
!> \brief Allocates 2D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n1 sizes given in terms of item-count (not bytes!)
!> \param n2 sizes given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_alloc_r8_2D (host_mem, n1, n2, stream)
    REAL(kind=real_8), DIMENSION(:,:), POINTER        :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes

    n_bytes = MAX(1,n1)*MAX(1,n2)*real_8_size
    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, n_bytes, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n1),MAX(1,n2) /))
#else
    MARK_USED(host_mem)
    MARK_USED(n1)
    MARK_USED(n2)
    MARK_USED(stream)
    CPABORT("acc_hostmem_alloc_r8_2D: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_alloc_r8_2D


! **************************************************************************************************
!> \brief Deallocates a 1D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_dealloc_r8 (host_mem, stream)
    REAL(kind=real_8), DIMENSION(:), POINTER          :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_r8', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1)), stream)
#else
    MARK_USED(host_mem)
    MARK_USED(stream)
    CPABORT("acc_hostmem_dealloc_r8: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_dealloc_r8


! **************************************************************************************************
!> \brief Deallocates a 2D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_dealloc_r8_2D (host_mem, stream)
    REAL(kind=real_8), DIMENSION(:,:), POINTER        :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_r8_2D', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1,1)), stream)
#else
    MARK_USED(host_mem)
    MARK_USED(stream)
    CPABORT("acc_hostmem_dealloc_r8: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_dealloc_r8_2D

# 149 "/data/isivkov/libdbcsr_svn18247/src/acc/acc_hostmem.F"

! **************************************************************************************************
!> \brief Allocates 1D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n size given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_alloc_c4 (host_mem, n, stream)
    COMPLEX(kind=real_4), DIMENSION(:), POINTER          :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr

    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, MAX(1,n)*2*real_4_size, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n) /))
#else
    MARK_USED(host_mem)
    MARK_USED(n)
    MARK_USED(stream)
    CPABORT("acc_hostmem_alloc_c4: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_alloc_c4



! **************************************************************************************************
!> \brief Allocates 2D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n1 sizes given in terms of item-count (not bytes!)
!> \param n2 sizes given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_alloc_c4_2D (host_mem, n1, n2, stream)
    COMPLEX(kind=real_4), DIMENSION(:,:), POINTER        :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes

    n_bytes = MAX(1,n1)*MAX(1,n2)*2*real_4_size
    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, n_bytes, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n1),MAX(1,n2) /))
#else
    MARK_USED(host_mem)
    MARK_USED(n1)
    MARK_USED(n2)
    MARK_USED(stream)
    CPABORT("acc_hostmem_alloc_c4_2D: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_alloc_c4_2D


! **************************************************************************************************
!> \brief Deallocates a 1D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_dealloc_c4 (host_mem, stream)
    COMPLEX(kind=real_4), DIMENSION(:), POINTER          :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_c4', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1)), stream)
#else
    MARK_USED(host_mem)
    MARK_USED(stream)
    CPABORT("acc_hostmem_dealloc_c4: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_dealloc_c4


! **************************************************************************************************
!> \brief Deallocates a 2D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_dealloc_c4_2D (host_mem, stream)
    COMPLEX(kind=real_4), DIMENSION(:,:), POINTER        :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_c4_2D', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1,1)), stream)
#else
    MARK_USED(host_mem)
    MARK_USED(stream)
    CPABORT("acc_hostmem_dealloc_c4: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_dealloc_c4_2D

# 149 "/data/isivkov/libdbcsr_svn18247/src/acc/acc_hostmem.F"

! **************************************************************************************************
!> \brief Allocates 1D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n size given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_alloc_c8 (host_mem, n, stream)
    COMPLEX(kind=real_8), DIMENSION(:), POINTER          :: host_mem
    INTEGER, INTENT(IN)                      :: n
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr

    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, MAX(1,n)*2*real_8_size, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n) /))
#else
    MARK_USED(host_mem)
    MARK_USED(n)
    MARK_USED(stream)
    CPABORT("acc_hostmem_alloc_c8: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_alloc_c8



! **************************************************************************************************
!> \brief Allocates 2D fortan-array as cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param n1 sizes given in terms of item-count (not bytes!)
!> \param n2 sizes given in terms of item-count (not bytes!)
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_alloc_c8_2D (host_mem, n1, n2, stream)
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER        :: host_mem
    INTEGER, INTENT(IN)                      :: n1, n2
    TYPE(acc_stream_type), INTENT(IN)        :: stream
#if defined (__ACC)
    TYPE(C_PTR)                              :: host_mem_c_ptr
    INTEGER                                  :: n_bytes

    n_bytes = MAX(1,n1)*MAX(1,n2)*2*real_8_size
    CALL acc_hostmem_alloc_raw(host_mem_c_ptr, n_bytes, stream)
    CALL C_F_POINTER (host_mem_c_ptr, host_mem, (/ MAX(1,n1),MAX(1,n2) /))
#else
    MARK_USED(host_mem)
    MARK_USED(n1)
    MARK_USED(n2)
    MARK_USED(stream)
    CPABORT("acc_hostmem_alloc_c8_2D: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_alloc_c8_2D


! **************************************************************************************************
!> \brief Deallocates a 1D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_dealloc_c8 (host_mem, stream)
    COMPLEX(kind=real_8), DIMENSION(:), POINTER          :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_c8', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1)), stream)
#else
    MARK_USED(host_mem)
    MARK_USED(stream)
    CPABORT("acc_hostmem_dealloc_c8: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_dealloc_c8


! **************************************************************************************************
!> \brief Deallocates a 2D fortan-array, which is cuda host-pinned memory.
!> \param host_mem pointer to array
!> \param stream ...
!> \author  Ole Schuett
! **************************************************************************************************
  SUBROUTINE acc_hostmem_dealloc_c8_2D (host_mem, stream)
    COMPLEX(kind=real_8), DIMENSION(:,:), POINTER        :: host_mem
    TYPE(acc_stream_type), INTENT(IN)        :: stream
    CHARACTER(len=*), PARAMETER :: routineN = 'acc_hostmem_dealloc_c8_2D', &
      routineP = moduleN//':'//routineN

    IF (SIZE (host_mem) == 0) RETURN
#if defined (__ACC)
    CALL acc_hostmem_dealloc_raw(C_LOC(host_mem(1,1)), stream)
#else
    MARK_USED(host_mem)
    MARK_USED(stream)
    CPABORT("acc_hostmem_dealloc_c8: ACC not compiled in.")
#endif
  END SUBROUTINE acc_hostmem_dealloc_c8_2D

# 251 "/data/isivkov/libdbcsr_svn18247/src/acc/acc_hostmem.F"

END MODULE acc_hostmem
