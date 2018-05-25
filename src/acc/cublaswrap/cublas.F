MODULE cublas
#if defined(__ACC_CUBLAS)
   USE ISO_C_BINDING, ONLY: C_INT, &
      C_PTR, &
      C_DOUBLE
#endif

   IMPLICIT NONE

#if defined(__ACC_CUBLAS)
   PUBLIC :: acc_dgemm_cu

   INTERFACE
      FUNCTION acc_dgemm_cu(m, n, k, transa, transb, a_data, a_lb, b_data, b_lb, c_data, c_lb, &
                            alpha, beta, stream_ptr) &
         RESULT(istat) &
         BIND(C, name="acc_dgemm")
         IMPORT
         INTEGER(KIND=C_INT), INTENT(IN), VALUE   :: m, n, k, a_lb, b_lb, c_lb, transa, transb
         TYPE(C_PTR), INTENT(IN), VALUE           :: a_data, b_data, c_data
         TYPE(C_PTR), VALUE                       :: stream_ptr
         REAL(KIND=C_DOUBLE), VALUE               :: alpha, beta
         INTEGER(KIND=C_INT)                      :: istat
      END FUNCTION acc_dgemm_cu
   END INTERFACE
#endif

END MODULE cublas
