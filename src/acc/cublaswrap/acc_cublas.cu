#if defined(__ACC_CUBLAS)
#include <cublas.h>
#endif

extern "C" int acc_dgemm (int m, int n, int k, int transa, int transb, 
                          void *a_data, int a_lb, void *b_data, int b_lb, void *c_data, int c_lb, 
			  double alpha, double beta, void *stream)
{
#if defined(__ACC_CUBLAS)
  cudaStream_t* custream = (cudaStream_t*) stream;
  double *a = (double*)a_data;
  double *b = (double*)b_data;
  double *c = (double*)c_data;

  // Check transposity
  char transa_char = 'N'; 
  int lda = m;
  if (transa!=0) {
    transa_char = 'T';
    lda = k;
  }
  char transb_char = 'N'; 
  int ldb = k;
  if (transb!=0) {
    transb_char = 'T';
    ldb = n;
  }

  cublasGetError(); // Clear error
  cublasSetKernelStream(*custream);  
  cublasDgemm(transa_char, transb_char, m, n, k, alpha, &a[a_lb], lda, &b[b_lb], ldb, beta, &c[c_lb], m);
  int status = cublasGetError();
  if (status != CUBLAS_STATUS_SUCCESS) return(-1);

  return(0); 
#else
  return(-1);
#endif
}

