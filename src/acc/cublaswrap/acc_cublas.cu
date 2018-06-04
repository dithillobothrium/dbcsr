#include <cuda_runtime.h>
#include <cuda.h>
#include <cuda_device_runtime_api.h>
#include <cublas_v2.h>
#include <stdio.h>

#ifdef __ACC_CUBLAS

extern "C" int acc_dgemm (cublasHandle_t *handle, int m, int n, int k, int transa, int transb, 
                          void *a_data, int a_lb, void *b_data, int b_lb, void *c_data, int c_lb, 
			  double alpha, double beta, void *stream)
{
  cudaStream_t* custream = (cudaStream_t*) stream;
  double *a = (double*)a_data;
  double *b = (double*)b_data;
  double *c = (double*)c_data;
/*  
  double *host;
  int sz = m*k > n*k ? n*k : m*k;
  int print_sz = 50;

  host = (double*)malloc(sizeof(double)*sz);
  
  cudaMemcpy(host,a,sizeof(double)*sz,cudaMemcpyDeviceToHost);
  for (int i=0; i < (print_sz < sz ? print_sz: sz) ; i++) printf("%f ",host[i]);
  printf("\n");
   
  cudaMemcpy(host,b,sizeof(double)*sz,cudaMemcpyDeviceToHost);
  for (int i=0; i < (print_sz < sz ? print_sz: sz) ; i++) printf("%f ",host[i]);
  printf("\n");
  
  free(host);

*/  
   // Check transposity
  char transa_char = 'N';
  cublasOperation_t trana = CUBLAS_OP_N;
  int lda = m;

  if (transa!=0) {
    transa_char = 'T';
    trana = CUBLAS_OP_T;
    lda = k;
  }

  char transb_char = 'N';
  cublasOperation_t tranb = CUBLAS_OP_N;
  int ldb = k;

  if (transb!=0) {
    transb_char = 'T';
    tranb = CUBLAS_OP_T;
    ldb = n;
  }
 
//  printf("in C-cublas, transa = %c transb = %c \n",transa_char,transb_char);

//  printf("m n k = %d %d %d, alpha beta = %f %f, a_lb b_lb c_lb = %d %d %d\n",m,n,k,alpha,beta,a_lb, b_lb, c_lb);

  cublasSetStream(*handle, *custream);  
   
  cublasStatus_t stat = cublasDgemm(*handle, trana, tranb, m, n, k, &alpha, &a[a_lb], lda, &b[b_lb], ldb, &beta, &c[c_lb], m);
    
  if (stat != CUBLAS_STATUS_SUCCESS) return(-1);

//  cudaDeviceSynchronize();

//  printf("cublas OK");
  return(0); 
}


// cublas interface ------------------
 
extern "C" int f_cublasCreate(cublasHandle_t **handle)
{
    *handle = (cublasHandle_t*)malloc(sizeof(cublasHandle_t));
    return cublasCreate(*handle);
}
 
extern "C" int f_cublasDgemm(cublasHandle_t *handle,
               cublasOperation_t transa, cublasOperation_t transb,
              int m, int n, int k,
              const double *alpha,
              const double *A, int lda,
              const double *B, int ldb,
              const double *beta,
              double *C, int ldc)
{
    return cublasDgemm(*handle,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc);
}
 
extern "C" int f_cublasDgemmBatched(cublasHandle_t *handle,
               cublasOperation_t transa, cublasOperation_t transb,
              int m, int n, int k,
              const double *alpha,
              const double **A, int lda,
              const double **B, int ldb,
              const double *beta,
              double **C, int ldc,
              int batch_count)
{
    return cublasDgemmBatched(*handle,transa,transb,m,n,k,alpha,A,lda,B,ldb,beta,C,ldc,batch_count);
}
 
extern "C" void f_cublasDestroy(cublasHandle_t *handle)
{
    cublasDestroy(*handle);
    free(handle);
}
 
extern "C" int f_cudaStreamCreate(cudaStream_t **stream)
{
    *stream = (cudaStream_t *) malloc(sizeof(cudaStream_t));
    return cudaStreamCreate(*stream);
}
 
extern "C" int f_cublasSetStream(cublasHandle_t *handle, cudaStream_t *streamid)
{
    return cublasSetStream(*handle, *streamid);
}
 
extern "C" void f_cudaStreamDestroy(cudaStream_t *stream)
{
    cudaStreamDestroy(*stream);
}

#endif

