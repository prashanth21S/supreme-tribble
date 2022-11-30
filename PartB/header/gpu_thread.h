#include <stdio.h>
#include <cuda.h>
#include <stdlib.h>

//Kernel methods
__global__ void matrixMul(const int *a, const int *b, int *c, int N){
  int row = blockIdx.y * blockDim.y + threadIdx.y;
  int col = blockIdx.x * blockDim.x + threadIdx.x;
  int rowA=row;
  int colB=col;

// Do these procedure for only even rows and columns respectly
if (((row %2)==0) && (col%2==0) && (row+1<N) && (col+1<N))
{
 int sum=0;
    for(int iter = 0; iter < N; iter++) 
     {
       sum += a[rowA * N + iter] * b[iter * N + colB];
       sum += a[(rowA+1) * N + iter] * b[iter * N + colB];
       sum += a[rowA * N + iter] * b[iter * N + (colB+1)];
       sum += a[(rowA+1) * N + iter] * b[iter * N + (colB+1)];
     }
     int rowC = rowA>>1;
     int colC = colB>>1;
     int indexC = rowC * (N>>1) + colC;
     c[indexC] = sum;
}
}


void gpuThread(int N, int *matA, int *matB, int *output)
{
  // allocate device variables
  int *d_a, *d_b, *d_c;
  cudaMalloc((void**)&d_a, N*N*sizeof(int));
  cudaMalloc((void**)&d_b, N*N*sizeof(int));
  cudaMalloc((void**)&d_c, N/2*N/2*sizeof(int));

  // Copy data to the device
  cudaMemcpy(d_a, matA, N*N*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_b, matB, N*N*sizeof(int), cudaMemcpyHostToDevice);

  // Number of threads per blocks
  int threads_per_block = 32;

  // Number of blocks
  int blocks_total = N/threads_per_block ;

  // Set dim3 tructure to use thr
   dim3 threads(threads_per_block, threads_per_block);
   dim3 blocks(blocks_total, blocks_total);

  // Call the Kernel
  matrixMul<<<blocks, threads>>>(d_a,d_b,d_c,N);

  // Copy back to the host
  cudaMemcpy(output, d_c, N/2*N/2*sizeof(int),cudaMemcpyDeviceToHost);

  // Free memory on device
  cudaFree(d_a);
  cudaFree(d_b);
  cudaFree(d_c);

}
