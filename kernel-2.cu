
#include "common.h"
#include "timer.h"

#define IN_TILE_DIM 32
#define OUT_TILE_DIM ((IN_TILE_DIM) - 2*(FILTER_RADIUS))

__constant__ float filter_c[FILTER_DIM][FILTER_DIM];

__global__ void convolution_tiled_kernel(float* input, float* output, unsigned int width, unsigned int height) {

     __shared__ float cov[IN_TILE_DIM][IN_TILE_DIM];
     unsigned int row = blockIdx.y*blockDim.y + threadIdx.y;
     unsigned int col = blockIdx.x*blockDim.x + threadIdx.x;
     float sum = 0.0f;
     for(unsigned int tile = 0; tile < (width + TILE_DIM -1 ) / TILE_DIM; ++tile) {
        if((row >= 0) && (row< height) && (col >= 0) && (col < width) ) {
        cov[threadIdx.y][threadIdx.x]= cov[row*width + tile*TILE_DIM + threadIdx.x];
    }
    else{
     cov[threadIdx.y][threadIdx.x]=0;
     }
     if(threadIdx.y < TILE_DIM && threadIdx.x < TILE_DIM){
        for(i = 0; i < FILTER_DIM; i++) {
            for(j = 0; j < FILTER_DIM; j++) { 
                output += filter_c_[i][j] * cov[i+threadIdx.y][j+threadIdx.x];
} }
     
    output[outRow*width + outCol] = sum;
     }
     
     









}

void copyFilterToGPU(float filter[][FILTER_DIM]) {

    // Copy filter to constant memory

    cudaMemcpyToSymbol(filter_c, filter, FILTER_DIM*FILTER_DIM*sizeof(float));

}

void convolution_tiled_gpu(float* input_d, float* output_d, unsigned int width, unsigned int height) {

    // Call kernel

    dim3 numThreadsPerBlock(OUT_TILE_DIM, OUT_TILE_DIM);
    dim3 numBlocks((width + OUT_TILE_DIM - 1)/OUT_TILE_DIM, (height + OUT_TILE_DIM - 1)/OUT_TILE_DIM);
    convolution_tiled_kernel <<< numBlocks, numThreadsPerBlock >>> (input_d, output_d, width, height);



}

