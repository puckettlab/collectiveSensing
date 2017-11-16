/*
GNU Public Licence Copyright (c) Colin Torney
Comments and questions to colin.j.torney@gmail.com

This code is provided freely, however when using this code you are asked to cite our related paper: 
Berdahl, A., Torney, C.J., Ioannou, C.C., Faria, J. & Couzin, I.D. (2013) Emergent sensing of complex environments by mobile animal groups, Science
*/
//#include <cutil_inline.h>
#include <helper_cuda.h>
#include <helper_cuda_gl.h>
#include <helper_functions.h>
#include <cuda.h>
#include <cutil_inline.h>
#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <cufft.h>
#include "noise_kernel.cu"

extern "C"
{


void generateRandomArray(int gridSize, float* d_rands, int* d_seeds, int numBlocks, int numThreads, int numRandBlocks, int numRandThreads, int numLoops)
{
    dim3 dimBlock(numThreads, numThreads);
    dim3 dimGrid(numBlocks, numBlocks);
    int* h_seeds = (int*)malloc(numRandBlocks*sizeof(int));
    float* d_gaussian;

    checkCudaErrors(cudaMalloc((void**)&d_gaussian, sizeof(float)*gridSize*gridSize*2));

    //srand(time(NULL));

    for (int i=0;i<numRandBlocks;i++)
    {
        h_seeds[i]=rand();
  //      printf("%d\n",h_seeds[i]);
    }
    checkCudaErrors(cudaMemcpy(d_seeds, h_seeds, numRandBlocks*sizeof(int), cudaMemcpyHostToDevice));
    gaussianRandsN<<< dimGrid, dimBlock>>>(gridSize,d_gaussian, d_seeds);
//    gaussianRandsD<<< dimGrid, dimBlock>>>(gridSize,d_rands, d_seeds);
	//    uniformRandsD<<< numRandBlocks, numRandThreads>>>(numLoops, d_gaussian, d_seeds);
//    cutilCheckMsg("Kernel execution failed: uniformRandsD");
    
    //convertGaussianD<<<dimGrid, dimBlock>>>(gridSize, d_gaussian);
//    cutilCheckMsg("Kernel execution failed: convertGaussianD");

    createRandArrayD<<<dimGrid, dimBlock>>>(gridSize, d_rands, d_gaussian);
//    cutilCheckMsg("Kernel execution failed: createRandArrayD");

    checkCudaErrors(cudaFree(d_gaussian));
    free (h_seeds);
}

void initializeArray(int gridSize, int numBlocks, int numThreads, float* d_fftArray, float* d_rands, float k0, float scaleFact)
{

    dim3 dimBlock(numThreads, numThreads);
    dim3 dimGrid(numBlocks, numBlocks);

    initializeArrayD<<<dimGrid, dimBlock>>>(gridSize, d_fftArray, d_rands, k0, scaleFact);
//    cutilCheckMsg("Kernel execution failed: initializeArrayD");

}

void advanceArray(int gridSize, int numBlocks, int numThreads, float* d_fftArray, float* d_rands, float k0, float scaleFact, float dT, float revertRate)
{
    dim3 dimBlock(numThreads, numThreads);
    dim3 dimGrid(numBlocks, numBlocks);

    advanceArrayD<<<dimGrid, dimBlock>>>(gridSize, d_fftArray, d_rands, k0, scaleFact, dT, revertRate);
//    cutilCheckMsg("Kernel execution failed: advanceArrayD");
}

void executeFFT(int gridSize, int numBlocks, int numThreads, float* d_fftArray, float* d_tsArray, cufftHandle fftPlan)
{

    dim3 dimBlock(numThreads, numThreads);
    dim3 dimGrid(numBlocks, numBlocks);

    float* d_complexT;

    checkCudaErrors(cudaMalloc((void**)&d_complexT, sizeof(float)*gridSize*gridSize*2)); //??
    cufftExecC2C(fftPlan, (cufftComplex *) d_fftArray, (cufftComplex *) d_complexT, CUFFT_INVERSE) ;  //no cufftsafecall

    convertComplex2RealD<<<dimGrid, dimBlock>>>(gridSize, d_complexT, d_tsArray);
//    cutilCheckMsg("Kernel execution failed: convertComplex2RealD");

    checkCudaErrors(cudaFree(d_complexT)); //???
}

void convertGradients(int gridSize, int numBlocks, int numThreads, float* d_fftArray, float* d_gradArray, cufftHandle fftPlan)
{
    dim3 dimBlock(numThreads, numThreads);
    dim3 dimGrid(numBlocks, numBlocks);

    float* d_xComplexT;
    float* d_yComplexT;
    float* d_xtsArray;
    float* d_ytsArray;

    checkCudaErrors(cudaMalloc((void**)&d_xComplexT, sizeof(float)*gridSize*gridSize*2));
    checkCudaErrors(cudaMalloc((void**)&d_yComplexT, sizeof(float)*gridSize*gridSize*2));

    convertGradD<<<dimGrid, dimBlock>>>(gridSize, d_fftArray, d_xComplexT, d_yComplexT);
//    cutilCheckMsg("Kernel execution failed: convertGradD");
    

    checkCudaErrors(cudaMalloc((void**)&d_xtsArray, sizeof(float)*gridSize*gridSize));
    executeFFT(gridSize, numBlocks, numThreads, d_xComplexT, d_xtsArray, fftPlan);
//    cutilSafeCall(cudaFree(d_xComplexT));

    checkCudaErrors(cudaMalloc((void**)&d_ytsArray, sizeof(float)*gridSize*gridSize));
    executeFFT(gridSize, numBlocks, numThreads, d_yComplexT, d_ytsArray, fftPlan);
    cudaFree(d_xComplexT);
    cudaFree(d_yComplexT);

    mergeVector<<<dimGrid, dimBlock>>>(gridSize, d_xtsArray, d_ytsArray, d_gradArray);
//    cutilCheckMsg("Kernel execution failed: mergeVector");

    checkCudaErrors(cudaFree(d_xtsArray));
    checkCudaErrors(cudaFree(d_ytsArray));
}
}

