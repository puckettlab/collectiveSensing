/*
GNU Public Licence Copyright (c) Colin Torney
Comments and questions to colin.j.torney@gmail.com

This code is provided freely, however when using this code you are asked to cite our related paper: 
Berdahl, A., Torney, C.J., Ioannou, C.C., Faria, J. & Couzin, I.D. (2013) Emergent sensing of complex environments by mobile animal groups, Science
*/
#include <cstdlib>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cstdio>
#include <string.h>

#include "particles_kernel.cu"
#include "cutil_inline.h"

#include <helper_cuda.h>
#include <helper_cuda_gl.h>
#include <helper_functions.h>
#include <helper_timer.h>


extern "C"
{

void copyArrayFromDevice(void* host, const void* device, int size)
{   
    cudaMemcpy(host, device, size, cudaMemcpyDeviceToHost);
}

void copyArrayToDevice(void* device, const void* host, int offset, int size)
{
    cudaMemcpy((char *) device + offset, host, size, cudaMemcpyHostToDevice);
}

void setParameters(SimParams *hostParams)
{
    // copy parameters to constant memory
	cudaMemcpyToSymbol(params, hostParams, sizeof(SimParams)) ;
     cudaMemcpyToSymbol(params, hostParams, sizeof(SimParams)) ;
}

void convert4to3DArray(float* input, float* output, int numParticles,  int numThreads, int numBlocks)
{
    // execute the kernel
    convert4to3DArrayD<<< numBlocks, numThreads >>>( (float4*)input, (float3*)output, numParticles);
    
    return;
}

void processInteractions(float* curPos, float* curVel, float* darkArray, float* dControl, float* dRands, float* gRands,uint numParticles, float xpatch, float ypatch, int numThreads, int numBlocks)
{
    // calls the interaction function and updates positions and velocities

    // process interactions
    interactD<<< numBlocks, numThreads >>>((float4*)curVel, (float4*)curPos, numParticles, darkArray, dControl, xpatch, ypatch, gRands); 
//    cudaMemcpy(dRands, m_dRands, numBlocks*sizeof(float), cudaMemcpyHostToDevice);// broke
    addRotationalNoiseD<<< numBlocks, numThreads >>>((float4*)curVel, numParticles, dRands);
//    cudaBindTexture(0, dRandTex, dRands, numParticles*sizeof(float));	

}


void generateRandomNumbers(float* d_rands, int* d_seeds, int numBlocks, int numThreads)
{
    int* h_seeds = (int*)malloc(numBlocks*sizeof(int));

    for (int i=0;i<numBlocks;i++)
        h_seeds[i]=rand();

    copyArrayToDevice(d_seeds, h_seeds, 0, numBlocks*sizeof(int));
    cudaMemcpy(d_seeds, h_seeds, numBlocks*sizeof(int), cudaMemcpyHostToDevice);

    gaussianRandsD<<< numBlocks, numThreads>>>(d_rands, d_seeds);

    free (h_seeds);
}
}   // extern "C"

