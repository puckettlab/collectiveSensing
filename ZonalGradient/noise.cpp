/*

This code is provided freely, however when using this code you are asked to cite these related paper: 
Puckett, J.G., Pokhrel, A.R. and Giannini, J.A. (2017) Collective gradient sensing in fish schools

This code is modified from original work by 
Berdahl, A., Torney, C.J., Ioannou, C.C., Faria, J. & Couzin, I.D. (2013) Emergent sensing of complex environments by mobile animal groups, Science

GNU Public Licence Copyright (c) Colin Torney
Comments and questions to colin.j.torney@gmail.com

*/
#include <iostream>
#include <string>
#include <noise.h>

/*
--------------------------------------------------------------------------------
	Code for the stochastic gradient field
--------------------------------------------------------------------------------
*/

using namespace std;
clsField::clsField(int gridSize_in, double dt_in, int seed_in, double visc_in)
{
    gridSize=gridSize_in;
    dt=dt_in;  
    seed = seed_in;//revertRate=visc

    numThreads=min(16,gridSize);
    numBlocks=gridSize/numThreads;

    numRandBlocks=32;
    numRandThreads=min(256,gridSize);
    numRandLoops=(2*gridSize*gridSize)/(numRandBlocks*numRandThreads);

    while (numRandLoops<1)
    {
        numRandBlocks/=2;
        numRandLoops=(2*gridSize*gridSize)/(numRandBlocks*numRandThreads);
    }

    checkCudaErrors(cudaMalloc((void**)&d_seeds, numRandBlocks*sizeof(int)));
    checkCudaErrors(cudaMalloc((void**)&d_rands, sizeof(float)*gridSize*gridSize*2));
    checkCudaErrors(cudaMalloc((void**)&d_fftArray, sizeof(float)*gridSize*gridSize*2));
    checkCudaErrors(cudaMalloc((void**)&d_tsArray, sizeof(float)*gridSize*gridSize));
    checkCudaErrors(cudaMalloc((void**)&d_gradArray, sizeof(float)*gridSize*gridSize*2));

    checkCudaErrors(cudaMemset((void*)d_seeds, 0, sizeof(int)*numRandBlocks));
    checkCudaErrors(cudaMemset((void*)d_rands, 0, sizeof(float)*gridSize*gridSize*2));
    checkCudaErrors(cudaMemset((void*)d_fftArray, 0, sizeof(float)*gridSize*gridSize*2));
    checkCudaErrors(cudaMemset((void*)d_tsArray, 0, sizeof(float)*gridSize*gridSize));
    checkCudaErrors(cudaMemset((void*)d_gradArray, 0, sizeof(float)*gridSize*gridSize*2));

    checkCudaErrors(cufftPlan2d(&fftPlan, gridSize, gridSize, CUFFT_C2C) );

    h_tsArray = new float[gridSize*gridSize];
    h_gradArray = new float[2*gridSize*gridSize];

    revertRate=visc_in;//visc
    k0=16.84;
    //k0=24;

    double ki,kj,ksq;
    double lambda_k=0.0;
    for (int i=1;i<gridSize;i++)
    {
        ki=2*M_PI*(i-(gridSize*0.5));
        for (int j=1;j<gridSize;j++)
        {
            kj=2*M_PI*(j-(gridSize*0.5));
            ksq=((ki*ki)+(kj*kj));
            lambda_k+=getSpectrum(k0,ksq);
        }
    }
    scaleFact=1.0/(lambda_k);
    //printf("lambda  %4.4f\n",lambda_k);
    //printf("scaleF  %4.4f\n",scaleFact);
    generateRandomArray(gridSize, d_rands, d_seeds, numBlocks, numThreads, numRandBlocks, numRandThreads, numRandLoops);
    initializeArray(gridSize, numBlocks, numThreads, d_fftArray, d_rands, k0, scaleFact);
    convertGradients(gridSize, numBlocks, numThreads, d_fftArray, d_gradArray, fftPlan);
}

void clsField::advanceTimestep()
{
    generateRandomArray(gridSize, d_rands, d_seeds, numBlocks, numThreads, numRandBlocks, numRandThreads, numRandLoops);
    advanceArray(gridSize, numBlocks, numThreads, d_fftArray, d_rands, k0, scaleFact, dt, revertRate);
    convertGradients(gridSize, numBlocks, numThreads, d_fftArray, d_gradArray, fftPlan);
    executeFFT(gridSize, numBlocks, numThreads, d_fftArray, d_tsArray, fftPlan);
	//printf("%4.4f\t%4.4f\n",&d_rands[0],&d_tsArray[0]);
}

void clsField::prepareSave()
{
    executeFFT(gridSize, numBlocks, numThreads, d_fftArray, d_tsArray, fftPlan);
    checkCudaErrors(cudaMemcpy(h_tsArray, d_tsArray, gridSize*gridSize*sizeof(float), cudaMemcpyDeviceToHost));
/* 
   float maxNoise=0;
    for (int i=0;i<gridSize;i++)
    {
        for (int j=0;j<gridSize;j++)
            if (maxNoise<h_tsArray[j*gridSize+i])
                maxNoise=h_tsArray[j*gridSize+i];

    }
*/
    //    h_tsArray[j*gridSize+i]=((float)i/(float)gridSize);
}

clsField::~clsField()
{
    delete [] h_tsArray;
    delete [] h_gradArray;

    checkCudaErrors(cudaFree(d_rands));
    checkCudaErrors(cudaFree(d_seeds));
    checkCudaErrors(cudaFree(d_fftArray));
    checkCudaErrors(cudaFree(d_tsArray));
    checkCudaErrors(cudaFree(d_gradArray));
    cufftDestroy(fftPlan);

}
double clsField::getSpectrum(double k0, double ksq)
{
    //float k=ksq/powf(k0,2);
    //return k*powf(1+k,-2.2);
    return (ksq/(pow(k0,4)))*exp(-1.0*ksq*pow((1.0/k0),2));

}

