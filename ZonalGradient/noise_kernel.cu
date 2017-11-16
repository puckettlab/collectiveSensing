/*
GNU Public Licence Copyright (c) Colin Torney
Comments and questions to colin.j.torney@gmail.com

This code is provided freely, however when using this code you are asked to cite our related paper: 
Berdahl, A., Torney, C.J., Ioannou, C.C., Faria, J. & Couzin, I.D. (2013) Emergent sensing of complex environments by mobile animal groups, Science
*/


#include <cutil_inline.h>
#include <curand_kernel.h>
#include <cstdlib>
#include <cstdio>
#include <string.h>
//#include <particleSystem.cuh>
#include <noise.cuh>

#include <mt19937_ref.cu>


extern "C"
{
__device__ float genSpectrum(float k0, float ksq)
{
 //   float k=ksq/powf(k0,2);
 //   float k=ksq/powf(k0,2);
    return (ksq/(powf(k0,4)))*exp(-1.0*ksq*powf((1.0/k0),2));

    //return k*powf(1+k, NOISE_TAIL);
}


__global__ void gaussianRandsN(int gridSize, float* d_rands, int* d_seeds)
{
    int indx = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    int indy = __mul24(blockIdx.y,blockDim.y) + threadIdx.y;

    int indx1=(2*gridSize*indy)+(2*indx);
    int indx2=(2*gridSize*indy)+(2*indx)+1;
    curandState state;
    curand_init(d_seeds[blockIdx.x], indx1, 0, &state);
    d_rands[indx1] = curand_normal(&state);
    curand_init(d_seeds[blockIdx.y], indx2, 0, &state);
    d_rands[indx2] = curand_normal(&state);
}
__global__ void createRandArrayD(int gridSize, float* d_rands, float* d_gaussian)
{
    int i = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    int j = __mul24(blockIdx.y,blockDim.y) + threadIdx.y;

    int ind1=(2*gridSize*j)+(2*i);
    int ind2=(2*gridSize*j)+(2*i)+1;


    int gridSize2 = (int)(0.5*gridSize);
    int ki,kj;
    int c_i,c_j,c_ki,c_kj;
    double real,imag;
    double stddev=sqrt(0.5);

    if (j<=gridSize2)
    {
        if (i<=gridSize2)
        {
            if ((i==0)||(j==0))
            {
                d_rands[ind1]=0.0;
                d_rands[ind2]=0.0;
            }
            else
            {
                ki=(i-gridSize2);
                kj=(j-gridSize2);
                c_ki=-1*ki;
                c_kj=-1*kj;
                c_i=c_ki+gridSize2;
                c_j=c_kj+gridSize2;
                if ((ki==c_ki)&&(kj==c_kj))
                {
                    d_rands[ind1]=d_gaussian[ind1];
                    d_rands[ind2]=0.0;
                }
                else
                {
                    real = stddev*d_gaussian[ind1];
                    imag = stddev*d_gaussian[ind2];
                    d_rands[ind1]=real;
                    d_rands[ind2]=imag;
                    d_rands[(c_j*gridSize*2)+(c_i*2)]=real;
                    d_rands[(c_j*gridSize*2)+(c_i*2)+1]=-imag;
                }
            }
        }
    }
    else
    {
        if (i<gridSize2)
        {
            if (i==0)
            {
                d_rands[ind1]=0.0;
                d_rands[ind2]=0.0;
            }
            else
            {
                ki=(i-gridSize2);
                kj=(j-gridSize2);
                c_ki=-1*ki;
                c_kj=-1*kj;
                c_i=c_ki+gridSize2;
                c_j=c_kj+gridSize2;

                real = stddev*d_gaussian[ind1];
                imag = stddev*d_gaussian[ind2];
                d_rands[ind1]=real;
                d_rands[ind2]=imag;
                d_rands[(c_j*gridSize*2)+(c_i*2)]=real;
                d_rands[(c_j*gridSize*2)+(c_i*2)+1]=-imag;
            }

        }
    }

    return;
}
__global__ void initializeArrayD(int gridSize, float* d_fftArray, float* d_rands, float k0, float scaleFact)
{
    int i = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    int j = __mul24(blockIdx.y,blockDim.y) + threadIdx.y;

    int ind1=(2*gridSize*j)+(2*i);
    int ind2=(2*gridSize*j)+(2*i)+1;

    float ki, kj, ksq, lambda_k, mag_ij;

    if ((i==0)||(j==0))
    {
        d_fftArray[ind1]=0.0;
        d_fftArray[ind2]=0.0;
    }
    else
    {
        ki=2*M_PI*(i-(gridSize*0.5));
        kj=2*M_PI*(j-(gridSize*0.5));
        ksq=((ki*ki)+(kj*kj));
        lambda_k=scaleFact*genSpectrum(k0,ksq);

        if (ksq==0.0)
            mag_ij=0.0;
        else
            mag_ij=sqrt(lambda_k/(2.0*ksq));

        d_fftArray[ind1]=mag_ij*d_rands[ind1];
        d_fftArray[ind2]=mag_ij*d_rands[ind2];
    }
}
__global__ void advanceArrayD(int gridSize, float* d_fftArray, float* d_rands, float k0, float scaleFact, float deltaT, float revertRate)
{
    int i = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    int j = __mul24(blockIdx.y,blockDim.y) + threadIdx.y;

    int ind1=(2*gridSize*j)+(2*i);
    int ind2=(2*gridSize*j)+(2*i)+1;

    double lambda_k, ki, kj, ksq;
    double alpha_ij, real_ij, imag_ij, mag_ij;

    if ((i==0)||(j==0))
    {
        d_fftArray[ind1]=0.0;
        d_fftArray[ind2]=0.0;
    }
    else
    {
        ki=2*M_PI*(i-(gridSize*0.5));
        kj=2*M_PI*(j-(gridSize*0.5));
        ksq=((ki*ki)+(kj*kj));
        lambda_k=scaleFact*revertRate*genSpectrum(k0,ksq);

        alpha_ij=revertRate*ksq;

        if (alpha_ij==0.0)
            mag_ij=0.0;
        else
            mag_ij=sqrt((lambda_k/(2.0*alpha_ij))*(1-exp(-2.0*alpha_ij*deltaT)));

        real_ij=mag_ij*d_rands[ind1];
        imag_ij=mag_ij*d_rands[ind2];
        d_fftArray[ind1]*=exp(-1.0*alpha_ij*deltaT);
        d_fftArray[ind1]+=real_ij;
        d_fftArray[ind2]*=exp(-1.0*alpha_ij*deltaT);
        d_fftArray[ind2]+=imag_ij;
    }
}
__global__ void convertComplex2RealD(int gridSize, float* input, float* output)
{
    int i = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    int j = __mul24(blockIdx.y,blockDim.y) + threadIdx.y;

    output[__mul24(j,gridSize)+i]=((((i+j)%2)==0)?1.0:-1.0)*input[__mul24(__mul24(2,gridSize),j)+__mul24(2,i)];
}
__global__ void mergeVector(int gridSize, float* xInput, float* yInput, float* output)
{
    int i = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    int j = __mul24(blockIdx.y,blockDim.y) + threadIdx.y;

    int indx1 = __mul24(j,gridSize)+i;
    int indx2 = __mul24(__mul24(2,gridSize),j)+__mul24(2,i);

    output[indx2]=xInput[indx1];
    output[indx2+1]=yInput[indx1];
}
__global__ void convertGradD(int gridSize, float* input, float* xOutput, float* yOutput)
{
    int i = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    int j = __mul24(blockIdx.y,blockDim.y) + threadIdx.y;

    int ind1=(2*gridSize*j)+(2*i);
    int ind2=(2*gridSize*j)+(2*i)+1;
    float ki=2.0*M_PI*((float)i-((float)gridSize*0.5));
    float kj=2.0*M_PI*((float)j-((float)gridSize*0.5));

    xOutput[ind1]=-ki*input[ind2];
    xOutput[ind2]=ki*input[ind1];
    yOutput[ind1]=-kj*input[ind2];
    yOutput[ind2]=kj*input[ind1];
}

}


