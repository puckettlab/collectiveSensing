#ifndef PARTICLES_KERNEL_H
#define PARTICLES_KERNEL_H

#include "vector_types.h"
typedef unsigned int uint;


// simulation parameters
struct SimParams {
    float noise;

    float deltaTime;

    float attractionRadius;
    float orientationRadius;
    float repulsionRadius;
    float turnSpeed;
    float maxAngle;
    float maxSwimSpeed;
    float minSwimSpeed;
    int numThreads;
    int numBlocks;

    int noiseSize;
    float noiseScale;
 float weight;
 float gradError;	

};

__device__ float2 rotateVector(float2 newVector, float2 oldVector, float maxAngle);
__device__ float interpolateArray1D(float x_pos, float y_pos, int gridSize, float* in_Array);
__device__ float convertAngle(float2 inVector);
__device__ float calculateDarkness(float x, float y, float sx, float sy, int noiseSize, float* noiseField);
__device__ float2 addGradient(float2 newDirection,float x,float y,float sx,float sy,int noiseSize,float* noiseField,float gRandtheta);
__device__ float calculateControl(float x, float y, int noiseSize, float* controlField);
__global__ void uniformRandsD(float* result, int* seeds);
__global__ void gaussianRandsD(float* result, int* seeds);
#endif

