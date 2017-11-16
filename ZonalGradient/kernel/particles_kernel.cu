/*

This code is provided freely, however when using this code you are asked to cite these related paper: 
Puckett, J.G., Pokhrel, A.R. and Giannini, J.A. (2017) Collective gradient sensing in fish schools

This code is modified from original work by 
Berdahl, A., Torney, C.J., Ioannou, C.C., Faria, J. & Couzin, I.D. (2013) Emergent sensing of complex environments by mobile animal groups, Science

GNU Public Licence Copyright (c) Colin Torney
Comments and questions to colin.j.torney@gmail.com

*/


/* 
 * Device code.
 */

#ifndef _PARTICLES_KERNEL_H_
#define _PARTICLES_KERNEL_H_

#include <stdio.h>
#include <math.h>
#include <cutil_math.h>
#include "math_constants.h"
#include "particles_kernel.cuh"
#include "device_functions.h"

#include "cuda_runtime.h"
//#include "device_launch_parameters.h"
#include <curand_kernel.h>

#define EPS 1e-16
using namespace std;

// textures for particle position and velocity
texture<float, 1, cudaReadModeElementType> dRandTex;

__device__ inline float signof(float a) { return (a<0 ? -1 : 1); }

__constant__ SimParams params;

__global__ void interactD(float4* newVel, float4* newPos, int numParticles, float* darkArray, float* dControl, float xpatch, float ypatch, float* gRands)
{
    int index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;
    int start_p = __mul24(blockIdx.x,blockDim.x);
    int stop_p = __mul24(blockIdx.x,blockDim.x)+blockDim.x;
    if (index<numParticles)
    {
        float4 pos = newPos[index];
        float4 vel = newVel[index];

        float2 newDirection=make_float2(0.0f);

        float2 repDirection = make_float2(0.0f);
        float2 attDirection = make_float2(0.0f);

        float maxDist = 0.5;;
        float worldSize = 1.0f;

        bool repulsionMode=false;

        // examine only individuals in group (block)
        for (int j=start_p;j<stop_p;j++)
        {
            // check not colliding with self
            if (j != index) 
            {             
                float4 posj = newPos[j];
                float4 velj = newVel[j];
                float2 dist = make_float2(posj.x-pos.x,posj.y-pos.y);

                if (dist.x>maxDist) dist.x-=worldSize;
                if (dist.x<-maxDist) dist.x+=worldSize;
                if (dist.y>maxDist) dist.y-=worldSize;
                if (dist.y<-maxDist) dist.y+=worldSize;

                float magDist=length(dist);
                if (magDist<params.repulsionRadius)
                {
                    repulsionMode=true;
                    if (magDist>EPS)
                        repDirection-=(dist/magDist);
                }
                if (!repulsionMode)
                {
                    if (magDist<(params.orientationRadius))
                    {
                        float2 vel2 = make_float2(velj.x,velj.y);
                        attDirection+= vel2/length(vel2);
                    }
                    if (magDist<(params.attractionRadius))
                    {
                        attDirection+=(dist/magDist);
                    }
                }

            }
        }

        // sync within the block so velocities are consistent across time frame
        __syncthreads();
        if (repulsionMode)
        {
            float magDir=length(repDirection);
            //check we're not too close to avoid nan
            if (magDir>EPS)
            {
                repDirection/=magDir;
                newDirection=rotateVector(repDirection, make_float2(vel.x, vel.y), params.maxAngle);
            }
            else
                newDirection = make_float2(vel.x, vel.y);
        }
        else
        {
            
            float magDir=length(attDirection);
            if (magDir>EPS) 
            {
                attDirection/=magDir;
            
                // write new velocity back to original unsorted location
		newDirection=attDirection;
                //newDirection=rotateVector(attDirection, make_float2(vel.x, vel.y), params.maxAngle);  //
		
            }
            else{
                newDirection = make_float2(vel.x, vel.y);
                }
	
	float th 		= gRands[index];
        float gRandtheta 	=  params.gradError*sqrt(params.deltaTime)*th;
	     float2 dirGradient=addGradient(newDirection,pos.x,pos.y,xpatch,ypatch,params.noiseSize,darkArray, gRandtheta);
       	    //newDirection = dirGradient;  
            newDirection=rotateVector(dirGradient, make_float2(vel.x, vel.y), params.maxAngle);
        }
        // current light level is stored in vel.z (0 = dark, 1=light)
	 float curSpeed = params.minSwimSpeed + vel.z*(params.maxSwimSpeed-params.minSwimSpeed);
       
        // new position = ol d position + velocity * deltaTime
        float newposx = pos.x+(params.deltaTime * curSpeed * newDirection.x );
        float newposy = pos.y+(params.deltaTime * curSpeed * newDirection.y );
        if (newposx>1.0f)
            newposx-=1.0f;

        if (newposx<0.0f)  
            newposx+=1.0f;

        if (newposy>1.0f)
            newposy-=1.0f;

        if (newposy<0.0)
            newposy+=1.0f;
	
	float darkness = calculateDarkness(newposx, newposy, xpatch, ypatch, params.noiseSize, darkArray);
	float control  = calculateControl(newposx, newposy, params.noiseSize, dControl);
      
        newVel[index] = make_float4(newDirection.x, newDirection.y, darkness, vel.w);
        newPos[index] = make_float4(newposx, newposy, pos.z+(1.0f-darkness), pos.w+(1.0f-control));
    }
}


__device__ float calculateDarkness(float x, float y, float sx, float sy, int noiseSize, float* noiseField)
{
    float amp = 1.5;
    float width =  0.275;
    float maxNoise = 0.1145;
    float noiseScale = params.noiseScale;
    // calculate distance to patch center
    float xdiff = x - sx;
    float ydiff = y - sy;

    if (xdiff > 0.5)
        xdiff = 1.0f - xdiff;
    if (xdiff < -0.5)
        xdiff = 1.0f + xdiff;
    if (ydiff > 0.5)
        ydiff = 1.0f - ydiff;
    if (ydiff < -0.5)
        ydiff = 1.0f + ydiff;

    // interpolate to find noise value at position
    float x_grid = y*(float)noiseSize;
    float y_grid = x*(float)noiseSize;

    int left=(int)floor(x_grid);
    int right=(int)ceil(x_grid);
    int down=(int)floor(y_grid);
    int up=(int)ceil(y_grid);
    
    if (right==noiseSize) right=0;
    if (left==noiseSize) left=0;
    if (down==noiseSize) down=0;
    if (up==noiseSize) up=0;
    
    int downLeft=down*noiseSize+left;
    int downRight=down*noiseSize+right;
    int upLeft=up*noiseSize+left;
    int upRight=up*noiseSize+right;
    
    float downLeftVal =  noiseField[downLeft];
    float downRightVal =  noiseField[downRight];
    float upLeftVal =  noiseField[upLeft];
    float upRightVal = noiseField[upRight];
    
    float x_w=x_grid-left;
    float y_w=y_grid-down;
    
    float noiseLevel = ((1-x_w)*(1-y_w)*downLeftVal)+((x_w)*(1-y_w)*downRightVal)+((1-x_w)*(y_w)*upLeftVal)+((x_w)*(y_w)*upRightVal);
	//printf("%4.4f\n",noiseLevel);
    return max(min(1.0f,1.0-amp*exp(-powf(pow(xdiff,2)+pow(ydiff,2),0.5)/width) + noiseScale*(noiseLevel/maxNoise)),0.0f);
}    

__device__ float calculateControl(float x, float y, int noiseSize, float* controlField)
{
    // interpolate to find noise value at position
    float x_grid = y*(float)noiseSize;
    float y_grid = x*(float)noiseSize;

    int left=(int)floor(x_grid);
    int right=(int)ceil(x_grid);
    int down=(int)floor(y_grid);
    int up=(int)ceil(y_grid);
    
    if (right==noiseSize) right=0;
    if (left==noiseSize) left=0;
    if (down==noiseSize) down=0;
    if (up==noiseSize) up=0;
    
    int downLeft=down*noiseSize+left;
    int downRight=down*noiseSize+right;
    int upLeft=up*noiseSize+left;
    int upRight=up*noiseSize+right;
    
    float downLeftVal =  controlField[downLeft];
    float downRightVal =  controlField[downRight];
    float upLeftVal =  controlField[upLeft];
    float upRightVal = controlField[upRight];
    
    float x_w=x_grid-left;
    float y_w=y_grid-down;
    
    return min(1.0f,((1-x_w)*(1-y_w)*downLeftVal)+((x_w)*(1-y_w)*downRightVal)+((1-x_w)*(y_w)*upLeftVal)+((x_w)*(y_w)*upRightVal));
    
} 
//calculate gradDir fucntion
__device__ float fixPBC(float x){
	if (x>0.5)
		x = 1.0f - x;
	if (x<-0.5)
		x = 1.0f + x;
	return x;
}
  //add gradeint sensing
__device__ float2 addGradient(float2 newDirec, float x, float y, float sx, float sy,int noiseSize, float* noiseField,  float gRandtheta)
{
    float noiseScale = params.noiseScale;
    float amp = 1.5;
    float width =  0.275; //0.275
    float maxNoise = 0.1145;
    
    float xdiff = x - sx;
    float ydiff = y - sy;
    xdiff = fixPBC(xdiff);
    ydiff = fixPBC(ydiff);

    float2 returnVector=make_float2(0.0f);

    // interpolate to find noise value at position
    float x_grid = y*(float)noiseSize;
    float y_grid = x*(float)noiseSize;

    int left =(int)floor(x_grid);
    int right=(int)ceil(x_grid);
    int down =(int)floor(y_grid);
    int up   =(int)ceil(y_grid);


    float sleft =(float)down/noiseSize - sx;
    float sright=(float)up/noiseSize - sx;
    float sdown =(float)left/noiseSize - sy;
    float sup   =(float)right/noiseSize - sy;

	sleft  = fixPBC(sleft);
	sright = fixPBC(sright);
	sdown  = fixPBC(sdown);
	sup   = fixPBC(sup);

    if (right==noiseSize) right=0;
    if (left==noiseSize) left=0;
    if (down==noiseSize) down=0;
    if (up==noiseSize) up=0;
//
    int downLeft=down*noiseSize+left;
    int downRight=down*noiseSize+right;
    int upLeft=up*noiseSize+left;
    int upRight=up*noiseSize+right;

	//
    float downLeftVal =  noiseField[downLeft];
    float downRightVal =  noiseField[downRight];
    float upLeftVal =  noiseField[upLeft];
    float upRightVal = noiseField[upRight];

	//printf("Values %2.2f %2.2f %2.2f %2.2f \n",downLeftVal,downRightVal,upLeftVal,upRightVal);
 	float downLeftI = max(min(1.0f,1.0-amp*exp(-powf(pow(sleft,2)+pow(sdown,2),0.5)/width) + noiseScale*(downLeftVal/maxNoise)),0.0f);
     	float downRightI = max(min(1.0f,1.0-amp*exp(-powf(pow(sright,2)+pow(sdown,2),0.5)/width) + noiseScale*(downRightVal/maxNoise)),0.0f);
    	float upLeftI = max(min(1.0f,1.0-amp*exp(-powf(pow(sleft,2)+pow(sup,2),0.5)/width) + noiseScale*(upLeftVal/maxNoise)),0.0f);
  	float upRightI = max(min(1.0f,1.0-amp*exp(-powf(pow(sright,2)+pow(sup,2),0.5)/width) +   noiseScale*(upRightVal/maxNoise)),0.0f);
	
   	float dx,dy;
	
	    
             dx = -(upRightI+downRightI-upLeftI-downLeftI);   
           
               dy = -(upRightI+upLeftI-downRightI-downLeftI); 
             
	float2 gradient = make_float2(dx,dy);
       
    	float gradientD=length(gradient);
	
	float gradientMin = 0.0001; 			
           
            if (gradientD>gradientMin) 
            {
			float2 gradientTmp= make_float2(gradient.x, gradient.y);
			gradient /= length(gradient);     //maxGradient;
        		gradient.x 	= cos(gRandtheta)*gradientTmp.x - sin(gRandtheta)*gradientTmp.y;
        		gradient.y 	= sin(gRandtheta)*gradientTmp.x + cos(gRandtheta)*gradientTmp.y;
	
                      
		     newDirec /= length(newDirec);
                   
		     newDirec += params.weight * gradient; 		
		     newDirec /= length(newDirec);
            }
            else{
                newDirec = newDirec;
		}
     return newDirec;
}   

__global__ void convert4to3DArrayD(float4* input, float3* output, int numParticles)
{

    int index = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

    if (index<numParticles)
    {
        output[index].x = input[index].x;
        output[index].y = input[index].y;
        output[index].z = input[index].z;
    }

}

__device__ float2 rotateVector(float2 newVector, float2 oldVector, float maxAngle)
{
    float newAngle=acos(newVector.x*oldVector.x+newVector.y*oldVector.y);
    float2 returnVector=make_float2(0.0f);
    if (newAngle>maxAngle)
    {
        float sign=signof(newVector.y*oldVector.x-newVector.x*oldVector.y);
        returnVector = make_float2(cos(maxAngle)*oldVector.x-sign*sin(maxAngle)*oldVector.y,sign*sin(maxAngle)*oldVector.x+cos(maxAngle)*oldVector.y);
    }
    else
        returnVector = newVector;

    return returnVector/length(returnVector);
}    



__global__ void uniformRandsD(float* result, int* seeds)
{
    //mt19937si(seeds[blockIdx.x]);
    //result[blockIdx.x * blockDim.x + threadIdx.x] = (mt19937sl()/(4294967296.0));
	//

	int i = blockDim.x * blockIdx.x + threadIdx.x;
    	curandState state;
    	curand_init(seeds[blockIdx.x], i, 0, &state);
    	result[i] = curand_uniform(&state);
	//printf("\t%3.3f\n",result[i]);
}

__global__ void gaussianRandsD(float* result, int* seeds)
{
	int i = blockDim.x * blockIdx.x + threadIdx.x;
    	curandState state;
    	curand_init(seeds[blockIdx.x], i, 0, &state);
	result[i] = curand_normal(&state);
	result[i+1] = curand_normal(&state);
}

__global__ void addRotationalNoiseD(float4* newVel, int numParticles, float* dRands)
{
    //add noise term to direction vector

    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i<numParticles)
    {
	float th = dRands[i];
        float rndn =  params.noise*sqrt(params.deltaTime)*th;
        float2 newDirection=make_float2(newVel[i].x, newVel[i].y);
        newVel[i].x = cos(rndn)*newDirection.x - sin(rndn)*newDirection.y;
        newVel[i].y = sin(rndn)*newDirection.x + cos(rndn)*newDirection.y;
    }
}


/*

__global__ void convertGaussianD(float* d_rands)
{
    int indx = __mul24(blockIdx.x,blockDim.x) + threadIdx.x;

    if ((indx%2)==0)
    {
        int indx1=indx;
        int indx2=indx+1;

        float x1=d_rands[indx1];
        float x2=d_rands[indx2];

        if (x1<1e-6) x1=1e-6;
        // box muller
        d_rands[indx1] = sqrtf(-2.0*logf(x1))*cosf(2*M_PI*x2);
        d_rands[indx2] = sqrtf(-2.0*logf(x1))*sinf(2*M_PI*x2);
    }
}
*/
#endif

