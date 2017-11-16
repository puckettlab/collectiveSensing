/*

This code is provided freely, however when using this code you are asked to cite these related paper: 
Puckett, J.G., Pokhrel, A.R. and Giannini, J.A. (2017) Collective gradient sensing in fish schools

This code is modified from original work by 
Berdahl, A., Torney, C.J., Ioannou, C.C., Faria, J. & Couzin, I.D. (2013) Emergent sensing of complex environments by mobile animal groups, Science

GNU Public Licence Copyright (c) Colin Torney
Comments and questions to colin.j.torney@gmail.com

*/
#include "particleSystem.h"
#include "particleSystem.cuh"
#include "particles_kernel.cuh"

#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuda_runtime.h>
#include <helper_functions.h>
#include <helper_cuda.h>
#include <iostream>
#include <assert.h>
#include <math.h>
#include <memory.h>
#include <cstdio>
#include <cstdlib>
#include <algorithm>


#ifndef CUDART_PI_F
#define CUDART_PI_F         3.141592654f
#endif
using namespace std;
ParticleSystem::ParticleSystem(uint numParticles, uint numBlocks, uint numThreads, float deltaT, int noiseSize, float noiseScale, float a_zone, float o_zone,float weight,float gradError) 
{
    m_params.numBlocks=numBlocks;
    m_params.numThreads=numThreads;
    float BL 			= 0.0357; //couzin:: 0.03125 //0.02250
    m_params.repulsionRadius  	= 0.50 * BL;
    m_params.orientationRadius  = o_zone * BL;
    m_params.attractionRadius  = a_zone * BL;
    m_params.maxSwimSpeed 	= 5.0f * BL;
    m_params.minSwimSpeed 	= 1.0f * BL; //real min i 0.1
    m_params.noise=0.01f;
    m_params.noiseSize = noiseSize; 	// size of noise image
    m_params.noiseScale = noiseScale;   //image noise scale
    m_params.deltaTime = deltaT;
    m_params.turnSpeed = 1.75;
    m_params.maxAngle = m_params.deltaTime*m_params.turnSpeed;
    m_params.weight=weight;
    m_params.gradError=gradError;
 
    _initialize(numParticles);
	//
}

ParticleSystem::~ParticleSystem()
{
    //destructor
    delete [] m_hPos;
    delete [] m_hVel;
    delete [] m_hcdfVel;
    delete [] m_hcdfPos;
    delete [] m_hcdfCol;
    cudaFree(m_dVel);
    cudaFree(m_dPos);
    cudaFree(m_dControl);
    cudaFree(m_dRands);
    cudaFree(m_gRands);
    cudaFree(m_dSeeds);
}

void ParticleSystem::_initialize(int numParticles)
{

    m_numParticles = numParticles;
    timecount=0;

    // allocate host storage
    m_hPos = new float[m_numParticles*4];
    m_hVel = new float[m_numParticles*4];
    memset(m_hPos, 0, m_numParticles*4*sizeof(float));
    memset(m_hVel, 0, m_numParticles*4*sizeof(float));

    m_hcdfPos = new float[m_numParticles*3];
    m_hcdfVel = new float[m_numParticles*3];
    m_hcdfCol = new float[m_numParticles];
    memset(m_hcdfPos, 0, m_numParticles*3*sizeof(float));
    memset(m_hcdfVel, 0, m_numParticles*3*sizeof(float));

    // allocate GPU data
    unsigned int memSize = sizeof(float) * 4 * m_numParticles;

    cudaMalloc((void**)&m_dVel, memSize);
    cudaMalloc((void**)&m_dPos, memSize);

    timecount=0;
    noiseLevel = 4; ///?????

    cudaMalloc((void**)&m_dSeeds, m_params.numBlocks*sizeof(int));
    cudaMalloc((void**)&m_dRands, sizeof(float)*m_numParticles);
    cudaMalloc((void**)&m_gRands, sizeof(float)*m_numParticles);

    // *****************************************
    // setup background and control image
    // *****************************************
   
    m_hControl = new float[m_params.noiseSize*m_params.noiseSize];
    cudaMalloc((void**)&m_dControl, m_params.noiseSize*m_params.noiseSize*sizeof(float));

    const char* fileName = "controlFile.dat";
    FILE* file = fopen(fileName, "rb");


    for(int i = 0; i< m_params.noiseSize*m_params.noiseSize; i++)
    {
        float f;
        fread(&f, sizeof(float), 1, file);
        m_hControl[i] = f;
    }

    fclose(file);

    copyArrayToDevice(m_dControl, m_hControl, 0, m_params.noiseSize*m_params.noiseSize*sizeof(float));



    setParameters(&m_params);
}


void ParticleSystem::advanceTimestep(float xpatch, float ypatch)
{
    setParameters(&m_params);

      //debug

        loadArray();
    
    timecount++;
    if (timecount>2879) timecount=0; 		

    copyArrayToDevice(m_dField, m_hField, 0, 1920*1080*sizeof(float));
    
    



    generateRandomNumbers(m_dRands, m_dSeeds, m_params.numBlocks, m_params.numThreads);
    generateRandomNumbers(m_gRands, m_dSeeds, m_params.numBlocks, m_params.numThreads);

    // process social interactions and move individuals
    processInteractions(m_dPos, m_dVel, m_dField, m_dControl, m_dRands, m_gRands, m_numParticles, xpatch, ypatch, m_params.numThreads, m_params.numBlocks);
    
    
    return;

}

void ParticleSystem::loadArray()
{
    // copy position data from GPU to CPU
    copyArrayFromDevice(m_hPos, m_dPos, m_numParticles*4*sizeof(float));
    copyArrayFromDevice(m_hVel, m_dVel, m_numParticles*4*sizeof(float));
}

void ParticleSystem::prepareSave()
{
    unsigned int memSize = sizeof(float) * 3 * m_numParticles;

    float *m_d3DPos;
    float *m_d3DVel;
    cudaMalloc((void**)&m_d3DVel, memSize);
    cudaMalloc((void**)&m_d3DPos, memSize);

    convert4to3DArray(m_dPos, m_d3DPos, m_numParticles, m_params.numThreads, m_params.numBlocks );
    convert4to3DArray(m_dVel, m_d3DVel, m_numParticles, m_params.numThreads, m_params.numBlocks );
    // copy position data from GPU to CPU
    copyArrayFromDevice(m_hcdfPos, m_d3DPos, m_numParticles*3*sizeof(float));
    copyArrayFromDevice(m_hcdfVel, m_d3DVel, m_numParticles*3*sizeof(float));

    for(uint i=0; i < m_numParticles; i++) 
    {
        m_hcdfCol[i]=m_hcdfVel[i*3+2];
        m_hcdfPos[i*3+2]=0.0f;
        m_hcdfVel[i*3+2]=0.0f;
    }

    cudaFree(m_d3DPos);
    cudaFree(m_d3DVel);
}
void ParticleSystem::setArray()
{
    // copy position & velocity data from CPU to GPU
    copyArrayToDevice(m_dPos, m_hPos, 0, m_numParticles*4*sizeof(float));
    copyArrayToDevice(m_dVel, m_hVel, 0, m_numParticles*4*sizeof(float));
}

inline float frand()
{
    return rand() / (float) RAND_MAX;
}

void ParticleSystem::reset()
{
    memset(m_hPos, 0, m_numParticles*4*sizeof(float));
    memset(m_hVel, 0, m_numParticles*4*sizeof(float));

    int p = 0, v = 0;
    for(int j=0; j < m_params.numBlocks; j++) 
    {
        float side = min(0.25,sqrt(powf(m_params.repulsionRadius,2.0)*2.0*3.142*m_params.numThreads));
        float centerY = 0.25+0.5f*frand();
        float centerX = 0.25+0.5f*frand();
        for(int i=0; i < m_params.numThreads; i++) 
        {
            float point[3];
            point[0] = centerX + side*(frand()-0.5);//0.25+0.4*frand();
            point[1] = centerY + side*(frand()-0.5);//0.25+0.4*frand();
            point[2] = frand();
            m_hPos[p++] = (point[0] );
            m_hPos[p++] = (point[1] );
            m_hPos[p++] = 0.0;
            m_hPos[p++] = 0.0;
            m_hVel[v++] = cos(2*M_PI*point[2]);
            m_hVel[v++] = sin(2*M_PI*point[2]);
            m_hVel[v++] = 0.0f;
            m_hVel[v++] = 0.0f;

        }
    }

    setArray();
}

void ParticleSystem::resetResource()
{
    loadArray();
    for(uint i=0; i < m_numParticles; i++) 
    {
        m_hPos[i*4+2] = 0.0;
        m_hPos[i*4+3] = 0.0;
    }
    copyArrayToDevice(m_dPos, m_hPos, 0, m_numParticles*4*sizeof(float));
}

float ParticleSystem::calculateResource(float totalTime)
{
    //calculate the average resource experienced by a particle
    loadArray();
    float avResource=0.0f;

    for(uint i=0; i < m_numParticles; i++) 
        avResource+=m_hPos[i*4+2];

    return avResource/((float)m_numParticles*totalTime);

}

float ParticleSystem::calculateControl(float totalTime)
{
    //calculate the average resource experienced by a particle
    loadArray();
    float avControl=0.0f;

    for(uint i=0; i < m_numParticles; i++) 
        avControl+=m_hPos[i*4+3];

    return avControl/((float)m_numParticles*totalTime);
}
 
void ParticleSystem::saveDataPsi(ofstream &myfile)
{
	for(int i=0; i < m_numParticles; i++)
    	{
    	int ind = 4*i;
    	float psi = m_hPos[ind+2]/m_hPos[ind+3];
              if (!isnan(psi))
		myfile<<psi<<'\n';  //text file
    	}
}
void ParticleSystem::saveData(ofstream &myfile,int t)
{
	for(uint i=0; i < m_numParticles; i++)
    	{
		myfile<<i<<' '<<t<<' '<<m_hPos[i*4+0]<<' '<<m_hPos[i*4+1]<<' '<<m_hPos[i*4+2]<<' '<<m_hPos[i*4+3]<<' ';  //text file
		myfile<<' '<<m_hVel[i*4+0]<<' '<<m_hVel[i*4+1]<<' '<<m_hVel[i*4+2]<<' '<<m_hVel[i*4+3]<<'\n';  //text file
    	}
}
