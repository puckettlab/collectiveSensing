#ifndef __PARTICLESYSTEM_H__
#define __PARTICLESYSTEM_H__

#define DEBUG_GRID 0
#define DO_TIMING 0

#include "particles_kernel.cuh"
#include "vector_functions.h"
#include <fstream>
#include <iostream>
#include <opencv2/highgui/highgui.hpp> 			//also included in noisehelper 
#include <opencv2/imgproc/imgproc.hpp> 
#include <opencv2/core/core.hpp>

// Particle system class
class ParticleSystem
{
public:
    ParticleSystem(uint numParticles, uint numBlocks, uint numThreads, float deltaT, int noiseSize, float noiseScale, float a_zone, float o_zone,float weight,float gradError);
    ~ParticleSystem();

    enum ParticleArray
    {
        POSITION,
        VELOCITY,
    };

    void advanceTimestep(float xpatch, float ypatch);
    void reset();
    void resetResource();

    void prepareSave();
    void loadArray();
    void setArray();

    int getNumParticles() const { return m_numParticles; }

    ///void setFieldArray(float* in_fieldArray) {m_dFieldArray=in_fieldArray;}


    float getParticleRadius() { return 2.5f*m_params.repulsionRadius; }
    float* getPositions() {return m_hPos;}
    float* getVelocity() {return m_hcdfVel;}
    float* getColour() {return m_hcdfCol;}
    float* getFieldArray() {return m_hField;}
    void setNoiseField(float* in_noise_field) {m_dField = in_noise_field;}
   

    float calculateResource(float totalTime);
    float calculateControl(float totalTime);

    float* get_dPos() {return m_dPos;}

    uint m_numParticles;

    float* m_hPos;              // particle positions
    float* m_hVel;              // particle velocities
    float* m_hcdfPos;              // particle positions
    float* m_hcdfVel;              // particle velocities
    float* m_hcdfCol;              // particle velocities

    
   void saveData(std::ofstream &myfile,int timeStep); //added
   void saveDataPsi(std::ofstream &myfile); //added 
  

protected: // methods
    ParticleSystem() {}

    void _initialize(int numParticles);
    void _finalize();
   

protected: // data

    // CPU data
    int timecount;
    int randseed;
    float* m_hMask;
    int noiseLevel;

    float noiseScale;
    // GPU data
    float* m_dPos;
    float* m_dVel;


    // order parameters for K-model
    float* m_dField;
    float* m_hField;
    float* m_dControl;
    float* m_hControl;
    float* m_hxGrad;
    float* m_hyGrad; 
    float* m_dxGrad;
    float* m_dyGrad; 

    float *m_dRands;
    float *m_gRands;
    int *m_dSeeds;


    // params
    SimParams m_params;
};

#endif // __PARTICLESYSTEM_H__

