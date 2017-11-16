#include <particles_kernel.cuh>
extern "C"
{


void convert4to3DArray(float* input, float* output, int numParticles, int numThreads, int numBlocks);
void copyArrayFromDevice(void* host, const void* device, int size);
void copyArrayToDevice(void* device, const void* host, int offset, int size);

void setParameters(SimParams *hostParams);

void processInteractions(float* curPos, float*  curVel, 
        float*  fieldArray,
        float*  dControl,
        float*  dRands,
	 float*  gRands,
        uint    numParticles,
        float   xpatch,
        float   ypatch,
        int     numThreads, 
        int     numBlocks
        );
 
void addNoise(float* curVel, float* dRands, uint numParticles, int numThreads, int numBlocks);
void generateRandomNumbers(float* d_rands, int* d_seeds, int numBlocks, int numThreads);//, bool
}

