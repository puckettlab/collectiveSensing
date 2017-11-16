#define NOISE_TAIL -2.0
extern "C"
{
void generateRandomArray(int gridSize, float* d_rands, int* d_seeds, int numBlocks, int numThreads, int numRandBlocks, int numRandThreads, int numLoops);
void initializeArray(int gridSize, int numBlocks, int numThreads, float* d_fftArray, float* d_rands, float k0, float scaleFact);
void advanceArray(int gridSize, int numBlocks, int numThreads, float* d_fftArray, float* d_rands, float k0, float scaleFact, float dT, float revertRate);
void executeFFT(int gridSize, int numBlocks, int numThreads, float* d_fftArray, float* d_tsArray, cufftHandle fftPlan);
void convertGradients(int gridSize, int numBlocks, int numThreads, float* d_fftArray, float* d_gradArray, cufftHandle fftPlan);
//void gaussianRandsN(int gridSize, float* d_rands, int* d_seeds);
}

