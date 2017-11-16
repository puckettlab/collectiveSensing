#ifndef _field_h_included_    
#define _field_h_included_
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_complex.h>
#include <math.h>
#include <cutil_inline.h>
#include <cufft.h>
#include <helper_cuda.h>
#include <noise.cuh>

class clsField 
{
    public:
        clsField(int N_in, double dt_in, int seed, double visc);
        ~clsField();
        void advanceTimestep();
        void prepareSave();
        double getSpectrum(double k0, double ksq);
        float* getDeviceGrads() { return d_gradArray; }
        float* getDeviceField() { return d_tsArray; }
        float* getHostGrads() { return h_gradArray; }
        float* getHostField() { return h_tsArray; }
        void setGradStrength(float in_gradStrength) {gradStrength=in_gradStrength;}


    private:
        double revertRate;
        double k0;
        double scaleFact;
        double dt;
        int gridSize;
        int numBlocks;
        int numThreads;
        int numRandBlocks;
        int numRandThreads;
        int numRandLoops;
        float gradStrength;
	int seed;//added
	double visc;//added

        
        float *d_fftArray;
        float *d_gradArray;
        float *d_tsArray;
        float *d_rands;
        int *d_seeds;
        float *h_gradArray;
        float *h_tsArray;
        cufftHandle fftPlan;
};
#endif

