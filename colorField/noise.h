
/*
% GNU Public Licence Copyright (c) Colin Torney
% Comments and questions to colin.j.torney@gmail.com

% This code is provided freely, however when using this code you are asked to cite our related paper: 
% Berdahl, A., Torney, C.J., Ioannou, C.C., Faria, J. & Couzin, I.D. (2013) Emergent sensing of complex environments by mobile animal groups, Science
*/

#ifndef _noise_h_included_
#define _noise_h_included_
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_complex.h>
#include <math.h>

class clsNoise
{
    public:
        clsNoise(int N_in, double dt_in, int seed,double visc_in);
        ~clsNoise();
        int generate_random_array(double *r_array);
        int advance_timestep();
        int prepare_save();
        double getSpectrum(double k0, double ksq);
        double* get_address();
        double* get_vel_address();
        double* get_x_vel();
        double* get_y_vel();
        double* get_x_second();
        double* get_y_second();
        double* get_xy_second();
    private:
	
        
        double k0;
        double l_scaling;
        double *sf_fourier;
        double *sf_normal;
        double *x_vel;
        double *y_vel;
        double *x_second;
        double *y_second;
        double *xy_second;
        const gsl_rng_type * T;
        gsl_rng * gsl_rand;
        double dt;
        int N;
	double visc;
};
#endif

