
/*
% GNU Public Licence Copyright (c) Colin Torney
% Comments and questions to colin.j.torney@gmail.com

% This code is provided freely, however when using this code you are asked to cite our related paper: 
% Berdahl, A., Torney, C.J., Ioannou, C.C., Faria, J. & Couzin, I.D. (2013) Emergent sensing of complex environments by mobile animal groups, Science
*/

#include <iostream>
#include <string>
#include "noise.h"

/*
   --------------------------------------------------------------------------------
   Code for the coloured noise
   --------------------------------------------------------------------------------
 */

using namespace std;
clsNoise::clsNoise(int N_in, double dt_in, int seed, double visc_in)
{
    /*
       Initialise random number generator, set up arrays and define the initial velocity distribution
Parameters: Number of grid points, time step
     */
    gsl_rng_env_setup();
    T = gsl_rng_default;
    gsl_rand = gsl_rng_alloc (T);
    gsl_rng_set(gsl_rand,seed);

    N=N_in;
    dt=dt_in;
	visc=visc_in;  

    sf_fourier = new double [N*N*2];
    sf_normal = new double [N*N];
    x_vel = new double [N*N];
    y_vel = new double [N*N];
    x_second = new double [N*N];
    y_second = new double [N*N];
    xy_second = new double [N*N];
    // smaller visc means longer time correlations/no evolution of flow
     //visc=0.01;
   //visc=0.001;
///Temporal correlation of noise,dampens the fluctuation in the field	
    k0=24;///Dominant wave number of noise



    // Create the energy spectrum
    double ki,kj,ksq;
    double lambda_k=0.0;
    for (int i=1;i<N;i++)
    {
        ki=2*M_PI*(i-(N*0.5));
        for (int j=1;j<N;j++)
        {
            kj=2*M_PI*(j-(N*0.5));
            ksq=((ki*ki)+(kj*kj));
            lambda_k+=getSpectrum(k0,ksq);
        }
    }
   //double flow_speed=0.00016;
    double flow_speed=1.000;
    l_scaling=flow_speed/(lambda_k);

    // Create the initial velocity
    double *c_rand = new double [N*N*2];
    double alpha_ij, real_ij, imag_ij, mag_ij;
    generate_random_array(c_rand);
    for (int i=0;i<N;i++)
        for (int j=0;j<N;j++)
        {
            if ((i==0)||(j==0))
            {
                sf_fourier[(N*2*j)+(i*2)]=0.0;
                sf_fourier[(N*2*j)+(i*2)+1]=0.0;
            }
            else
            {
                ki=2*M_PI*(i-(N*0.5));
                kj=2*M_PI*(j-(N*0.5));
                ksq=((ki*ki)+(kj*kj));
                lambda_k=(l_scaling*visc*getSpectrum(k0,ksq));
                alpha_ij=visc*ksq;
                if (alpha_ij==0.0)
                    mag_ij=0.0;
                else
                    mag_ij=sqrt((lambda_k/(2.0*alpha_ij)));
                real_ij=mag_ij*c_rand[(j*N*2)+(i*2)];
                imag_ij=mag_ij*c_rand[(j*N*2)+(i*2)+1];
                sf_fourier[(N*2*j)+(i*2)]=real_ij;
                sf_fourier[(N*2*j)+(i*2)+1]=imag_ij;
            }
        }

    prepare_save();
    delete [] c_rand;
}

int clsNoise::generate_random_array(double *r_array)
{
    //random number generator
    int N2 = (int)(0.5*N);
    int i,j,ki,kj;
    int c_i,c_j,c_ki,c_kj;
    double real,imag;
    double stddev=sqrt(0.5);

    for (j=0;j<=N2;j++)
        for (i=0;i<=N2;i++)
        {
            if ((i==0)||(j==0))
            {
                r_array[(j*N*2)+(i*2)]=0.0;
                r_array[(j*N*2)+(i*2)+1]=0.0;
            }
            else
            {
                ki=(i-N2);
                kj=(j-N2);
                c_ki=-1*ki;
                c_kj=-1*kj;
                c_i=c_ki+N2;
                c_j=c_kj+N2;
                if ((ki==c_ki)&&(kj==c_kj))
                {
                    real = gsl_ran_gaussian (gsl_rand, 1);
                    r_array[(j*N*2)+(i*2)]=real;
                    r_array[(j*N*2)+(i*2)+1]=0.0;
                }
                else
                {
                    real = gsl_ran_gaussian (gsl_rand, stddev);
                    imag = gsl_ran_gaussian (gsl_rand, stddev);
                    r_array[(j*N*2)+(i*2)]=real;
                    r_array[(j*N*2)+(i*2)+1]=imag;
                    r_array[(c_j*N*2)+(c_i*2)]=real;
                    r_array[(c_j*N*2)+(c_i*2)+1]=-imag;
                }
            }
        }

    for (j=(N2+1);j<(N);j++)
        for (i=0;i<N2;i++)
        {
            if (i==0)
            {
                r_array[(j*N*2)+(i*2)]=0.0;
                r_array[(j*N*2)+(i*2)+1]=0.0;
            }
            else
            {
                ki=(i-N2);
                kj=(j-N2);
                c_ki=-1*ki;
                c_kj=-1*kj;
                c_i=c_ki+N2;
                c_j=c_kj+N2;

                real = gsl_ran_gaussian (gsl_rand, stddev);
                imag = gsl_ran_gaussian (gsl_rand, stddev);
                r_array[(j*N*2)+(i*2)]=real;
                r_array[(j*N*2)+(i*2)+1]=imag;
                r_array[(c_j*N*2)+(c_i*2)]=real;
                r_array[(c_j*N*2)+(c_i*2)+1]=-imag;
            }

        }

    return 0;
}
int clsNoise::advance_timestep()
{
    // this routine called from main code to advance to next flow realisation
    // using the mean reverting random process it advances the velocity in fourier space to the next time step
    double *nsf_fourier = new double [N*N*2];
    double lambda_k;
    double ki,kj,ksq;
    double i_dbl,j_dbl;

    double *c_rand = new double [N*N*2];
    double alpha_ij, real_ij, imag_ij, mag_ij;
    generate_random_array(c_rand);

    for (int i=0;i<N;i++)
        for (int j=0;j<N;j++)
        {
            if ((i==0)||(j==0))
            {
                nsf_fourier[(N*2*j)+(i*2)]=0.0;
                nsf_fourier[(N*2*j)+(i*2)+1]=0.0;
            }
            else
            {
                ki=2*M_PI*(i-(N*0.5));
                kj=2*M_PI*(j-(N*0.5));
                ksq=((ki*ki)+(kj*kj));
                lambda_k=l_scaling*visc*getSpectrum(k0,ksq);
                alpha_ij=visc*ksq;


                if (alpha_ij==0.0)
                    mag_ij=0.0;
                else
                    mag_ij=sqrt((lambda_k/(2.0*alpha_ij))*(1-exp(-2.0*alpha_ij*dt)));
                real_ij=mag_ij*c_rand[(j*N*2)+(i*2)];
                imag_ij=mag_ij*c_rand[(j*N*2)+(i*2)+1];
                nsf_fourier[(j*N*2)+(i*2)]=(sf_fourier[(j*N*2)+(i*2)]*exp(-1.0*alpha_ij*dt))+real_ij;
                nsf_fourier[(j*N*2)+(i*2)+1]=(sf_fourier[(j*N*2)+(i*2)+1]*exp(-1.0*alpha_ij*dt))+imag_ij;
            }
        }


    delete [] sf_fourier;
    delete [] c_rand;
    sf_fourier = nsf_fourier;

    prepare_save();
    return 0;
}
int clsNoise::prepare_save()
{
    // this routine converts the fourier space flow velocity to physical space and stores results in arrays
    gsl_complex_packed_array data;
    double *xtrans_array = new double [N*N*2];
    double *ytrans_array = new double [N*N*2];
    double *x2trans_array = new double [N*N*2];
    double *y2trans_array = new double [N*N*2];
    double *xytrans_array = new double [N*N*2];
    double ki,kj;
    int left,right,up,down;

    double *trans_array = new double [N*N*2];

    // remove for efficiency, next lines convert stream function to real values
    memcpy(trans_array,sf_fourier,N*N*2*sizeof(double));

    for (int i=0;i<N;i++)
    {   
        data=&trans_array[i*2];
        gsl_fft_complex_radix2_backward(data,N,N);
    }   
    for (int i=0;i<N;i++)
    {   
        data=&trans_array[i*2*N];
        gsl_fft_complex_radix2_backward(data,1,N);
    }   
    for (int i=0;i<N;i++)
        for (int j=0;j<N;j++)
            sf_normal[(j*N)+i]=((((i+j)%2)==0)?1.0:-1.0)*trans_array[(j*N*2)+(i*2)];
    delete [] trans_array;


    delete [] xtrans_array;
    delete [] ytrans_array;
    delete [] x2trans_array;
    delete [] y2trans_array;
    delete [] xytrans_array;

    return 0;
}
// external access to velocity fields
double* clsNoise::get_address()
{
    //returns real space amps
    return sf_normal;
}
double* clsNoise::get_x_vel()
{
    //returns x direction velocity in real space
    return x_vel;
}
double* clsNoise::get_y_vel()
{
    //returns y direction velocity in real space
    return y_vel;
}
//derivatives of velocities
double* clsNoise::get_x_second()
{
    return x_second;
}
double* clsNoise::get_y_second()
{
    return y_second;
}
double* clsNoise::get_xy_second()
{
    return xy_second;
}
clsNoise::~clsNoise()
{
    delete [] sf_fourier; delete [] sf_normal; 
    delete [] x_vel; delete [] y_vel; delete [] x_second; delete [] y_second; delete [] xy_second;
    gsl_rng_free (gsl_rand);
}

double clsNoise::getSpectrum(double k0, double ksq)
{
    float k=ksq/powf(k0,2);
    //return k*powf(1+k,-2.2);
    return (ksq/(pow(k0,4)))*exp(-1.0*ksq*pow((1.0/k0),2)); 
}

