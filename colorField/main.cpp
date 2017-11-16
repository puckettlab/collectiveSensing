
/*
% GNU Public Licence Copyright (c) Colin Torney
% Comments and questions to colin.j.torney@gmail.com

% This code is provided freely, however when using this code you are asked to cite our related paper: 
% Berdahl, A., Torney, C.J., Ioannou, C.C., Faria, J. & Couzin, I.D. (2013) Emergent sensing of complex environments by mobile animal groups, Science


%
%
%JGP update to code to write images a binary files which can be read my PsychoToolbox / Matlab at 30Hz
%
%This code is provided freely, however when using this code you are asked to cite these related paper: 
Puckett, J.G., Pokhrel, A.R. and Giannini, J.A. (2017) Collective gradient sensing in fish schools
*/

#include <pngwriter.h>
#include "writeImageBinary.cpp"
#include <math.h>
#include <iostream>
#include <sstream>
#include "noise.h"
#include "noise.cpp"
#include <fstream>
#include <string>
#include <cstdlib>
// basic file operations
using namespace std;

int main(int argc, char** argv) 
{
/*
input arguments; how to use
./myprogram randint nx ny scale_noise 

order of argument
1 randint
2 nx
3 ny
4 scale_noise
5 length of movie tmax

*/

/// check if enough arguments given; if (argc>4)
int randseed=1;
int seedoffset=321;
int nx;
int ny;
float scale_noise;
int tmax;
float dt;
float sigma;
float amp;
int noiseType;
double visc;
float speed;
int Lborder;
string outputFolder;


if (argc<13){///speed of the bob
	cout<<"Not enough argments\n";
	cout<<"Usage: ./myprogram outputFolder randseed nx ny scale_noise tmax dt sigma amp noiseType visc Lborder speed\n";
	cout<<"Usage: ./myprogram /media/data2/share/noise/noise1 1 960 540 0.25 100 0.025 50 1.5 3 0.002 50 1\n";
	cout<<"randseed: 1     int\n";
	cout<<"nx:       1940  int\n";	
	return(0);
	}
else{
	outputFolder=argv[1];
	randseed=atoi(argv[2]);
	nx=atoi(argv[3]);
 	ny=atoi(argv[4]);
	scale_noise=atof(argv[5]);
	tmax=atoi(argv[6]);
	dt=atof(argv[7]);
	sigma=atof(argv[8]);
	amp=atof(argv[9]);
	noiseType=atoi(argv[10]);
	visc=atof(argv[11]);
	Lborder=atoi(argv[12]);
	speed=atoi(argv[13]);
 }

    
    cout<<"outputFolder: "<<outputFolder<<"\n";
    cout<<"randseed:     "<<randseed<<"\n";
    cout<<"nx:           "<<nx<<"\n";
    cout<<"ny:           "<<ny<<"\n";
    cout<<"scale_noise:  "<<scale_noise<<"\n";
    cout<<"tmax:         "<<tmax<<"\n";
    cout<<"dt:           "<<dt<<"\n";
    cout<<"sigma:        "<<sigma<<"\n";
    cout<<"amp:         "<<amp<<"\n";

    cout<<"noiseType:    "<<noiseType<<"\n";
    cout<<"visc:         "<<visc<<"\n";
    cout<<"Lborder:      "<<Lborder<<"\n";
    cout<<"speed:        "<<speed<<"\n";
    //create output folders
    string mkdirCMD;
    mkdirCMD = "mkdir -p "+outputFolder;
    cout<<mkdirCMD<<'\n';
    int dir_err = system(mkdirCMD.c_str());
    if (-1 == dir_err)
    {
        printf("Error creating directory!n");
        exit(1);
    }
    mkdirCMD = "mkdir -p "+outputFolder+"/output/";
    dir_err = system(mkdirCMD.c_str());
    mkdirCMD = "mkdir -p "+outputFolder+"/img/";
    dir_err = system(mkdirCMD.c_str());
   

    //
    srand(randseed+seedoffset);

    float xmin=Lborder;///decrease from the left edge
   
    float xmax=nx-Lborder; /// decrease  from the right edge 
    
    float ymin=Lborder;///decrease  from  top edge
   
    float ymax=ny-Lborder;///decrease  from  bottom edge


    int N=1024;//2048;//powf(2,ceil(log2(ny))); ///Number of grid points

    clsNoise noise(N,dt,randseed+seedoffset,visc); ///Number of grid points,time step,seed
    
    //float scale_noise=0.1;//

    double* n_array = noise.get_address();
	

    //float sigma=50; // size of gaussian blob  ///controlled by width ,Decay length of the patch 
   // float amp=1.5;  // strong of signal  ///
    float x;
    float y;

    char fontname[35];
    sprintf(fontname, "/home/aawaz/Downloads/create_field/fonts/FreeMonoBold.ttf");
    

   // float omega=0.01;
   // float speed=1;	//  controls the speed of the bob


   //int tmax=2880;//length of the movie ////

    float width = 1.0/(2.0*pow(sigma,1));///
    
    int xrange = xmax-xmin;
    int yrange = ymax-ymin;
	std::cout<<xrange<<'\n';
    int switch_time=0;
    float heading=0;

    int offset0=100;// location of the bob initial, starting
    int offset=75;// location of the bob
    x = (xmax-xmin)/2 + ((offset0) * (rand()/(RAND_MAX+1.0)));
    y = (ymax-ymin)/2 + ((offset0) * (rand()/(RAND_MAX+1.0)));
	 float nmax=0.0;
    for (int i=0;i<N;i++)
        if (nmax<n_array[i])
            nmax=n_array[i];
	
	//
    char maskname[11];

	//open stream for writing to the file, x,y
  	ofstream fout;
	string fname;
	string fNotesname;
	fNotesname = outputFolder + "/Notes.txt";
	cout<<fNotesname<<'\n';
	fout.open(fNotesname.c_str());
    	



	 for (int t=0;t<tmax;t++)
    	{
		std::cout<<x<<'\n';
		//update new heading
		if (t==switch_time)
		{
		    float nextx =(xmin+offset) + ((xrange-offset) * (rand()/(RAND_MAX+1.0)));
		    float nexty =(ymin+offset) + ((yrange-offset) * (rand()/(RAND_MAX+1.0)));

			     // cout<<x<<":"<<y<<endl;
		    //           cout<<nextx<<":"<<nexty<<":";
		    switch_time = t+1+(int)(sqrt(pow(x-nextx,2)+pow(y-nexty,2))/speed);
		    //         cout<<switch_time<<endl;
		    heading = atan2(nexty-y, nextx-x);
		    //           cout<<heading<<endl<<";;;"<<endl;
		}
		char fnameB[100];
		char fnameZ[500];
		char fnameI[500];
		sprintf(fnameB, "/img/f%05d.png", t);
		strcpy(fnameI,outputFolder.c_str());
		strcat(fnameI,fnameB);
		pngwriter png(nx,ny,1.0,fnameI);
		
 		sprintf(fnameB, "/output/f%05d.png", t);
		strcpy(fnameZ,outputFolder.c_str());
		strcat(fnameZ,fnameB);
		ofstream bin(fnameZ, ios::binary);
		bin.write(reinterpret_cast<char*>(&nx), sizeof(int));
		bin.write(reinterpret_cast<char*>(&ny), sizeof(int));

		x+=speed*cos(heading);
		y+=speed*sin(heading);
		// write x,y data to file	
	//	fout<<x<<" "<<y<<endl;
			
		float nav=0.0;
		int nc = 0;
		for(int i = xmin; i < xmax;i++)
		    for(int j = ymin; j < ymax;j++)
		    {
			nav+=n_array[N*j+i];
			nc++;
		    }
		nav/=(nmax*(float)(nc));
		
		 //     y=ymin+(ymax-ymin)*(1.0+offset+cos(omega*t) + offset*cos(omega*t*0.5) )/(2.0*(1.0+offset));
		float grey;
		unsigned char greyi;
		for(int i = 0; i < nx;i++){
		    for(int j = 0; j < ny;j++)
		    {
			if ((i<xmin)|(i>xmax)|(j<ymin)|(j>ymax))
				grey = 1;
			else{
				float x_d=i-x;
				float y_d=j-y;

				
				if (noiseType==1){
				   grey = 1.0-amp*exp(-(powf(x_d,2)+pow(y_d,2))*width);
				}
				 if (noiseType==2){	
				    grey = 1.0-amp*exp(-powf(pow(x_d,2)+pow(y_d,2),0.5)*width) + scale_noise*powf(n_array[N*((int)floor(1.0*j))+(int)floor(1.0*i)]/nmax,3);
				}		
				 if (noiseType==3){
				     grey = 1.0-amp*exp(-powf(pow(x_d,2)+pow(y_d,2),0.5)*width) + scale_noise*((n_array[N*j+i]/nmax)-nav);
				}
				if (noiseType==4){
				      grey = 0.5-scale_noise*powf(n_array[N*((int)floor(j))+(int)floor(i)]/(0.5*nmax),1);
				      
				}
				png.plot(i,j, grey, grey,grey);
			}
			grey = 255*grey;
			if (grey>255)
				grey = 255;
			if (grey<0)
				grey = 0;
			greyi =(unsigned char)(grey);
			bin.write(reinterpret_cast<const char *>(&greyi),sizeof(greyi));	
		   }
		}	
	        bin.close();
		png.close();
	//mask.close();
        	noise.advance_timestep();
	
 
    }
    
    fout.close();

    return 0;
}

