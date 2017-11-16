#ifndef __NOISEHELPER_H__
#define __NOISEHELPER_H__

#include "noise.h" 			//

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp> 	//cv::Mat
#include <math.h>

class NoiseHelper
{
	private: 			//or public???
		float width;///
		
		float amp; //
		
		int xrange;
		int yrange;
		
		int xmin;
		int ymin;
		int xmax;
		int ymax;		
		
		float nmax;
		float navg;
		int   switch_time;
		float speed;
		float heading;
		float* n_array; 	//noise address
		bool debug;
		int offset;// location of the bob
		int N; 			//2048;//powf(2,ceil(log2(ny))); ///Number of grid points
	public:
		float x; //  xy of dark spot; 0..1
		float y; //	
		int nx;
		int ny;
 		cv::Mat noiseImg; 	//our noise image
 		float noiseScale; //0.25
 		float sigma;///
 		int Lborder;
 		cv::Math noiseImg2;
 		NoiseHelper(int Ni, int nxi, int nyi, float sigmai, float speedi, int Lborderi, float noiseScalei, float nmaxi, float navgi, float* n_arrayi, bool debugi)  ;
		void update(int t);
     
 
};


NoiseHelper::NoiseHelper(int Ni, int nxi, int nyi, float widthi, float speedi, int Lborderi, float noiseScalei, float nmaxi, float navgi, float * n_arrayi, bool debugi) 
{
	nx 		= nxi;
	ny 		= nyi;
	Lborder = Lborderi;	
	amp 	= 1.5; 					// no need to change
	noiseScale = noiseScalei; 	//
	xmin   	= Lborder;
	ymin   	= Lborder;
	xmax 	= nx-Lborder; 		/// decrease  from the right edge 
    ymax 	= ny-Lborder;		///decrease  from  bottom edge
	xrange 	= nx-xmin;
	yrange 	= ny-ymin;
	width 	= widthi; 			//1.0/(2.0*pow(sigma,1));/// 
	speed 	= speedi; 			//	speed = 0.00178; 		// in px, set????????????????????	
	switch_time=0;
	heading = 0;
	offset  = 20; 			// location of the bob
	debug 	= debugi;
    x 		= (rand()/(RAND_MAX+1.0)); //initial pos; 0-1 
    y 		= (rand()/(RAND_MAX+1.0));
	nmax 	= nmaxi;
	navg 	= navgi;	
	N 		= Ni;//2048;//powf(2,ceil(log2(ny))); ///Number of grid points
	n_array = n_arrayi;
	noiseImg= cv::Mat(nx,ny, CV_32F, 0.0);
   
}


void NoiseHelper::update(int t)
{	
	//update new heading
	float LL = 1.0;//*nx;
	if (t==switch_time)
	{
// 		float nextx =(xmin+offset) + ((xrange-offset) * (rand()/(RAND_MAX+1.0)));
// 		float nexty =(ymin+offset) + ((yrange-offset) * (rand()/(RAND_MAX+1.0)));
		float nextx = LL*(rand()/(RAND_MAX+1.0));
		float nexty = LL*(rand()/(RAND_MAX+1.0));
		//
		switch_time = t+1+(int)(sqrt(pow(x-nextx,2)+pow(y-nexty,2))/speed);
		heading 	= atan2(nexty-y, nextx-x);
	}
	x+=speed*cos(heading);
	y+=speed*sin(heading);
	if (x>LL)
		x-=LL;
	if (y>LL)
		y-=LL;
	if (x<0)
		x+=LL;
	if (y<0)
		y+=LL;
	// below is only for visual//debugging
	if (debug)
	{
		float xN = x*nx;
		float yN = y*ny;	
		float grey;
		for(int i = 0; i < nx;i++)
		{
			for(int j = 0; j < ny;j++)
			{
				if ((i<xmin)|(i>xmax)|(j<ymin)|(j>ymax))
					grey = 1;
				else
				{
					float x_d=(i - xN)/nx; 	//
					float y_d=(j - yN)/ny; 	//
			
					//only do type 3
					grey = 1.0-amp*exp(-powf(pow(x_d,2)+pow(y_d,2),0.5)/width) + noiseScale*((n_array[N*j+i]/nmax) );//-navg);
	// 				std::cout<<' '<<grey<<'\n';
				
				}
				noiseImg.at<float>(cv::Point(i,j))	= grey;		 				
			}
		}	
		noiseImg.setTo(0, noiseImg<0); 	//min is 0
		noiseImg.setTo(1, noiseImg>1);	//max is 1	
	}

}

#endif // __NOISEHELPER_H__
