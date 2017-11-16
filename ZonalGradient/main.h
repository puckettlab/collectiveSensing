#ifndef __MAIN_H__
#define __MAIN_H__

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/core/core.hpp> 	//cv::Mat
#include <math.h>
#include <stdio.h>
#include <sys/stat.h>


float scale = 0.95;

void getRandColorArray(int numParticles, std::vector<cv::Scalar> &color)
{
	cv::RNG rng(12345); 				//random
	color.clear();
	for (int i=0;i<numParticles;i++)
	{
		float b = (float)rng.uniform(0,255) / 256;
		float g = (float)rng.uniform(0,255) / 256;
		float r = (float)rng.uniform(0,255) / 256;				
		color.push_back(cvScalar(b, g,r,0));
	}
}

void plotParticles(float *pos, int numParticles, int nx, std::vector<cv::Scalar> &color, cv::Mat &output)
{
	cv::Mat outputrgb;
// 	output = output*256;
	cv::cvtColor(output, outputrgb, CV_GRAY2RGB);
	
//	cv::Scalar color = cv::Scalar(200, 0, 0,25);
	for (int i=0; i<numParticles; i++) // run our simulation// 
	{
		float posx = (float)(nx)*pos[4*i];
		float posy = (float)(nx)*pos[4*(i)+1];	
         //std::cout<<pos[4*i]<<' '<<pos[4*(i)+1]<<'\n';	
		
		cv::circle(outputrgb, cv::Point(posx,posy), 3, color[i], -1,8,0 );		
	
	//		std::cout<<'x'<<' '<<posx<<' '<<posy<<'\n';		
	}
// 	cv::resize(outputrgb,outputrgb, cv::Size(),scale, scale, cv::INTER_NEAREST);	 		
	//std::cout<<"h11";
	cv::imshow("a", outputrgb);	
	 cv::waitKey(10);
		
}




#endif // __NOISEHELPER_H__
