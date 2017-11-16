#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os,  sys
#input into cpp:
#outputfolder randSeed  nx ny  scale_noise  tmax  dt sigma amp noiseType visc Lborder  speed
##
# ./my_program /media/data2/noise/noise1 1 960 540 0.25 100 0.025 50 1.5 3 0.002 50 1
##
randSeed=1
nx = 960
ny = 540
scale_noise = 0.25
dt   = 0.010#noise fluctuation speed
visc = 0.002
sigma = 100
amp = 1.5
noiseType=3
tmax=5000
speed=1
Lborder=25
fps=25#movie frame rate only
if len(sys.argv) < 3:
    print "You must set argument!!!"
    print "Usage: python callCreateField.py 1 0.1 1000 0.15 "
    print "Inputs:  randSeed scale_noise tmax dt"
    exit(0)
if len(sys.argv)>1:
    randSeed = sys.argv[1]
if len(sys.argv)>2:
    scale_noise  = sys.argv[2]
if len(sys.argv)>3:
    tmax = sys.argv[3]
if len(sys.argv)>3:
   dt = sys.argv[4]
#----
cmdPath="/media/data1/jgp/code/c/colorField/my_program"
ffmpegIncl="-framerate "+str(fps)
ffmpegOptsCompress="-qp 0 -threads 0 -r "+str(fps)+" -pix_fmt yuv420p -c:v libx264 -preset faster -crf 19 -threads 0 -f mp4 -loglevel 0 -an -y"
baseDir ="/media/data2/share/noise"
runname='noise='+str(scale_noise)+'-randSeed='+str(randSeed)
outputFolder=baseDir+'/'+runname+'/'
print '------' 
print outputFolder
print runname
print scale_noise
print visc
print sigma
print amp
print noiseType
print tmax
print '------' 
#------------------------
#initialize functions
#------------------------
def makeField():
	# ./my_program /media/data2/noise/noise1 1 960 540 0.25 100 0.025 50 1.5 3 0.002 50 1
	fcmd1=cmdPath+" "+outputFolder+" "+str(randSeed)+" "+str(nx)+" "+str(ny)+" "+str(scale_noise)+" "+str(tmax)+" "+str(dt)+" "+str(sigma)+" "+str(amp)+" "+str(noiseType)+" "+str(visc)+" "+str(Lborder)+" "+str(speed)
	print(fcmd1)
	os.system(fcmd1)
	ffmpegcmd1="ffmpeg "+ffmpegIncl+" " +"-i "+outputFolder+"/img/f%05d.png "+ffmpegOptsCompress+" "+outputFolder+"/"+runname+".mp4"
	print(ffmpegcmd1)
	os.system(ffmpegcmd1)
	chmodcmd="chmod 766 -R " +baseDir+'/'+runname
	sys.exit(0)

makeField()





# 1 float grey1 = 1.0-amp*exp(-(powf(x_d,2)+pow(y_d,2))*width);
# 2 float grey2 = 1.0-amp*exp(-powf(pow(x_d,2)+pow(y_d,2),0.5)*width) + scale_noise*powf(n_array[N*((int)floor(1.0*j))+(int)floor(1.0*i)]/nmax,3);
# 3 float grey3 = 1.0-amp*exp(-powf(pow(x_d,2)+pow(y_d,2),0.5)*width) + scale_noise*((n_array[N*j+i]/nmax)-nav);	
# 4 float grey4 = 0.5-scale_noise*powf(n_array[N*((int)floor(j))+(int)floor(i)]/(0.5*nmax),1);



'''
outputFolder=argv[1];
	randseed=atoi(argv[2]);
	nx=atoi(argv[3]);
 	ny=atoi(argv[4]);
	scale_noise=atof(argv[5]);
	tmax=atoi(argv[6]);
	dt=atof(argv[7]);
	sigma=atof(argv[8]);
	noiseType=atoi(argv[9]);
	visc=atof(argv[10]);
	Lborder=atoi(argv[11]);
	speed=atoi(argv[12]);
'''
