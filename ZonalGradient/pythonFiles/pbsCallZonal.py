import sys, os
import numpy as np
#get input argument to pass to callCreateField.py
#read in params.txt
#get array index for easy queueing
#./zonalApp 128 0.1 1000 1 /media/data2/sim/20160623a 2 6
#------
#index=int(os.getenv('PBS_ARRAYID'))
#index starts at 0
INPUTFILE=sys.argv[1]
print(INPUTFILE)
f=open(INPUTFILE)

def alreadyFinished(dirL):
	q = os.path.exists(dirL)
	
	T = 0
	if (q):
		fname = dirL+'/particles.dat'
		with open(fname, 'rb') as fh:
			lines = fh.readlines()
		    
	
	return q

while True:
	line=f.readline()
	if not line.rstrip(): break
	N=line.split(' ')[0]
	noiseScale=line.split(' ')[1]   
	numSteps=line.split(' ')[2]
	randSeed=line.split(' ')[3]
	runName=line.split(' ')[4]
	ozone=line.split(' ')[5]
	azone=line.split(' ')[6]
        weight=line.split(' ')[7]
	gradError=line.split(' ')[8]
	#check if already exist
	q = alreadyFinished(runName)
	runCMD='./zonalApp '+str(N)+' '+str(noiseScale)+' '+str(numSteps)+' '+str(randSeed)+' '+runName+' '+str(ozone)+' '+str(azone)+' '+str (weight)+' '+str (gradError)
	if (q==False):
		print "=======\n"
		print runCMD
		os.system(runCMD)
		
		#sys.exit(0)
	#runCMDpy='python postProcessMakeMovies.py  '+runName
	#print runCMDpy
	#os.system(runCMDpy)
