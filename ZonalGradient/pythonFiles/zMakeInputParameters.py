import os,sys
import numpy as np
# 16 0.25 2000 1 ./ 3.0 5.5 0
noise = 0.25
gradError=np.array([0.0])
particles=np.array([32])
Tsteps = 10000; 
startSeed=0
endSeed=1
dirL = '/media/data1/aawaz/test/'
print dirL
azones     = np.array([5.5])
ozones     = np.array([3.0])
#weight=np.logspace(-2,2,num=49)  #multiple of 5+1 , 25 and 49 works 
weight=np.logspace(-2,3,num=31)  #multiple of 5+1 , 25 and 49 works 
weight= np.insert(weight, 0, 0.0)
print weight
Cnt = 0
Cnt0 = -1
#weight=([0.0]) ## unccomment this out to take data without weight
print particles
count=0
for p in particles:
	count=0
	for r in range(startSeed,endSeed):
		for err in gradError: 	
			Cnt=0
			for  w in weight:
				paramsFile="/media/data1/aawaz/testing/ultimateCode/ultimate6.1/params-%.2f/params-" %noise +str(p)+"/params"+str(count)+'.txt'
				print paramsFile
				f = open(paramsFile, 'w')
				line = str(p)+' '+str(noise)+' '+'10000'+' '+str(r)+' '+dirL+'seed'+str(r) +'/'+str (Cnt).zfill(6)+' '+'3.0'+' '+'5.5'+' ' + str(w) +' '+str(err)+'\n'
				f.write(line)
				f.close()
				Cnt=Cnt+1
 	        		count=count+1

