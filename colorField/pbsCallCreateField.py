import sys, os
#get inpute argument to pass to callCreateField.py
#read in params.txt
#get array index for easy queueing
#------
index=int(os.getenv('PBS_ARRAYID'))
INPUTFILE='params.txt'
f=open(INPUTFILE)
lines=f.readlines()
line = lines[index]
scale_noise=line.split(' ')[0]
visc=line.split(' ')[1]
sigma=line.split(' ')[2]
runCMD='python callCreateField.py '+str(scale_noise)+' '+str(visc)+' '+str(sigma)
os.system(runCMD)
