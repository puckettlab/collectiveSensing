import sys, os
import numpy as np
import os.path
    

dirL = '/media/data2/sim/final/Model_Couzin_GError/inputParams/params-0.25/params-32/'

for i in range(0,6720):
	runCmd= 'python /media/data1/aawaz/testing/ultimateCode/ultimate6.1/pbsCallZonal.py ' + dirL+'params'+str(i)+'.txt'
	paramFile = dirL+'params'+str(i)+'.txt'
	print runCmd
	read = np.genfromtxt(paramFile, delimiter=' ',dtype='str')
	dirS =  read[4]
	fileOut = dirS + '/psi.dat'
	print fileOut	
	if (not os.path.exists(fileOut) ): #does file exist??
		os.system(runCmd)
