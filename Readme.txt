/*
This code is provided freely, however when using this code you are asked to cite these related paper: 
Puckett, J.G., Pokhrel, A.R. and Giannini, J.A. (2017) Collective gradient sensing in fish schools

This code is modified from original work by 
Berdahl, A., Torney, C.J., Ioannou, C.C., Faria, J. & Couzin, I.D. (2013) Emergent sensing of complex environments by mobile animal groups, Science






 *****Collective gradient sensing in fish schools**???
*/
1)To compile the code go the  code folder in command line and type in $makeall. Make sure you are using the right nvcc flag on make file while compiling the code.
2)The input to the compiled .zonalapp field is given as

 ./zonalApp 32 0.1 500 1 /media/data2/sim/final/Model_Couzin_Gradient/codeClean/test/00000 3 5.5 0.1 0.0 1
 where
	Number of particles=32
     	NoiseLevel=0.1
	Timestep=500
	Randseed=1
	The directory where you want the output file to be=/media/data2/sim...
	Radius of orientation=3 
	Radius of attraction=5.5
	Weight of the gradient=0.1
	Error in gradient sensing=0.0
	Visualisation (1 to enable Open CV)=1

3Parameter files are made using zMakeInputParameter.py and then fed into pbsCallZonal.py.
4)A big simulation is ran by submitting the params file through a submitPython.py that calls this pbsCallZonal.py