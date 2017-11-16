#!/bin/bash
##use 1 machine, 1 processor
##not using this option, artifact
##PBS -l nodes=1:ppn=1 		
#PBS -l nodes=hydra3:ppn=4	
#PBS -t 0-1999
#PBS -N zonalSimulation
#PBS -S /bin/bash
#PBS -l walltime=50:00:00
#PBS -M 3135aawaz@gmail.com
#PBS -m e
#PBS -d /media/data1/aawaz/testing/parameterBasedCouzin
#PBS -j oe
#PBS -o /media/data1/aawaz/testing/parameterBasedCouzin/output/${PBS_JOBID}.out
PBS_O_WORKDIR='/media/data1/aawaz/testing/parameterBasedCouzin' ###uncomment when debugging
cd ${PBS_O_WORKDIR}
echo "Job ID is            ${PBS_JOBID}"
echo "Timestamp is         $(date +%F_%T)"
echo "Directory is         $(pwd)"
echo "Running on host      $(hostname)"
echo "Working directory is ${PBS_O_WORKDIR}"
echo "Starting the program"
echo "-----------"
##scale noise, visc sigma amp noise_type
#index=1 			##if NOT using array qsub 
#index=${PBS_ARRAYID}
#python pbsCallZonal.py params1.txt 
python pbsCallZonal.py /media/data1/aawaz/testing/parameterBasedCouzin/params-0.4/params-16/params${PBS_ARRAYID}.txt

