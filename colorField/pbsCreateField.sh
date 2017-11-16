#!/bin/bash
##to change index in vim::  23gg, then edit
#PBS -l nodes=1:ppn=1
##PBS -l nodes=hydra1:ppn=1
#PBS -t 1-3
#PBS -N colorField
#PBS -S /bin/bash
#PBS -l walltime=100:00:00
#PBS -M 3135aawaz@gmail.com
#PBS -m abe
#PBS -d /media/data1/jgp/code/c/colorField
#PBS -j oe
#PBS -o ${PBS_JOBID}.out
PBS_O_WORKDIR='/media/data1/jgp/code/c/colorField'
cd ${PBS_O_WORKDIR}
echo "Job ID is            ${PBS_JOBID}"
echo "Timestamp is         $(date +%F_%T)"
echo "Directory is         $(pwd)"
echo "Running on host      $(hostname)"
echo "Working directory is ${PBS_O_WORKDIR}"
echo "Starting the program"
echo "-----------"
##scale noise, visc sigma amp noise_type
index=${PBS_ARRAYID}
params="${PBS_O_WORKDIR}/params.txt"
paramsCMD="$index"p
input=$(sed -n $paramsCMD  $params)
randSeed=`echo $input |  cut -d ' ' -f 1`
noise=`echo $input |  cut -d ' ' -f 2`
tmax=`echo $input |  cut -d ' ' -f 3`
dt=`echo $input |  cut -d ' ' -f 4`
python callCreateFieldSingle.py $randSeed $noise $tmax $dt

