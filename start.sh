#!/bin/bash

### parameters for start in queue. Not used for 

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -q  plasma@en001.binp.gcf
##$ -pe mpi 1
#$ -j y
#$ -l h_rt=999:01:00

## start compile and install program to WorkDir
python3 configure.py $1

WorkDir=$(cat workdir.tmp)      # WORK DIRECTORY, setted in configure.py

if [$WorkDir == ""]
then
    exit
fi
rm workdir.tmp
cd $WorkDir

### run program
OMP_NUM_THREADS=8 OMP_PLACES=cores OMP_PROC_BIND=true  ./beren3d.exe
