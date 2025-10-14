#!/bin/bash

### parameters for start in queue. Not used for 

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -q  plasma@en067.binp.gcf
##$ -pe mpi 1
#$ -pe smp 256
#$ -j y
#$ -l h_rt=999:01:00

## start compile and install program to WorkDir
python3 configure.py $1

WorkDir=$(cat workdir.tmp)      # WORK DIRECTORY, setted in configure.py

if [$WorkDir == ""]
then
    exit
fi

np=$(cat proc.tmp)      # WORK DIRECTORY, setted in configure.py
name=$(cat name.tmp)      # WORK DIRECTORY, setted in configure.py

if [$WorkDir == ""]
then
    exit
fi
rm workdir.tmp
proc.tmp
rm name.tmp

cd $WorkDir

### run program
OMP_NUM_THREADS=$np OMP_PLACES=cores OMP_PROC_BIND=true numactl --interleave=all  .\/$name
#numactl --cpunodebind=1 --localalloc OMP_NUM_THREADS=64 OMP_PROC_BIND=close OMP_PLACES=cores   .\/$name

#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/master/soft/install/petsc/lib
#mpirun -np 8 ./beren3d.exe
