#! /bin/bash

#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -q test_mpi@en001.binp.gcf
#$ -pe mpi 1
#$ -j y
#$ -l h_rt=360:01:00

#!/bin/bash
CurrentDir=$(pwd)
SourceDir="./srcBeren/simulation"
EigenPath=/home/berendeev/bpi/Progs/eigen-3.4.0/ #
#EigenPath=~/sf_C/Work/eigen-3.3.2/
name=Beren3D.ex

BuildDir=build

mkdir $BuildDir
cp $SourceDir/Makefile_cpu ./$BuildDir

cd $CurrentDir

rm *.tmp
rm *.cfg
python set_params.py



echo "......Reading of auxilary variables......"
WorkDir=$(cat workdir.tmp)      # WORK DIRECTORY
proc=$(cat proc.tmp)      # Number of procs

if [[ -z $WorkDir ]]
then
  exit
fi


mv *.h ./srcBeren/constants
mv *.par ./srcBeren/constants

echo "......Compile......"

cd ./$BuildDir


# New gcc compiler
export PATH="/ceph/groups/plasma/bpi/Progs/GCC-10.2/bin:$PATH"
export LD_LIBRARY_PATH="/ceph/groups/plasma/bpi/Progs/GCC-10.2/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/ceph/groups/plasma/bpi/Progs/GCC-10.2/lib64:$LD_LIBRARY_PATH"

# New OpenMPI
export PATH="/ceph/groups/plasma/bpi/Progs/OpenMPI/bin:$PATH"
export LD_LIBRARY_PATH="/ceph/groups/plasma/bpi/Progs/OpenMPI/lib:$LD_LIBRARY_PATH"




make -f Makefile_cpu EIGEN_PATH=$EigenPath
echo "************ End compile **************"
if [[ $(ls | grep $name) != $name ]]
then
echo -en "************ ${WARNING}COMPILATION FAILED${NORMAL} **************\n"
fi
cd $CurrentDir


echo "......Copy files to work directory......"   
   rm *.tmp
   rm -r $WorkDir
   mkdir $WorkDir
   cp -r PlotScripts ./$WorkDir
   cp -r srcBeren ./$WorkDir
   cp -r Scripts ./$WorkDir
   cp -r Clusters ./$WorkDir
   cp set_params.py ./$WorkDir
   mv *.cfg ./$WorkDir
   cp *.sh ./$WorkDir
   mv $BuildDir/$name ./$WorkDir
   cd ./$WorkDir


OMP_NUM_THREADS=128 OMP_PLACES=cores OMP_PROC_BIND=true ./$name

