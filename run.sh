#!/bin/bash
CurrentDir=$(pwd) 
SourceDir="./srcBeren/Beren3D"
EigenPath=~/eigen-3.4.0/ #
#EigenPath=/home/eberendeev/eigen-3.4.0/
#EigenPath=~/sf_C/Work/eigen-3.3.2/
name=Beren3D.ex

WARNING='\033[37;1;41m'
BOLD='\033[1m'       #  ${BOLD}      # жирный шрифт (интенсивный цвет)
# Цвет текста:
BLACK='\033[0;30m'     #  ${BLACK}    # чёрный цвет знаков
RED='\033[0;31m'       #  ${RED}      # красный цвет знаков
GREEN='\033[0;32m'     #  ${GREEN}    # зелёный цвет знаков
YELLOW='\033[0;33m'     #  ${YELLOW}    # желтый цвет знаков
BLUE='\033[0;34m'       #  ${BLUE}      # синий цвет знаков
MAGENTA='\033[0;35m'     #  ${MAGENTA}    # фиолетовый цвет знаков
CYAN='\033[0;36m'       #  ${CYAN}      # цвет морской волны знаков
GRAY='\033[0;37m'       #  ${GRAY}      # серый цвет знаков

NORMAL='\033[0m'      #  ${NORMAL}    # все атрибуты по умолчанию

BuildDir=build

mkdir $BuildDir
cp $SourceDir/Makefile_cpu ./$BuildDir
cd ./$BuildDir

if [[ $1 == "clean" ]]; then
  echo "Clean object files..."
  make -f Makefile_cpu clean
  rm ../srcBeren/Constants/const.h
  rm ../srcBeren/Constants/defines.h
  echo "Object files has been removed!"
  cd $CurrentDir
  exit
fi

cd $CurrentDir

rm *.tmp
rm *.cfg
python set_params.py

echo "......Reading of auxilary variables......"
WorkDir=$(cat workdir.tmp)      # WORK DIRECTORY
proc=$(cat proc.tmp)      # Number of procs
queue=$(cat queue.tmp)      # queue type
clu=$(cat cluster.tmp)      # cluster

if [[ -z $WorkDir ]]
then
  exit
fi


   mv *.h ./srcBeren/Constants
   mv *.par ./srcBeren/Constants

echo "......Compile......"
cd ./$BuildDir

if [[ $clu = nks1p ]]
then
module purge
module load intel/2017.4.196 parallel/mpi.intel.broadwell/2017.4.196 compilers/intel/2017.4.196
fi 

#make -f Makefile_cpu clean
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

if [[ $queue == home ]]
then

OMP_NUM_THREADS=4 OMP_PLACES=cores OMP_PROC_BIND=true  ./$name

else

#./Clusters/submit."$queue"."$clu".sh  $proc $name
  if [[ $clu == nks1p ]]
  then
    sbatch ./start.sh
  fi
#qsub submit.sh
fi
