### parameters for start in queue. Not used for loca run
#!/usr/bin/env bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -q  plasma@en067.binp.gcf
##$ -pe mpi 1
#$ -pe smp 256
#$ -j y
#$ -l h_rt=999:01:00
set -euo pipefail

EIGEN_PATH=~/soft/eigen-3.4.0/
AMGCL_PATH=~/soft/amgcl/
#EIGEN_PATH="/home/berendeev/bpi/Progs/eigen-3.4.0/"
#AMGCL_PATH="/home/berendeev/bpi/Progs/amgcl/"

np=8

export EIGEN_PATH
export AMGCL_PATH

python3 build.py "$@"

WORKDIR_FILE="workdir.tmp"

if [ ! -f "$WORKDIR_FILE" ]; then
    echo "Ошибка: файл $WORKDIR_FILE не найден" >&2
    exit 1
fi

# Читаем первую строку (или всё содержимое) и удаляем пробельные символы
RUN_DIR=$(cat "$WORKDIR_FILE" | tr -d '[:space:]')
rm -f workdir.tmp 

if [ -z "$RUN_DIR" ]; then
    echo "Ошибка: файл $WORKDIR_FILE пуст" >&2
    rm -f "$WORKDIR_FILE"
    exit 1
fi


if [ ! -x "$RUN_DIR/beren3d" ]; then
    echo "Ошибка: исполняемый файл не найден в $RUN_DIR" >&2
    exit 1
fi

cd $RUN_DIR
OMP_NUM_THREADS=$np OMP_PLACES=cores OMP_PROC_BIND=true numactl --interleave=all  ./beren3d