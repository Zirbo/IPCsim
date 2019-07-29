#!/bin/sh

#SBATCH -J sIPC_39_42_48_50
#SBATCH -N 1
#SBATCH -p mem_0064
#SBATCH --tasks-per-node=16
#SBATCH --mail-user=silvano.ferrari@ymail.com
#SBATCH --time=0-18:00:00 

module purge
module load gcc/7.2 cmake/3.9.6
pushd ~/gitrepos/IPCsim/
    mkdir -p IPCsim
    cd IPCsim
    cmake ../src/2patches/
    make -j6
popd
cp ~/gitrepos/IPCsim/IPCsim/IPCsim .

for dir in symmetric/*
    do
    cp IPCsim $dir/IPCsim
    pushd $dir
    echo $dir
    touch info.txt
    OMP_NUM_THREADS=4 ./IPCsim oldstaged >> info.txt 2>&1 &
    popd
done

wait

module purge
