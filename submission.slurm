#!/bin/sh

#SBATCH -J 40_5c_r567
#SBATCH -N 1
#SBATCH -p mem_0064
#SBATCH --tasks-per-node=16
#SBATCH --mail-user=silvano.ferrari@ymail.com
#SBATCH --mail-type=BEGIN,END
#SBATCH --time=1-12:00:00 

module purge
module load gcc/7.2 cmake/3.9.6

IPCsim=~silfer/gitrepos/IPCsim/bld/IPCsim/IPCsim

runset=( rho55/40_5c rho65/40_5c rho75/40_5c )

for dir in ${runset[@]}; do
    echo $dir
    cp $IPCsim $dir/IPCsim
    pushd $dir
    OMP_NUM_THREADS=4 ./IPCsim oldstaged > stdouterr.txt 2>&1 &
    popd
done

wait

module purge
